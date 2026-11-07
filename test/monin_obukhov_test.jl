using Microclimate
using Microclimate: calc_obukhov_length, monin_obukhov_iteration, ObukhovLength, InverseObukhovLength, stability_parameter_at_height,
    momentum_stability_correction, heat_stability_correction, heat_gradient_stability_factor, calc_air_volumetric_heat_capacity, atmospheric_surface_profile!, allocate_profile, dry_air_properties
using Unitful
using Test

include("helpers.jl")

bl = MoninObukhov()
z0 = 0.01u"m"
z = 10.0u"m"
ρcpTκg = u"J*s^2/m^4"(6.003e-8u"cal*minute^2/cm^4")
atmospheric_pressure = 101325.0u"Pa"

function iteration_arguments(dT)  # dT = reference_temp - surface_temp
    reference_temp = 293.15u"K"
    surface_temp = reference_temp - dT
    ΔT = reference_temp - surface_temp
    mean_temp = (reference_temp + surface_temp) / 2
    air_volumetric_heat_capacity = calc_air_volumetric_heat_capacity(mean_temp)
    kinematic_viscosity = dry_air_properties(mean_temp, atmospheric_pressure).kinematic_viscosity
    return (reference_temp, surface_temp, 2.0u"m/s", z0, z, ρcpTκg, bl.karman_constant, ΔT, air_volumetric_heat_capacity, kinematic_viscosity)
end
stability_keywords = (; bl.stability_functions)

obukhov(dT) = calc_obukhov_length(iteration_arguments(dT)...; stability_keywords..., initial_obukhov_length=-0.3u"m")
inverse_obukhov(dT) = monin_obukhov_iteration(InverseObukhovLength(), iteration_arguments(dT)...; stability_keywords...,
    initial_value=inv(-0.3u"m"))

# Stability parameter at the reference height for an Obukhov length
stability_parameter_for(obukhov_length) = stability_parameter_at_height(ObukhovLength(), z, obukhov_length)

@testset "calc_obukhov_length: sign and continuity across the stable/unstable transition" begin
    unstable = obukhov(-2.0u"K")  # surface warmer than reference
    neutral = obukhov(0.0u"K")
    stable = obukhov(2.0u"K")     # surface cooler than reference
    @test ustrip(u"m", unstable.obukhov_length) < 0
    @test isinf(ustrip(u"m", neutral.obukhov_length))
    @test ustrip(u"m", stable.obukhov_length) > 0
    @test unstable.convective_heat_flux < 0.0u"W/m^2"
    @test stable.convective_heat_flux > 0.0u"W/m^2"

    # L should vary smoothly (no jump/oscillation) as dT sweeps the crossover.
    dTs = -3.0u"K":0.25u"K":3.0u"K"
    Ls = ustrip.(u"m", [obukhov(dT).obukhov_length for dT in dTs])
    @test all(isfinite, Ls[Ls .!= Inf .&& Ls .!= -Inf])  # no NaN anywhere
    # sign(L) matches sign(dT) throughout (dT>0 => surface cooler => stable => L>0)
    @test all(zip(dTs, Ls)) do (dT, L)
        ustrip(u"K", dT) == 0.0 || sign(L) == sign(ustrip(u"K", dT))
    end
end

@testset "inverse Obukhov length iteration matches the Obukhov length iteration" begin
    for dT in (-3.0u"K", -1.0u"K", -0.25u"K", 0.25u"K", 1.0u"K", 3.0u"K")
        direct = obukhov(dT)
        inverse = inverse_obukhov(dT)
        @test isapprox(stability_parameter_for(direct.obukhov_length), stability_parameter_at_height(InverseObukhovLength(), z, inverse.value); atol = 0.05)
        @test isapprox(direct.friction_velocity, inverse.friction_velocity; rtol = 0.01)
        @test isapprox(direct.convective_heat_flux, inverse.convective_heat_flux; rtol = 0.01, atol = 0.1u"W/m^2")
    end
    # Neutral: finite, and the stability parameter converges to zero
    @test isapprox(stability_parameter_at_height(InverseObukhovLength(), z, inverse_obukhov(0.0u"K").value), 0.0; atol = 0.01)
end

@testset "$(nameof(typeof(f))) momentum_stability_correction/heat_stability_correction/heat_gradient_stability_factor: type-stable and zero-allocation, both regimes" for f in
        (BusingerDyer(), ChengBrutsaert())
    for L in (-5.0u"m", 5.0u"m", Inf * u"m")
        @inferred momentum_stability_correction(f, stability_parameter_for(L))
        @inferred heat_stability_correction(f, stability_parameter_for(L))
        @inferred heat_gradient_stability_factor(f, stability_parameter_for(L))
        @test allocations(momentum_stability_correction, f, stability_parameter_for(L)) == 0
        @test allocations(heat_stability_correction, f, stability_parameter_for(L)) == 0
        @test allocations(heat_gradient_stability_factor, f, stability_parameter_for(L)) == 0
    end
    # Φ_h = 1 at neutral; stable Φ_h > 1 (suppressed mixing), unstable Φ_h < 1 (enhanced).
    @test heat_gradient_stability_factor(f, stability_parameter_for(Inf * u"m")) == 1.0
    @test heat_gradient_stability_factor(f, stability_parameter_for(5.0u"m")) > 1.0
    @test heat_gradient_stability_factor(f, stability_parameter_for(-5.0u"m")) < 1.0
    # Stays bounded for very small L (calm, strongly stable), keeping
    # RaupachLTheoryAirProfile's T_L bounded
    for L in (10.0u"m", 1.0u"m", 0.1u"m", 0.01u"m", 0.0001u"m")
        @test 1.0 <= heat_gradient_stability_factor(f, stability_parameter_for(L)) <= 10.0
    end
end

@testset "ChengBrutsaert: smooth and strictly decreasing through the stable range" begin
    f = ChengBrutsaert()
    stability_parameters = range(0.1, 50.0; length=2001)
    for ψ in (momentum_stability_correction, heat_stability_correction)
        values = ψ.(Ref(f), stability_parameters)
        @test all(<(0), diff(values))
        # No kinks: the slope changes gradually everywhere (BusingerDyer's clamp jumps it by ~6)
        slopes = diff(values) ./ step(stability_parameters)
        @test maximum(abs, diff(slopes)) < 0.5
    end
    # Φ_h is the heat gradient of ψ_h, 1 - ζ dψ_h/dζ
    for stability_parameter in (0.01, 0.5, 2.0, 20.0)
        h = 1e-6
        slope = (heat_stability_correction(f, stability_parameter + h) - heat_stability_correction(f, stability_parameter - h)) / 2h
        @test heat_gradient_stability_factor(f, stability_parameter) ≈ 1 - stability_parameter * slope rtol = 1e-6
    end
end

@testset "$(nameof(typeof(bl.stability_functions))) atmospheric_surface_profile!: no crash, sensible output across regimes" for bl in
        (MoninObukhov(), MoninObukhov(; stability_functions=BusingerDyer()))
    site = Site(; latitude=0.0u"°", longitude=0.0u"°", elevation=0.0u"m", slope=0.0u"°", aspect=0.0u"°",
        horizon_angles=fill(0.0u"°", 24), sky_view_fraction=1.0, albedo=0.15, roughness_height=0.01u"m",
        atmospheric_pressure=101325.0u"Pa")
    heights = [0.5, 1.0, 2.0, 5.0, 10.0]u"m"
    buffers = allocate_profile(bl, heights)
    for (reference_temperature, zenith_angle) in (
        (300.0u"K", 30.0u"°"),   # warm reference, daytime
        (280.0u"K", 30.0u"°"),   # cold reference, daytime
        (280.0u"K", 170.0u"°"),  # night, surface colder than reference (typical stable night)
        (300.0u"K", 170.0u"°"),  # night, surface *warmer* than reference (atypical, e.g. warm advection)
    )
        environment_instant = (; atmospheric_pressure=101325.0u"Pa", reference_temperature,
            reference_wind_speed=2.0u"m/s", reference_humidity=0.6, zenith_angle)
        result = atmospheric_surface_profile!(bl, buffers; site, environment_instant, surface_temperature=290.0u"K")
        @test all(isfinite, ustrip.(u"K", result.air_temperature))
        @test all(isfinite, ustrip.(u"m/s", result.wind_speed))
        @test all(t -> 0.0 <= t <= 1.0, result.relative_humidity)
        @test all(>=(0.0u"m/s"), result.wind_speed)
    end
end
