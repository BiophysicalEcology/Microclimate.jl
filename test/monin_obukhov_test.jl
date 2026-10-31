using Microclimate
using Microclimate: calc_Obukhov_length, monin_obukhov_iteration, ObukhovLength, InverseObukhovLength, stability_parameter,
    calc_ψ_m, calc_ψ_h, calc_Φ_h, calc_ρ_cp, atmospheric_surface_profile!, allocate_profile, dry_air_properties
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
    ρ_cp = calc_ρ_cp(mean_temp)
    kinematic_viscosity = dry_air_properties(mean_temp, atmospheric_pressure).kinematic_viscosity
    return (reference_temp, surface_temp, 2.0u"m/s", z0, z, ρcpTκg, bl.karman_constant, ΔT, ρ_cp, kinematic_viscosity)
end
stability_keywords = (; γ=bl.dyer_constant, stable_beta=bl.stable_beta, turbulent_prandtl_number=bl.turbulent_prandtl_number)

obukhov(dT) = calc_Obukhov_length(iteration_arguments(dT)...; stability_keywords..., initial_obukhov_length=-0.3u"m")
inverse_obukhov(dT) = monin_obukhov_iteration(InverseObukhovLength(), iteration_arguments(dT)...; stability_keywords...,
    initial_value=inv(-0.3u"m"))

# Stability parameter at the reference height for an Obukhov length
ζ(obukhov_length) = stability_parameter(ObukhovLength(), z, obukhov_length)

@testset "calc_Obukhov_length: sign and continuity across the stable/unstable transition" begin
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
        @test isapprox(ζ(direct.obukhov_length), stability_parameter(InverseObukhovLength(), z, inverse.value); atol = 0.05)
        @test isapprox(direct.friction_velocity, inverse.friction_velocity; rtol = 0.01)
        @test isapprox(direct.convective_heat_flux, inverse.convective_heat_flux; rtol = 0.01, atol = 0.1u"W/m^2")
    end
    # Neutral: finite, and the stability parameter converges to zero
    @test isapprox(stability_parameter(InverseObukhovLength(), z, inverse_obukhov(0.0u"K").value), 0.0; atol = 0.01)
end

@testset "calc_ψ_m/calc_ψ_h/calc_Φ_h: type-stable and zero-allocation, both regimes" begin
    for L in (-5.0u"m", 5.0u"m", Inf * u"m")
        @inferred calc_ψ_m(ζ(L), bl.dyer_constant, bl.stable_beta)
        @inferred calc_ψ_h(ζ(L), bl.dyer_constant, bl.stable_beta, bl.turbulent_prandtl_number)
        @inferred calc_Φ_h(ζ(L), bl.dyer_constant, bl.stable_Φ_h_coefficient, bl.min_stable_Φ_h, bl.max_stable_Φ_h)
        @test allocations(calc_ψ_m, ζ(L), bl.dyer_constant, bl.stable_beta) == 0
        @test allocations(calc_ψ_h, ζ(L), bl.dyer_constant, bl.stable_beta, bl.turbulent_prandtl_number) == 0
        @test allocations(calc_Φ_h, ζ(L), bl.dyer_constant, bl.stable_Φ_h_coefficient, bl.min_stable_Φ_h, bl.max_stable_Φ_h) == 0
    end
end

@testset "calc_Φ_h: stable branch increases with height, unstable branch decreases" begin
    # Φ_h = 1 at neutral; stable Φ_h > 1 (suppressed mixing), unstable Φ_h < 1 (enhanced).
    @test calc_Φ_h(ζ(Inf * u"m"), bl.dyer_constant, bl.stable_Φ_h_coefficient, bl.min_stable_Φ_h, bl.max_stable_Φ_h) == 1.0
    @test calc_Φ_h(ζ(5.0u"m"), bl.dyer_constant, bl.stable_Φ_h_coefficient, bl.min_stable_Φ_h, bl.max_stable_Φ_h) > 1.0
    @test calc_Φ_h(ζ(-5.0u"m"), bl.dyer_constant, bl.stable_Φ_h_coefficient, bl.min_stable_Φ_h, bl.max_stable_Φ_h) < 1.0
end

@testset "calc_Φ_h: stable branch stays bounded even for very small L (calm, strongly-stable)" begin
    # Regression test: a linear (unclamped) stable Φ_h let very small L (common
    # on calm clear nights, real full-year data) inflate RaupachLTheoryAirProfile's
    # T_L by orders of magnitude -- R's own dphih (sensibleheat.R) is a
    # saturating form for exactly this reason, matched here.
    for L in (10.0u"m", 1.0u"m", 0.1u"m", 0.01u"m", 0.0001u"m")
        Φ_h = calc_Φ_h(ζ(L), bl.dyer_constant, bl.stable_Φ_h_coefficient, bl.min_stable_Φ_h, bl.max_stable_Φ_h)
        @test bl.min_stable_Φ_h <= Φ_h <= bl.max_stable_Φ_h
    end
end

@testset "atmospheric_surface_profile!: no crash, sensible output across regimes" begin
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
