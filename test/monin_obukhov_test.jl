using Microclimate
using Microclimate: calc_Obukhov_length, monin_obukhov_iteration, ObukhovLength, InverseObukhovLength, stability_parameter,
    calc_ψ_m, calc_ψ_h, calc_Φ_h, calc_ρ_cp, atmospheric_surface_profile!, allocate_profile, dry_air_properties
using Unitful
using Test

include("helpers.jl")

bl = MoninObukhov()
z0 = 0.01u"m"
z = 10.0u"m"
atmospheric_pressure = 101325.0u"Pa"

function iteration_arguments(dT)  # dT = reference_temp - surface_temp
    reference_temp = 293.15u"K"
    surface_temp = reference_temp - dT
    ΔT = reference_temp - surface_temp
    mean_temp = (reference_temp + surface_temp) / 2
    ρ_cp = calc_ρ_cp(mean_temp)
    kinematic_viscosity = dry_air_properties(mean_temp, atmospheric_pressure).kinematic_viscosity
    return (reference_temp, surface_temp, 2.0u"m/s", z0, z, ΔT, ρ_cp, kinematic_viscosity)
end

obukhov(dT) = calc_Obukhov_length(bl, iteration_arguments(dT)...; initial_obukhov_length=-0.3u"m")
inverse_obukhov(dT) = monin_obukhov_iteration(bl, InverseObukhovLength(), iteration_arguments(dT)...;
    initial_value=inv(-0.3u"m"))

# Stability parameter at the reference height for an Obukhov length
ζ(obukhov_length) = stability_parameter(ObukhovLength(), z, obukhov_length)

@testset "calc_Obukhov_length: sign and continuity across the stable/unstable transition" begin
    unstable = obukhov(-2.0u"K")  # surface warmer than reference
    neutral = obukhov(0.0u"K")
    stable = obukhov(2.0u"K")     # surface cooler than reference
    @test ustrip(u"m", unstable.obukhov_length) < 0
    @test isapprox(ζ(neutral.obukhov_length), 0.0; atol = 0.01)
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

@testset "bulk Stanton number stays positive in very stable air" begin
    for stability in (0.5, 5.0, 20.0, 100.0)
        heat_correction = calc_ψ_h(bl, stability)
        @test Microclimate._bulk_stanton(bl, log(z / z0), heat_correction, stability, z0, z) > 0
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

stability_pairs = (
    MoninObukhov(),
    MoninObukhov(; stable_stability=BeljaarsHoltslag()),
    MoninObukhov(; stable_stability=BusingerLinear()),
    MoninObukhov(; unstable_stability=NicheMapRUnstable(), stable_stability=NicheMapRStable()),
)
pair_name(model) = "$(nameof(typeof(model.unstable_stability)))/$(nameof(typeof(model.stable_stability)))"

@testset "$(pair_name(model)) calc_ψ_m/calc_ψ_h/calc_Φ_h: type-stable, zero-allocation and continuous at neutral" for model in stability_pairs
    for L in (-5.0u"m", 5.0u"m", Inf * u"m")
        @inferred calc_ψ_m(model, ζ(L))
        @inferred calc_ψ_h(model, ζ(L))
        @inferred calc_Φ_h(model, ζ(L))
        @test allocations(calc_ψ_m, model, ζ(L)) == 0
        @test allocations(calc_ψ_h, model, ζ(L)) == 0
        @test allocations(calc_Φ_h, model, ζ(L)) == 0
    end
    # Φ_h = 1 at neutral; stable Φ_h > 1 (suppressed mixing), unstable Φ_h < 1 (enhanced).
    @test calc_Φ_h(model, ζ(Inf * u"m")) == 1.0
    @test calc_Φ_h(model, ζ(5.0u"m")) > 1.0
    @test calc_Φ_h(model, ζ(-5.0u"m")) < 1.0
    # Both sides meet at neutral
    for ψ in (calc_ψ_m, calc_ψ_h)
        @test abs(ψ(model, -1e-9)) < 1e-7 && abs(ψ(model, 1e-9)) < 1e-7
    end
    @test calc_Φ_h(model, -1e-9) ≈ calc_Φ_h(model, 1e-9) atol = 1e-7
end

@testset "$(pair_name(model)): Φ_h stays bounded in strongly stable air, keeping RaupachLTheoryAirProfile's T_L bounded" for model in
        (MoninObukhov(), MoninObukhov(; unstable_stability=NicheMapRUnstable(), stable_stability=NicheMapRStable()))
    for L in (10.0u"m", 1.0u"m", 0.1u"m", 0.01u"m", 0.0001u"m")
        @test 1.0 <= calc_Φ_h(model, ζ(L)) <= 10.0
    end
end

@testset "$(nameof(typeof(f))): Φ_h is the heat gradient of ψ_h, 1 - ζ dψ_h/dζ" for f in
        (ChengBrutsaert(), BeljaarsHoltslag(), BusingerLinear(), PaulsonDyer())
    for ζ₀ in (0.01, 0.5, 2.0, 20.0)
        stability = f isa PaulsonDyer ? -ζ₀ : ζ₀
        h = 1e-6
        slope = (calc_ψ_h(f, stability + h) - calc_ψ_h(f, stability - h)) / 2h
        @test calc_Φ_h(f, stability) ≈ 1 - stability * slope rtol = 1e-6
    end
end

@testset "$(nameof(typeof(f))): smooth and strictly decreasing through the stable range" for f in (ChengBrutsaert(), BeljaarsHoltslag())
    ζs = range(0.1, 50.0; length=2001)
    for ψ in (calc_ψ_m, calc_ψ_h)
        values = ψ.(Ref(f), ζs)
        @test all(<(0), diff(values))
        # No kinks: the slope changes gradually everywhere (NicheMapRStable's clamp jumps it by ~6)
        slopes = diff(values) ./ step(ζs)
        @test maximum(abs, diff(slopes)) < 0.5
    end
end

@testset "$(pair_name(bl)) atmospheric_surface_profile!: no crash, sensible output across regimes" for bl in stability_pairs
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

@testset "MomentumCorrectionLimit bounds the stability where the momentum correction reaches its fraction of the log term" begin
    limit = MomentumCorrectionLimit()
    limits = Microclimate.inverse_obukhov_length_limits(limit, bl, z, z0)
    correction(inverse_length) = calc_ψ_m(bl, z0 * inverse_length) - calc_ψ_m(bl, z * inverse_length)
    @test limits.unstable < 0.0u"m^-1" < limits.stable
    @test isfinite(limits.unstable) && isfinite(limits.stable)
    @test correction(limits.stable) ≈ limit.fraction * log(z / z0) rtol = 1e-6
    @test -correction(limits.unstable) ≈ limit.fraction * log(z / z0) rtol = 1e-6
    for inverse_length in (-10.0, -0.01, 0.01, 10.0) .* u"m^-1"
        limited = Microclimate.limit_inverse_obukhov_length(inverse_length, limits)
        @test limits.unstable <= limited <= limits.stable
        @test sign(limited) == sign(inverse_length)
    end
    unlimited = Microclimate.inverse_obukhov_length_limits(NoStabilityLimit(), bl, z, z0)
    @test Microclimate.limit_inverse_obukhov_length(10.0u"m^-1", unlimited) == 10.0u"m^-1"
end

@testset "atmospheric_surface_profile!: calm air keeps the log-law wind profile" begin
    site = Site(; latitude=0.0u"°", longitude=0.0u"°", elevation=0.0u"m", slope=0.0u"°", aspect=0.0u"°",
        horizon_angles=fill(0.0u"°", 24), sky_view_fraction=1.0, albedo=0.15, roughness_height=0.01u"m",
        atmospheric_pressure=101325.0u"Pa")
    heights = [0.5, 10.0]u"m"
    buffers = allocate_profile(bl, heights)
    environment_instant = (; atmospheric_pressure=101325.0u"Pa", reference_temperature=293.0u"K",
        reference_wind_speed=0.2u"m/s", reference_humidity=0.6, zenith_angle=30.0u"°")
    result = atmospheric_surface_profile!(bl, buffers; site, environment_instant, surface_temperature=293.0u"K")
    @test result.friction_velocity == bl.min_friction_velocity
    # The iteration stops within its tolerance of neutral, 0.01 in z/L
    @test result.wind_speed[1] ≈ 0.2u"m/s" * log(0.5 / 0.01) / log(10.0 / 0.01) rtol = 1e-2
end
