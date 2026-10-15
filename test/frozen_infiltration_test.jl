using Microclimate
using Unitful
using Test
using Microclimate: allocate_soil_water_balance, infiltration_step!,
    RainfallEntryStep, rainfall_entry!, soil_water_storage,
    allocate_phase_transition, frozen_water_content!,
    ice_impeded_conductivity, ice_free_capacity, NoIce,
    ICE_IMPEDANCE_MIN_POROSITY, ICE_CONDUCTIVITY_FLOOR_FACTOR, LATENT_HEAT_FUSION

function _test_profile()
    profile = example_soil_profile(; bulk_density=1.3u"Mg/m^3", mineral_density=2.65u"Mg/m^3")
    num_layers = length(profile.hydraulics.campbell_b_parameter)
    depths = fill(0.15u"m", num_layers) .* (1:num_layers)
    sat = 1 - profile.bulk_density[1] / profile.mineral_density[1]
    (; profile, depths, num_layers, sat)
end

_entry_step(mode, profile, depths, pool, dt, frozen_water_content) =
    RainfallEntryStep(mode; soil_profile=profile, depths, pool, moisture_timestep=dt, niter_moist=10,
        max_surface_pool=1.0e4u"kg/m^2", frozen_water_content)

@testset "ice_impeded_conductivity" begin
    K = 5.0u"kg*s/m^3"
    sat = 0.5

    @test ice_impeded_conductivity(K, 0.0, sat) === K  # NoIce regression: exact identity

    # Monotonic non-increasing as ice content rises.
    ice_contents = 0.0:0.02:sat
    factors = [ustrip(ice_impeded_conductivity(K, ic, sat) / K) for ic in ice_contents]
    @test issorted(factors, rev=true)

    # Below/at/above the Bloomsburg & Wang threshold.
    just_above = sat - ICE_IMPEDANCE_MIN_POROSITY - 1e-6
    at = sat - ICE_IMPEDANCE_MIN_POROSITY
    just_below = sat - ICE_IMPEDANCE_MIN_POROSITY + 1e-6
    @test ustrip(ice_impeded_conductivity(K, just_above, sat) / K) > ICE_CONDUCTIVITY_FLOOR_FACTOR
    @test ice_impeded_conductivity(K, at, sat) ≈ K * ICE_CONDUCTIVITY_FLOOR_FACTOR
    @test ice_impeded_conductivity(K, just_below, sat) ≈ K * ICE_CONDUCTIVITY_FLOOR_FACTOR

    # Ice content exceeding porosity (shouldn't happen, but must stay finite).
    @test isfinite(ustrip(ice_impeded_conductivity(K, 2*sat, sat)))
    @test ice_impeded_conductivity(K, 2*sat, sat) ≈ K * ICE_CONDUCTIVITY_FLOOR_FACTOR

    @test ice_free_capacity(sat, 0.0) == sat
    @test ice_free_capacity(sat, 2*sat) == 0.0  # clamped, not negative
end

@testset "frozen_water_content!" begin
    model = PhaseTransitionLatentHeat()
    n = 4
    buffers = allocate_phase_transition(model, n)
    soil_moisture = [0.10, 0.20, 0.30, 0.40]  # distinct per layer -- verifies no index shift
    buffers.layer_mass .= [1.0, 2.0, 3.0, 4.0]u"kg"

    @testset "fraction scales soil_moisture per layer, no NaN when dry" begin
        accumulated = [1.0, 0.5, 0.0, 1.0] .* (LATENT_HEAT_FUSION .* buffers.layer_mass)
        buffers.layer_mass[3] = 0.0u"kg"  # dry layer: max_latent_heat == 0
        fwc = frozen_water_content!(model, buffers, accumulated, soil_moisture)
        @test fwc ≈ [0.10, 0.10, 0.0, 0.40]
    end

    @testset "clamped to [0,1] even if accumulated exceeds the budget" begin
        buffers.layer_mass .= [1.0, 2.0, 3.0, 4.0]u"kg"
        accumulated = 1.5 .* (LATENT_HEAT_FUSION .* buffers.layer_mass)
        fwc = frozen_water_content!(model, buffers, accumulated, soil_moisture)
        @test fwc ≈ soil_moisture
    end

    @testset "freeze then thaw updates content back down" begin
        accumulated = copy(LATENT_HEAT_FUSION .* buffers.layer_mass)  # fully frozen
        fwc_frozen = frozen_water_content!(model, buffers, accumulated, soil_moisture)
        @test fwc_frozen ≈ soil_moisture
        accumulated .*= 0.0  # thawed
        fwc_thawed = frozen_water_content!(model, buffers, accumulated, soil_moisture)
        @test all(iszero, fwc_thawed)
    end
end

@testset "frozen soil suppresses the Darcy solve" begin
    (; profile, depths, num_layers, sat) = _test_profile()
    input_soil_temperature = fill(293.0u"K", num_layers)
    model = example_soil_hydraulic_model()
    buffers = allocate_soil_water_balance(model, num_layers)

    function run(frozen_water_content; soil_moisture=fill(0.30, num_layers))
        infiltration_step!(buffers, model;
            soil_profile=profile, depths, atmospheric_pressure=101325.0u"Pa",
            local_relative_humidity=0.5, leaf_area_index=1.0u"Mg/m^3",
            soil_moisture, evapotranspiration=0.0u"kg/m^2/s", input_soil_temperature,
            moisture_timestep=360.0u"s", moisture_tolerance=1e-6u"kg/m^2/s", moisture_max_iterations=200,
            frozen_water_content,
        )
    end

    @testset "all layers frozen: finite, no NaN" begin
        out = run(fill(sat - 0.02, num_layers))
        @test !any(isnan, ustrip.(out.soil_water_potential))
        @test !any(isnan, out.soil_moisture)
        @test isfinite(ustrip(out.drainage))
    end

    @testset "deepest real layer frozen suppresses drainage; ghost boundary untouched by design" begin
        # Physical index num_layers-1 is the deepest *real* solved layer
        # (buffer index num_layers); physical index num_layers is the fixed
        # saturated free-drainage ghost boundary and is never ice-impeded.
        fwc = zeros(num_layers)
        fwc[num_layers-1] = sat - 0.05
        unfrozen = run(NoIce())
        frozen = run(fwc)
        @test !any(isnan, frozen.soil_moisture)
        @test frozen.drainage < unfrozen.drainage / 1000
    end
end

@testset "frozen soil blocks rainfall entry" begin
    (; profile, depths, num_layers, sat) = _test_profile()
    rain = 8.0u"kg/m^2"
    dt = 360.0u"s"
    wet = fill(sat - 0.10, num_layers)  # frozen_water_content=wet leaves ~0.10 ice-free porosity, below threshold

    @testset "$mode: near-saturated + fully frozen leaves rain in the pool" for mode in (PoolCapacityRainfall(), RateLimitedFrontRainfall())
        soil_moisture_unfrozen = copy(wet)
        remaining_unfrozen = rainfall_entry!(soil_moisture_unfrozen, mode, rain, _entry_step(mode, profile, depths, rain, dt, NoIce()))

        soil_moisture_frozen = copy(wet)
        remaining_frozen = rainfall_entry!(soil_moisture_frozen, mode, rain, _entry_step(mode, profile, depths, rain, dt, wet))

        @test remaining_frozen > remaining_unfrozen
        @test ustrip(u"kg/m^2", rain - remaining_frozen) < ustrip(u"kg/m^2", rain) * 1e-3  # essentially none infiltrated
    end

    @testset "$mode: mass conservation holds while frozen" for mode in (PoolCapacityRainfall(), RateLimitedFrontRainfall())
        step = _entry_step(mode, profile, depths, rain, dt, wet)
        soil_moisture = copy(wet)
        storage_before = soil_water_storage(soil_moisture, step)
        remaining_pool = rainfall_entry!(soil_moisture, mode, rain, step)
        storage_after = soil_water_storage(soil_moisture, step)
        @test isapprox(ustrip(u"kg/m^2", storage_after - storage_before + remaining_pool - rain), 0.0; atol=1e-6)
    end

    @testset "ImplicitFluxRainfall: ice impedes the boundary flux, rainfall_entry! is a no-op" begin
        mode = ImplicitFluxRainfall()
        soil_moisture = copy(wet)
        pool = rainfall_entry!(soil_moisture, mode, rain, _entry_step(mode, profile, depths, rain, dt, wet))
        @test pool == rain
        @test soil_moisture == wet
        unfrozen = _entry_step(mode, profile, depths, rain, dt, NoIce())
        frozen = _entry_step(mode, profile, depths, rain, dt, wet)
        @test frozen.rainfall_flux < unfrozen.rainfall_flux * 1e-3
    end
end
