using Microclimate
using Unitful
using Test
using Microclimate: allocate_soil_water_balance, infiltration_step!,
    RainfallEntryStep, rainfall_entry!, soil_water_storage,
    allocate_phase_transition, frozen_water_content!,
    conductivity_factor, conductivity_blocked, LATENT_HEAT_FUSION

function _test_profile()
    profile = example_soil_profile(; bulk_density=1.3u"Mg/m^3", mineral_density=2.65u"Mg/m^3")
    num_layers = length(profile.hydraulics.campbell_b_parameter)
    depths = fill(0.15u"m", num_layers) .* (1:num_layers)
    saturation = 1 - profile.bulk_density[1] / profile.mineral_density[1]
    (; profile, depths, num_layers, saturation)
end

_entry_step(mode, profile, depths, pool, dt, frozen_water_content; frozen_conductivity=BloomsburgWangFrozenConductivity()) =
    RainfallEntryStep(example_soil_hydraulic_model(; rainfall_entry_mode=mode, frozen_conductivity); soil_profile=profile,
        depths, pool, moisture_timestep=dt, niter_moist=10, max_surface_pool=1.0e4u"kg/m^2", frozen_water_content)

@testset "BloomsburgWangFrozenConductivity" begin
    model = BloomsburgWangFrozenConductivity()
    (; min_ice_free_porosity, conductivity_floor) = model
    porosity = 0.5

    @test conductivity_factor(model, 0.0, porosity) === 1.0

    # Monotonic non-increasing as ice content rises.
    ice_contents = 0.0:0.02:porosity
    @test issorted([conductivity_factor(model, ice, porosity) for ice in ice_contents], rev=true)

    # Below/at/above the Bloomsburg & Wang threshold.
    just_above = porosity - min_ice_free_porosity - 1e-6
    at = porosity - min_ice_free_porosity
    just_below = porosity - min_ice_free_porosity + 1e-6
    @test conductivity_factor(model, just_above, porosity) > conductivity_floor
    @test conductivity_factor(model, at, porosity) == conductivity_floor
    @test conductivity_factor(model, just_below, porosity) == conductivity_floor
    @test !conductivity_blocked(model, just_above, porosity)
    @test conductivity_blocked(model, at, porosity)

    # Ice content exceeding porosity (shouldn't happen, but must stay finite).
    @test conductivity_factor(model, 2porosity, porosity) == conductivity_floor
    @test conductivity_blocked(model, 2porosity, porosity)
end

@testset "NoFrozenConductivity" begin
    @test conductivity_factor(NoFrozenConductivity(), 0.4, 0.5) == 1.0
    @test !conductivity_blocked(NoFrozenConductivity(), 0.5, 0.5)
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
        ice = frozen_water_content!(model, buffers, accumulated, soil_moisture)
        @test ice ≈ [0.10, 0.10, 0.0, 0.40]
    end

    @testset "clamped to [0,1] even if accumulated exceeds the budget" begin
        buffers.layer_mass .= [1.0, 2.0, 3.0, 4.0]u"kg"
        accumulated = 1.5 .* (LATENT_HEAT_FUSION .* buffers.layer_mass)
        ice = frozen_water_content!(model, buffers, accumulated, soil_moisture)
        @test ice ≈ soil_moisture
    end

    @testset "freeze then thaw updates content back down" begin
        accumulated = copy(LATENT_HEAT_FUSION .* buffers.layer_mass)  # fully frozen
        frozen = frozen_water_content!(model, buffers, accumulated, soil_moisture)
        @test frozen ≈ soil_moisture
        accumulated .*= 0.0  # thawed
        thawed = frozen_water_content!(model, buffers, accumulated, soil_moisture)
        @test all(iszero, thawed)
    end
end

@testset "frozen soil suppresses the Darcy solve" begin
    (; profile, depths, num_layers, saturation) = _test_profile()
    input_soil_temperature = fill(293.0u"K", num_layers)
    function run(frozen_water_content; frozen_conductivity=BloomsburgWangFrozenConductivity(), soil_moisture=fill(0.30, num_layers))
        model = example_soil_hydraulic_model(; frozen_conductivity)
        buffers = allocate_soil_water_balance(model, num_layers)
        infiltration_step!(buffers, model;
            soil_profile=profile, depths, atmospheric_pressure=101325.0u"Pa",
            local_relative_humidity=0.5, leaf_area_index=1.0u"Mg/m^3",
            soil_moisture, evapotranspiration=0.0u"kg/m^2/s", input_soil_temperature,
            moisture_timestep=360.0u"s", moisture_tolerance=1e-6u"kg/m^2/s", moisture_max_iterations=200,
            frozen_water_content,
        )
    end

    @testset "all layers frozen: finite, no NaN" begin
        out = run(fill(saturation - 0.02, num_layers))
        @test !any(isnan, ustrip.(out.soil_water_potential))
        @test !any(isnan, out.soil_moisture)
        @test isfinite(ustrip(out.drainage))
    end

    @testset "deepest real layer frozen suppresses drainage; ghost boundary untouched by design" begin
        # Physical index num_layers-1 is the deepest *real* solved layer
        # (buffer index num_layers); physical index num_layers is the fixed
        # saturated free-drainage ghost boundary and is never ice-impeded.
        frozen_water_content = zeros(num_layers)
        frozen_water_content[num_layers-1] = saturation - 0.05
        unfrozen = run(frozen_water_content; frozen_conductivity=NoFrozenConductivity())
        frozen = run(frozen_water_content)
        @test !any(isnan, frozen.soil_moisture)
        @test frozen.drainage < unfrozen.drainage / 1000
    end
end

@testset "frozen soil blocks rainfall entry" begin
    (; profile, depths, num_layers, saturation) = _test_profile()
    rain = 8.0u"kg/m^2"
    dt = 360.0u"s"
    wet = fill(saturation - 0.10, num_layers)  # frozen_water_content=wet leaves ~0.10 ice-free porosity, below threshold

    @testset "$mode: near-saturated + fully frozen leaves rain in the pool" for mode in (PoolCapacityRainfall(), RateLimitedFrontRainfall())
        soil_moisture_unfrozen = copy(wet)
        remaining_unfrozen = rainfall_entry!(soil_moisture_unfrozen, mode, rain, _entry_step(mode, profile, depths, rain, dt, wet; frozen_conductivity=NoFrozenConductivity()))

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
        unfrozen = _entry_step(mode, profile, depths, rain, dt, wet; frozen_conductivity=NoFrozenConductivity())
        frozen = _entry_step(mode, profile, depths, rain, dt, wet)
        @test frozen.rainfall_flux < unfrozen.rainfall_flux * 1e-3
    end
end
