using Microclimate
using Unitful
using Test
using Microclimate: allocate_soil_water_balance, infiltration_step!,
    rainfall_flux_for_step, apply_rainfall_entry!, post_infiltration_pool_update

# rainfall_entry_mode's interface functions are internal, same rationale as
# infiltration_algorithm_comparison.jl's own white-box testing.

function _test_profile()
    profile = example_soil_profile(; bulk_density=1.3u"Mg/m^3", mineral_density=2.65u"Mg/m^3")
    num_layers = length(profile.hydraulics.campbell_b_parameter)
    depths = fill(0.05u"m", num_layers) .* (1:num_layers)  # 0.05, 0.10, ..., dense near-surface spacing
    (; profile, depths, num_layers)
end

_layer_thickness(depths, half_thickness, i, n) =
    i == 1 ? half_thickness : i == n ? (depths[n] - depths[n-1]) / 2 : (depths[i+1] - depths[i-1]) / 2

# One hour's worth of rainfall_entry_mode + infiltration_step! sub-stepping,
# mirroring soil_water_balance!'s own loop but with evaporation held at zero
# so the mass ledger only has rain/pool/storage/runoff terms to track.
function run_rainfall_hour(mode; soil_moisture, pool, dt=360.0u"s", niter_moist=10)
    (; profile, depths, num_layers) = _test_profile()
    input_soil_temperature = fill(293.0u"K", num_layers)
    model = example_soil_hydraulic_model(; rainfall_entry_mode=mode)
    buffers = allocate_soil_water_balance(model, num_layers)
    sat = 1 - profile.bulk_density[1] / profile.mineral_density[1]
    half_thickness = depths[1] / 2

    rainfall_flux_rate = rainfall_flux_for_step(mode, pool, dt, niter_moist)
    pool = apply_rainfall_entry!(mode, soil_moisture, pool, sat, half_thickness, rainfall_flux_rate, dt; depths, soil_profile=profile)
    local out
    for _ in 1:niter_moist
        out = infiltration_step!(buffers, model;
            soil_profile=profile, depths, atmospheric_pressure=101325.0u"Pa",
            local_relative_humidity=0.5, leaf_area_index=1.0u"Mg/m^3",
            soil_moisture, evapotranspiration=0.0u"kg/m^2/s", input_soil_temperature,
            moisture_timestep=dt, moisture_tolerance=1e-6u"kg/m^2/s", moisture_max_iterations=200,
            rainfall_flux_rate,
        )
        soil_moisture = out.soil_moisture
        water_flux = max(0.0u"kg/m^2", out.surface_water_flux)
        pool = post_infiltration_pool_update(mode, pool, rainfall_flux_rate, dt, water_flux, 0.0u"kg/m^2", 1.0e4u"kg/m^2")
        pool = apply_rainfall_entry!(mode, soil_moisture, pool, sat, half_thickness, rainfall_flux_rate, dt; depths, soil_profile=profile)
    end
    storage = sum(i -> soil_moisture[i] * _layer_thickness(depths, half_thickness, i, num_layers), 1:num_layers) * 1000.0u"kg/m^3"
    return (; soil_moisture, pool, storage)
end

@testset "rainfall entry modes" begin
    modes = (PoolCapacityRainfall(), ImplicitFluxRainfall(), RateLimitedFrontRainfall())
    n = _test_profile().num_layers

    @testset "sign convention: rain increases moisture, not decreases it" begin
        for mode in modes
            dry = fill(0.10, n)
            wet_result = run_rainfall_hour(mode; soil_moisture=copy(dry), pool=5.0u"kg/m^2")
            dry_result = run_rainfall_hour(mode; soil_moisture=copy(dry), pool=0.0u"kg/m^2")
            @test !any(isnan, wet_result.soil_moisture)
            @test sum(wet_result.soil_moisture) > sum(dry_result.soil_moisture)
        end
    end

    # apply_rainfall_entry! alone, no infiltration_step! -- isolates the new
    # rainfall-entry mechanism from the fixed-water-table boundary's own
    # (intentional, pre-existing) storage contribution.
    @testset "apply_rainfall_entry! conserves mass" begin
        (; profile, depths, num_layers) = _test_profile()
        sat = 1 - profile.bulk_density[1] / profile.mineral_density[1]
        half_thickness = depths[1] / 2
        rain = 8.0u"kg/m^2"
        dt = 360.0u"s"
        for mode in (PoolCapacityRainfall(), RateLimitedFrontRainfall())
            soil_moisture = fill(0.10, num_layers)
            storage_before = sum(i -> soil_moisture[i] * _layer_thickness(depths, half_thickness, i, num_layers), 1:num_layers) * 1000.0u"kg/m^3"
            remaining_pool = apply_rainfall_entry!(mode, soil_moisture, rain, sat, half_thickness, 0.0u"kg/m^2/s", dt; depths, soil_profile=profile)
            storage_after = sum(i -> soil_moisture[i] * _layer_thickness(depths, half_thickness, i, num_layers), 1:num_layers) * 1000.0u"kg/m^3"
            @test isapprox(ustrip(u"kg/m^2", storage_after - storage_before + remaining_pool - rain), 0.0; atol=1e-6)
        end
        # ImplicitFluxRainfall never mutates soil_moisture here -- its flux
        # rate integrates back to the original pool over the hour instead.
        rate = rainfall_flux_for_step(ImplicitFluxRainfall(), rain, dt, 10)
        @test isapprox(ustrip(u"kg/m^2", rate * 10 * dt), ustrip(u"kg/m^2", rain); atol=1e-9)
    end

    @testset "end-to-end smoke test per mode" begin
        for mode in modes
            soil_hydraulic_model = example_soil_hydraulic_model(; rainfall_entry_mode=mode)
            config = MicroConfig(; soil_moisture_strategy=DynamicSoilMoisture())
            soil_profile = example_soil_profile(; bulk_density=1.3u"Mg/m^3", mineral_density=2.65u"Mg/m^3")
            out = solve(example_microclimate_problem(; soil_hydraulic_model, config, soil_profile))
            @test !any(isnan, ustrip.(out.soil_temperature))
            @test !any(isnan, out.soil_moisture)
            @test all(0 .<= out.soil_moisture .<= 1)
            @test !any(isnan, ustrip.(out.runoff))
            @test all(out.runoff .>= 0u"kg/m^2")
        end
    end
end

@testset "runoff tracking" begin
    @testset "default max_pond_depth: runoff stays exactly zero (backward compat)" begin
        out = solve(example_microclimate_problem(; config=MicroConfig(; soil_moisture_strategy=DynamicSoilMoisture())))
        @test all(iszero, out.runoff)
    end

    @testset "small max_pond_depth generates nonzero, non-decreasing runoff" begin
        config = MicroConfig(; soil_moisture_strategy=DynamicSoilMoisture(), max_pond_depth=0.01u"kg/m^2")
        out = solve(example_microclimate_problem(; config))
        @test out.runoff[end] > 0u"kg/m^2"
        @test all(diff(ustrip.(u"kg/m^2", out.runoff)) .>= -1e-9)  # cumulative, monotonically non-decreasing
    end

    @testset "runoff still accumulates under PrescribedSoilMoisture (no infiltration solve)" begin
        config = MicroConfig(; max_pond_depth=0.01u"kg/m^2")  # soil_moisture_strategy defaults to PrescribedSoilMoisture
        out = solve(example_microclimate_problem(; config))
        @test out.runoff[end] > 0u"kg/m^2"
    end

    @testset "PoolCapacityRainfall() explicit == default (regression guard)" begin
        config = MicroConfig(; soil_moisture_strategy=DynamicSoilMoisture())
        soil_profile = example_soil_profile(; bulk_density=1.3u"Mg/m^3", mineral_density=2.65u"Mg/m^3")
        out_default = solve(example_microclimate_problem(; soil_hydraulic_model=example_soil_hydraulic_model(), config, soil_profile))
        out_explicit = solve(example_microclimate_problem(; soil_hydraulic_model=example_soil_hydraulic_model(; rainfall_entry_mode=PoolCapacityRainfall()), config, soil_profile))
        @test out_default.soil_moisture == out_explicit.soil_moisture
        @test out_default.soil_temperature == out_explicit.soil_temperature
    end
end
