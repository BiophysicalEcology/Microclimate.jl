using Microclimate
using Unitful
using Test
using Microclimate: allocate_soil_water_balance, infiltration_substeps!,
    rainfall_entry!, rainfall_flux_for_step, soil_water_storage

# rainfall_entry_mode's interface functions are internal, same rationale as
# infiltration_algorithm_comparison.jl's own white-box testing.

function _test_profile()
    profile = example_soil_profile(; bulk_density=1.3u"Mg/m^3", mineral_density=2.65u"Mg/m^3")
    num_layers = length(profile.hydraulics.campbell_b_parameter)
    depths = fill(0.05u"m", num_layers) .* (1:num_layers)  # 0.05, 0.10, ..., dense near-surface spacing
    (; profile, depths, num_layers)
end

# One hour of rainfall entry and infiltration with evaporation held at zero,
# so the only mass terms are rain, pool and storage.
_moisture_mode() = DynamicSoilMoisture(; moisture_timestep=360.0u"s", moisture_tolerance=1e-6u"kg/m^2/s", moisture_max_iterations=200)

function run_rainfall_hour(mode; soil_moisture, pool)
    (; profile, depths, num_layers) = _test_profile()
    model = example_soil_hydraulic_model(; rainfall_entry_mode=mode)
    moisture_mode = _moisture_mode()
    buffers = allocate_soil_water_balance(model, num_layers)
    rainfall_flux = rainfall_flux_for_step(mode, moisture_mode, pool)
    (; soil_moisture, pool) = infiltration_substeps!(buffers, model, moisture_mode;
        soil_profile=profile, depths, soil_moisture, pool, rainfall_flux, frozen_water_content=zeros(num_layers),
        atmospheric_pressure=101325.0u"Pa", local_relative_humidity=0.5, leaf_area_index=1.0u"Mg/m^3",
        evapotranspiration=0.0u"kg/m^2/s", input_soil_temperature=fill(293.0u"K", num_layers),
    )
    return (; soil_moisture, pool, storage=soil_water_storage(soil_moisture, depths))
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

    # rainfall_entry! alone, no infiltration_step! -- isolates the
    # rainfall entry from the fixed-water-table boundary's own storage change.
    @testset "rainfall_entry! conserves mass" begin
        (; profile, depths, num_layers) = _test_profile()
        rain = 8.0u"kg/m^2"
        moisture_mode = _moisture_mode()
        for mode in (PoolCapacityRainfall(), RateLimitedFrontRainfall())
            soil_moisture = fill(0.10, num_layers)
            storage_before = soil_water_storage(soil_moisture, depths)
            remaining_pool = rainfall_entry!(soil_moisture, mode, example_soil_hydraulic_model(; rainfall_entry_mode=mode), moisture_mode, rain;
                soil_profile=profile, depths, frozen_water_content=zeros(num_layers))
            storage_after = soil_water_storage(soil_moisture, depths)
            @test isapprox(ustrip(u"kg/m^2", storage_after - storage_before + remaining_pool - rain), 0.0; atol=1e-6)
        end
        # ImplicitFluxRainfall never mutates soil_moisture here -- its flux
        # rate integrates back to the original pool over the hour instead.
        rate = rainfall_flux_for_step(ImplicitFluxRainfall(), moisture_mode, rain)
        @test isapprox(ustrip(u"kg/m^2", rate * 1u"hr"), ustrip(u"kg/m^2", rain); atol=1e-9)
    end

    @testset "end-to-end smoke test per mode" begin
        for mode in modes
            soil_hydraulic_model = example_soil_hydraulic_model(; rainfall_entry_mode=mode)
            soil_moisture_strategy = DynamicSoilMoisture()
            soil_profile = example_soil_profile(; bulk_density=1.3u"Mg/m^3", mineral_density=2.65u"Mg/m^3")
            out = solve(example_microclimate_problem(; soil_hydraulic_model, soil_moisture_strategy, soil_profile))
            @test !any(isnan, ustrip.(out.soil_temperature))
            @test !any(isnan, out.soil_moisture)
            @test all(0 .<= out.soil_moisture .<= 1)
            @test !any(isnan, ustrip.(out.runoff))
            @test all(out.runoff .>= 0u"kg/m^2")
        end
    end
end

@testset "runoff tracking" begin
    @testset "default max_surface_pool: runoff stays exactly zero" begin
        out = solve(example_microclimate_problem(; soil_moisture_strategy=DynamicSoilMoisture()))
        @test all(iszero, out.runoff)
    end

    @testset "small max_surface_pool generates nonzero, non-decreasing runoff" begin
        soil_moisture_strategy = DynamicSoilMoisture()
        max_surface_pool = 0.01u"kg/m^2"
        out = solve(example_microclimate_problem(; soil_moisture_strategy, max_surface_pool))
        @test out.runoff[end] > 0u"kg/m^2"
        @test all(diff(ustrip.(u"kg/m^2", out.runoff)) .>= -1e-9)  # cumulative, monotonically non-decreasing
    end

    @testset "runoff still accumulates under PrescribedSoilMoisture (no infiltration solve)" begin
        max_surface_pool = 0.01u"kg/m^2"  # soil_moisture_strategy defaults to PrescribedSoilMoisture
        out = solve(example_microclimate_problem(; max_surface_pool))
        @test out.runoff[end] > 0u"kg/m^2"
    end

    @testset "PoolCapacityRainfall() explicit == default (regression guard)" begin
        soil_moisture_strategy = DynamicSoilMoisture()
        soil_profile = example_soil_profile(; bulk_density=1.3u"Mg/m^3", mineral_density=2.65u"Mg/m^3")
        out_default = solve(example_microclimate_problem(; soil_hydraulic_model=example_soil_hydraulic_model(), soil_moisture_strategy, soil_profile))
        out_explicit = solve(example_microclimate_problem(; soil_hydraulic_model=example_soil_hydraulic_model(; rainfall_entry_mode=PoolCapacityRainfall()), soil_moisture_strategy, soil_profile))
        @test out_default.soil_moisture == out_explicit.soil_moisture
        @test out_default.soil_temperature == out_explicit.soil_temperature
    end
end
