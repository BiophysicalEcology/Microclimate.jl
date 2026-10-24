using Microclimate
using Microclimate: allocate_canopy, canopy_energy_balance!
using Unitful
using Test

heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"

@testset "MultilayerCanopy wired into the main solve loop: finite, physically sane" begin
    canopy_model = example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 2.0)
    problem = example_microclimate_problem(; heights, canopy_model)
    output = solve(problem)

    @test all(isfinite, ustrip.(u"K", output.soil_temperature))
    @test all(t -> 200.0u"K" < t < 340.0u"K", output.soil_temperature)
    @test all(isfinite, output.soil_moisture)

    @test size(output.canopy.leaf_temperature, 2) == 10
    @test all(isfinite, ustrip.(u"K", output.canopy.leaf_temperature))
    @test all(t -> 200.0u"K" < t < 340.0u"K", output.canopy.leaf_temperature)
    @test all(isfinite, ustrip.(u"K", output.canopy.air_temperature))
    @test all(isfinite, ustrip.(u"m/s", output.canopy.wind_speed))
    @test all(>=(0.0u"W/m^2"), output.canopy.ground_absorbed_shortwave)
    @test all(>(0), output.canopy.iterations)

    @test size(output.canopy.boundary_downward_shortwave, 2) == 11
    @test size(output.canopy.boundary_downward_longwave, 2) == 11
    @test all(>=(0.0u"W/m^2"), output.canopy.boundary_downward_shortwave)
    @test all(>=(0.0u"W/m^2"), output.canopy.boundary_upward_shortwave)
    @test all(>=(0.0u"W/m^2"), output.canopy.boundary_downward_longwave)
    @test all(>=(0.0u"W/m^2"), output.canopy.boundary_upward_longwave)
    # canopy top (boundary 1) receives at least as much downward shortwave as the ground (boundary end)
    @test all(output.canopy.boundary_downward_shortwave[:, 1] .>= output.canopy.boundary_downward_shortwave[:, end])
end

@testset "NoCanopy leaves canopy output zero-columned" begin
    output = solve(example_microclimate_problem(; heights))
    @test size(output.canopy.leaf_temperature) == (size(output.soil_temperature, 1), 0)
    @test all(iszero, output.canopy.ground_absorbed_shortwave)
    @test size(output.canopy.boundary_downward_shortwave) == (size(output.soil_temperature, 1), 0)
    @test size(output.canopy.boundary_upward_longwave) == (size(output.soil_temperature, 1), 0)
end

@testset "MultilayerCanopy changes soil temperature relative to NoCanopy" begin
    canopy_model = example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 3.0)
    canopy_output = solve(example_microclimate_problem(; heights, canopy_model))
    bare_output = solve(example_microclimate_problem(; heights))

    @test canopy_output.soil_temperature != bare_output.soil_temperature
end

@testset "MultilayerCanopy leaf temperature responds to soil moisture via Campbell coupling" begin
    canopy_model = example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 3.0)
    soil_moisture_strategy = DynamicSoilMoisture()
    n = length(Microclimate.DEFAULT_DEPTHS)

    wet_output = solve(example_microclimate_problem(; heights, canopy_model, soil_moisture_strategy,
        initial_soil_moisture = fill(0.35, n)))
    dry_output = solve(example_microclimate_problem(; heights, canopy_model, soil_moisture_strategy,
        initial_soil_moisture = fill(0.03, n)))

    @test all(isfinite, ustrip.(u"K", wet_output.canopy.leaf_temperature))
    @test all(isfinite, ustrip.(u"K", dry_output.canopy.leaf_temperature))
    # Different initial soil moisture propagates through Campbell's demand/supply
    # solve into a different leaf_water_potential and hence leaf temperature. The
    # *direction* of that difference over a multi-day run also depends on ground-
    # temperature feedback and how quickly soil moisture re-equilibrates (verified
    # separately, directly on canopy_energy_balance!, that more negative
    # leaf_water_potential warms the leaf — see canopy_energy_balance_picard_test.jl).
    @test wet_output.canopy.leaf_temperature != dry_output.canopy.leaf_temperature
end

@testset "SurfaceNodePrediction tracks the iterated canopy-soil coupling" begin
    canopy_model = example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 2.0)
    run(coupling) = solve(example_microclimate_problem(; heights, canopy_model,
        config = MicroConfig(; canopy_soil_convergence = coupling)))
    predicted = run(SurfaceNodePrediction())
    iterated = run(IterationToleranceConvergence(; tolerance = 0.1u"K", max_iterations_per_day = 80))
    difference = ustrip.(u"K", predicted.soil_temperature[:, 1] .- iterated.soil_temperature[:, 1])
    @test sqrt(sum(abs2, difference) / length(difference)) < 0.3
end
