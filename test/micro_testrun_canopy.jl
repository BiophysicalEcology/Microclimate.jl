using Microclimate
using Microclimate: allocate_canopy, canopy_energy_balance!
using Unitful
using Test
using ADTypes: AutoFiniteDiff
using OrdinaryDiffEqSDIRK: ImplicitEuler

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

@testset "An implicit soil solver follows the default, with or without a canopy" begin
    implicit_soil = SoilHeatTransport1D(; ode_solver = ImplicitEuler(; autodiff = AutoFiniteDiff()),
        ode_kwargs = (; adaptive = false, dt = 60.0u"minute"))
    for canopy_model in (NoCanopy(), example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 2.0))
        surface_temperature(soil_energy_model) = solve(example_microclimate_problem(; heights, canopy_model,
            soil_energy_model)).soil_temperature[:, 1]
        @test isapprox(surface_temperature(implicit_soil), surface_temperature(SoilHeatTransport1D()); rtol = 2e-3)
    end
end

@testset "solve! allocates nothing, whatever the canopy and soil solver" begin
    function solve_allocations(cache)
        solve!(cache)
        solve!(cache)
        return @allocated solve!(cache)
    end
    implicit_soil = SoilHeatTransport1D(; ode_solver = ImplicitEuler(; autodiff = AutoFiniteDiff()),
        ode_kwargs = (; adaptive = false, dt = 60.0u"minute"))
    for canopy_model in (NoCanopy(), example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 2.0)),
        soil_energy_model in (SoilHeatTransport1D(), implicit_soil)
        cache = init(example_microclimate_problem(; heights, canopy_model, soil_energy_model))
        @test solve_allocations(cache) == 0
    end
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

@testset "reinit! with a canopy: same result as a fresh solve" begin
    problem = example_microclimate_problem(; heights,
        canopy_model = example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 3.0))
    fresh = solve(problem)
    cache = init(problem)
    solve!(cache)
    reinit!(cache, problem.inputs)
    reused = solve!(cache)
    @test reused.soil_temperature == fresh.soil_temperature
    @test reused.canopy.leaf_temperature == fresh.canopy.leaf_temperature
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
    peak_sun = argmax(wet_output.global_radiation)
    @test all(dry_output.canopy.leaf_temperature[peak_sun, :] .> wet_output.canopy.leaf_temperature[peak_sun, :])
end
