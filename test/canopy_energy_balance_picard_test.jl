using Microclimate
using Microclimate: allocate_canopy, canopy_energy_balance!
using Unitful
using Test

include("helpers.jl")

canopy_height = 1.0u"m"
plant_area_index = 3.0
heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"
boundary_layer_model = MoninObukhov()
site = Microclimate.example_site()

canopy(; plant_area_index = plant_area_index,
    convergence = IterationToleranceConvergence(; tolerance = 0.01u"K", max_iterations_per_day = 15)) =
    Microclimate.example_multilayer_canopy(; canopy_height, plant_area_index,
        convergence_model = PicardCanopyConvergence(; convergence))

function solve_canopy(model, inputs)
    buffers = allocate_canopy(model, heights, boundary_layer_model;
        site = Microclimate.example_site(), vapour_pressure_equation = Microclimate.GoffGratch(), soil_hydraulic_model = nothing)
    result = canopy_energy_balance!(buffers, model, boundary_layer_model, inputs; site)
    return result, buffers
end

mean_leaf_temperature(buffers) = sum(buffers.leaf.leaf_temperature) / length(buffers.leaf.leaf_temperature)

@testset "midday: converges, physically reasonable" begin
    model = canopy()
    result, buffers = solve_canopy(model, canopy_inputs(model))
    @test result.iterations < Microclimate.max_iterations(model.convergence_model.convergence)
    @test all(t -> 250.0u"K" < t < 340.0u"K", buffers.leaf.leaf_temperature)
    @test all(isfinite, buffers.air_profile.air_temperature)
    @test result.net_absorbed_shortwave > 0.0u"W/m^2"
    @test result.canopy_potential_transpiration >= sum(buffers.leaf.evaporation_mass_flow) / 1.0u"m^2"
end

@testset "new hourly inputs drive day -> night on the same buffers" begin
    model = canopy()
    buffers = allocate_canopy(model, heights, boundary_layer_model;
        site = Microclimate.example_site(), vapour_pressure_equation = Microclimate.GoffGratch(), soil_hydraulic_model = nothing)
    @test canopy_energy_balance!(buffers, model, boundary_layer_model, canopy_inputs(model); site).net_absorbed_shortwave > 0.0u"W/m^2"

    night_inputs = canopy_inputs(model; zenith_angle = 100.0u"°", reference_temperature = 283.0u"K",
        direct_horizontal_irradiance = 0.0u"W/m^2", diffuse_horizontal_irradiance = 0.0u"W/m^2",
        ground_temperature = 280.0u"K")
    night = canopy_energy_balance!(buffers, model, boundary_layer_model, night_inputs; site)
    @test night.net_absorbed_shortwave == 0.0u"W/m^2"
    @test all(t -> 250.0u"K" < t < 310.0u"K", buffers.leaf.leaf_temperature)
end

@testset "water stress warms canopy leaves (less evaporative cooling)" begin
    model = canopy(; convergence = FixedIterationConvergence(20))
    leaf_temperature(leaf_water_potential) = mean_leaf_temperature(last(solve_canopy(model,
        canopy_inputs(model; reference_temperature = 298.0u"K", ground_temperature = 295.0u"K", leaf_water_potential))))
    @test leaf_temperature(-1.0e6u"J/kg") > leaf_temperature(0.0u"J/kg")
end

@testset "FixedIterationConvergence runs exactly its iteration count" begin
    model = canopy(; convergence = FixedIterationConvergence(3))
    @test first(solve_canopy(model, canopy_inputs(model))).iterations == 3
end

@testset "zero allocations" begin
    for plant_area_index in (plant_area_index, fill(plant_area_index / 10, 10), [0.6, 0.5, 0.45, 0.4, 0.35, 0.3, 0.2, 0.15, 0.1, 0.05])
        model = canopy(; plant_area_index, convergence = FixedIterationConvergence(3))
        buffers = allocate_canopy(model, heights, boundary_layer_model;
        site = Microclimate.example_site(), vapour_pressure_equation = Microclimate.GoffGratch(), soil_hydraulic_model = nothing)
        @test allocations(canopy_energy_balance!, buffers, model, boundary_layer_model, canopy_inputs(model); site) == 0
    end
end

@testset "uniform vector plant_area_index matches the equivalent scalar" begin
    scalar_model = canopy(; convergence = FixedIterationConvergence(3))
    vector_model = canopy(; plant_area_index = fill(plant_area_index / 10, 10), convergence = FixedIterationConvergence(3))
    scalar = last(solve_canopy(scalar_model, canopy_inputs(scalar_model))).leaf.leaf_temperature
    vector = last(solve_canopy(vector_model, canopy_inputs(vector_model))).leaf.leaf_temperature
    @test all(isapprox.(scalar, vector; atol = 1e-2u"K"))
end

@testset "reset_canopy_scratch! clears Aitken warm-start state, keeps omega" begin
    model = canopy()
    _, buffers = solve_canopy(model, canopy_inputs(model))
    @test buffers.leaf.aitken_have_prev[]
    @test buffers.air_profile.aitken_have_prev[]

    leaf_omega, air_omega = buffers.leaf.aitken_omega[], buffers.air_profile.aitken_omega[]
    Microclimate.reset_canopy_scratch!(model, buffers)
    @test !buffers.leaf.aitken_have_prev[]
    @test !buffers.air_profile.aitken_have_prev[]
    @test !buffers.air_profile.aitken_have_prev_latent[]
    @test buffers.leaf.aitken_omega[] == leaf_omega
    @test buffers.air_profile.aitken_omega[] == air_omega
end
