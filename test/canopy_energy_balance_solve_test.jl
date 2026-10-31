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

canopy(; plant_area_index = plant_area_index) = Microclimate.example_multilayer_canopy(; canopy_height, plant_area_index)

allocate(model) = allocate_canopy(model, heights, boundary_layer_model;
    site, vapour_pressure_equation = Microclimate.GoffGratch(), soil_hydraulic_model = nothing)

function solve_canopy(model, inputs)
    buffers = allocate(model)
    result = canopy_energy_balance!(buffers, model, boundary_layer_model, inputs; site)
    return result, buffers
end

mean_leaf_temperature(buffers) = sum(buffers.leaf.leaf_temperature) / length(buffers.leaf.leaf_temperature)

# Every residual row within the tolerance for its kind
function residual_within_tolerance(buffers)
    (; residual, tolerances) = buffers.convergence
    return all(eachindex(residual)) do row
        abs(residual[row]) <= Microclimate.canopy_row_tolerance(tolerances, buffers, row)
    end
end

@testset "midday: residual solved, physically reasonable" begin
    model = canopy()
    result, buffers = solve_canopy(model, canopy_inputs(model))
    @test result.iterations < model.convergence_model.maxiters
    @test residual_within_tolerance(buffers)
    @test all(t -> 250.0u"K" < t < 340.0u"K", buffers.leaf.leaf_temperature)
    @test all(isfinite, buffers.air_profile.air_temperature)
    @test result.net_absorbed_shortwave > 0.0u"W/m^2"
    @test result.canopy_potential_transpiration >= sum(buffers.leaf.evaporation_mass_flow) / 1.0u"m^2"
end

@testset "new hourly inputs drive day -> night on the same buffers" begin
    model = canopy()
    buffers = allocate(model)
    @test canopy_energy_balance!(buffers, model, boundary_layer_model, canopy_inputs(model); site).net_absorbed_shortwave > 0.0u"W/m^2"

    night_inputs = canopy_inputs(model; zenith_angle = 100.0u"°", reference_temperature = 283.0u"K",
        direct_horizontal_irradiance = 0.0u"W/m^2", diffuse_horizontal_irradiance = 0.0u"W/m^2",
        ground_temperature = 280.0u"K")
    night = canopy_energy_balance!(buffers, model, boundary_layer_model, night_inputs; site)
    @test night.net_absorbed_shortwave == 0.0u"W/m^2"
    @test residual_within_tolerance(buffers)
    @test all(t -> 250.0u"K" < t < 310.0u"K", buffers.leaf.leaf_temperature)
end

@testset "water stress warms canopy leaves (less evaporative cooling)" begin
    model = canopy()
    leaf_temperature(leaf_water_potential) = mean_leaf_temperature(last(solve_canopy(model,
        canopy_inputs(model; reference_temperature = 298.0u"K", ground_temperature = 295.0u"K", leaf_water_potential))))
    @test leaf_temperature(-1.0e6u"J/kg") > leaf_temperature(0.0u"J/kg")
end

@testset "failed solve errors" begin
    model = Microclimate.example_multilayer_canopy(; canopy_height, plant_area_index,
        convergence_model = NonlinearCanopyConvergence(; maxiters = 1))
    buffers = allocate(model)
    @test_throws r"did not converge" canopy_energy_balance!(buffers, model, boundary_layer_model,
        canopy_inputs(model); site)
    @test !buffers.convergence.has_jacobian[]
end

@testset "zero allocations" begin
    for plant_area_index in (plant_area_index, fill(plant_area_index / 10, 10), [0.6, 0.5, 0.45, 0.4, 0.35, 0.3, 0.2, 0.15, 0.1, 0.05])
        model = canopy(; plant_area_index)
        buffers = allocate(model)
        @test allocations(canopy_energy_balance!, buffers, model, boundary_layer_model, canopy_inputs(model); site) == 0
    end
end

@testset "AndersonCanopyConvergence matches the Newton solve" begin
    newton_model = canopy()
    anderson_model = Microclimate.example_multilayer_canopy(; canopy_height, plant_area_index,
        convergence_model = AndersonCanopyConvergence())
    _, newton = solve_canopy(newton_model, canopy_inputs(newton_model))
    _, anderson = solve_canopy(anderson_model, canopy_inputs(anderson_model))
    @test residual_within_tolerance(anderson)
    @test maximum(abs, newton.leaf.leaf_temperature .- anderson.leaf.leaf_temperature) < 1.0e-4u"K"
    # One pass can't converge
    failing_model = Microclimate.example_multilayer_canopy(; canopy_height, plant_area_index,
        convergence_model = AndersonCanopyConvergence(; maxiters = 1))
    @test_throws r"did not converge" canopy_energy_balance!(allocate(failing_model), failing_model,
        boundary_layer_model, canopy_inputs(failing_model); site)
    @test allocations(canopy_energy_balance!, anderson, anderson_model, boundary_layer_model,
        canopy_inputs(anderson_model); site) == 0
    # A layer without leaf area maps its leaf to the layer's air temperature
    gappy_model = Microclimate.example_multilayer_canopy(; canopy_height,
        plant_area_index = [0.6, 0.5, 0.0, 0.4, 0.35, 0.3, 0.2, 0.15, 0.1, 0.05],
        convergence_model = AndersonCanopyConvergence())
    _, gappy = solve_canopy(gappy_model, canopy_inputs(gappy_model))
    (; state, residual, mapped, forcing) = gappy.convergence
    rows = Microclimate.canopy_state_rows(10)
    state[rows.leaf_temperature[3]] += 1.0
    Microclimate.canopy_fixed_point!(mapped, residual, state, gappy, gappy_model, boundary_layer_model, forcing[];
        site, vapour_pressure_equation = Microclimate.GoffGratch())
    @test mapped[rows.leaf_temperature[3]] ≈ state[rows.air_temperature[3]]
end

@testset "uniform vector plant_area_index matches the equivalent scalar" begin
    scalar_model = canopy()
    vector_model = canopy(; plant_area_index = fill(plant_area_index / 10, 10))
    scalar = last(solve_canopy(scalar_model, canopy_inputs(scalar_model))).leaf.leaf_temperature
    vector = last(solve_canopy(vector_model, canopy_inputs(vector_model))).leaf.leaf_temperature
    @test all(isapprox.(scalar, vector; atol = 1e-6u"K"))
end

@testset "lowest canopy layer must be above the ground roughness height" begin
    @test_throws ArgumentError allocate_canopy(canopy(), heights, boundary_layer_model;
        site = Microclimate.example_site(; roughness_height = 0.1u"m"),
        vapour_pressure_equation = Microclimate.GoffGratch(), soil_hydraulic_model = nothing)
end

@testset "canopy layer heights must be distinct and above the ground" begin
    @test_throws ArgumentError Microclimate.canopy_layer_heights([1.0, 0.5, 0.5] .* u"m", canopy_height, 3)
    @test_throws ArgumentError Microclimate.canopy_layer_heights([1.0, 0.5, 0.0] .* u"m", canopy_height, 3)
end
