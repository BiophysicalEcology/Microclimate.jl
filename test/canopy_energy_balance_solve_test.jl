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

# Every residual within the tolerance for its kind
function residual_within_tolerance(buffers)
    (; residual, tolerances) = buffers.convergence
    return Microclimate.canopy_scaled_residual(residual, tolerances) <= 1
end

@testset "midday: residual solved, physically reasonable" begin
    model = canopy()
    result, buffers = solve_canopy(model, canopy_inputs(model))
    @test result.iterations < model.convergence_model.maxiters
    @test residual_within_tolerance(buffers)
    @test all(t -> 250.0u"K" < t < 340.0u"K", buffers.leaf.leaf_temperature)
    @test all(isfinite, buffers.air_profile.air_temperature)
    @test result.net_absorbed_shortwave > 0.0u"W/m^2"
    @test result.canopy_potential_transpiration >= sum(buffers.leaf.evaporation_mass_flux)
end

@testset "the heat leaving the canopy top is the heat the leaves and the ground give the air" begin
    model = canopy()
    result, buffers = solve_canopy(model, canopy_inputs(model))
    (; heat_conductance) = result.ground_air
    ground_sensible_heat_flux = heat_conductance * (290.0u"K" - last(buffers.air_profile.air_temperature))
    @test result.ground_sensible_heat_flux ≈ ground_sensible_heat_flux atol = heat_conductance * model.convergence_model.air_temperature_tolerance
    @test result.canopy_top_sensible_heat_flux ≈ sum(buffers.leaf.sensible_heat_flux) + ground_sensible_heat_flux atol =
        heat_conductance * model.convergence_model.air_temperature_tolerance
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

@testset "zero allocations" begin
    for plant_area_index in (plant_area_index, fill(plant_area_index / 10, 10), [0.6, 0.5, 0.45, 0.4, 0.35, 0.3, 0.2, 0.15, 0.1, 0.05])
        model = canopy(; plant_area_index)
        buffers = allocate(model)
        @test allocations(canopy_energy_balance!, buffers, model, boundary_layer_model, canopy_inputs(model); site) == 0
    end
end

@testset "AndersonCanopyConvergence failures and leafless layers" begin
    # One pass can't converge
    failing_model = Microclimate.example_multilayer_canopy(; canopy_height, plant_area_index,
        convergence_model = AndersonCanopyConvergence(; maxiters = 1))
    @test_throws r"did not converge" canopy_energy_balance!(allocate(failing_model), failing_model,
        boundary_layer_model, canopy_inputs(failing_model); site)
    # ... unless it takes the state with the smallest residual
    best_model = Microclimate.example_multilayer_canopy(; canopy_height, plant_area_index,
        convergence_model = AndersonCanopyConvergence(; maxiters = 1, on_failure = SmallestResidual()))
    best_buffers = allocate(best_model)
    best = @test_logs (:warn, r"did not converge") canopy_energy_balance!(best_buffers, best_model,
        boundary_layer_model, canopy_inputs(best_model); site)
    @test best.iterations == 1
    @test all(isfinite, best_buffers.leaf.leaf_temperature)
    # A layer without leaf area maps its leaf to the layer's air temperature
    gappy_model = Microclimate.example_multilayer_canopy(; canopy_height,
        plant_area_index = [0.6, 0.5, 0.0, 0.4, 0.35, 0.3, 0.2, 0.15, 0.1, 0.05],
        convergence_model = AndersonCanopyConvergence())
    gappy_inputs = canopy_inputs(gappy_model)
    gappy_result, gappy = solve_canopy(gappy_model, gappy_inputs)
    hour_conditions = Microclimate.canopy_hour_conditions(gappy, gappy_model, boundary_layer_model, gappy_inputs, gappy_result;
        site, vapour_pressure_equation = Microclimate.GoffGratch(), soil_hydraulic_model = nothing)
    (; state, step, next_state) = gappy.convergence
    (; ground_temperature, inverse_obukhov_length) = Microclimate.unpack_canopy_state!(gappy, state)
    gappy.leaf.leaf_temperature[3] += 1.0u"K"
    starting_air_temperature = copy(gappy.air_profile.current_air_temperature)
    Microclimate.pack_canopy_state!(state, gappy.leaf.leaf_temperature, starting_air_temperature,
        gappy.air_profile.current_vapor_density, ground_temperature, inverse_obukhov_length)
    Microclimate.canopy_fixed_point!(next_state, step, state, gappy, gappy_model, boundary_layer_model, gappy_inputs, hour_conditions;
        site, vapour_pressure_equation = Microclimate.GoffGratch(),
        ground_temperature_after_hour = (ground_temperature, _, _) -> ground_temperature)
    Microclimate.unpack_canopy_state!(gappy, next_state)
    @test gappy.leaf.leaf_temperature[3] ≈ starting_air_temperature[3]
end

@testset "AndersonCanopyConvergence solves the ground temperature with the canopy" begin
    model = canopy()
    buffers = allocate(model)
    calls = Ref(0)
    # A soil that moves halfway towards 290 K, whatever the canopy does
    function ground_temperature_after_hour(ground_temperature, ground_absorbed_shortwave, ground_absorbed_longwave, ground_air)
        calls[] += 1
        @test ground_absorbed_shortwave > 0.0u"W/m^2"
        @test ground_absorbed_longwave > 0.0u"W/m^2"
        @test ground_air.heat_conductance > 0.0u"W/m^2/K"
        return (ground_temperature + 290.0u"K") / 2
    end
    canopy_energy_balance!(buffers, model, boundary_layer_model, canopy_inputs(model; ground_temperature = 300.0u"K");
        site, ground_temperature_after_hour)
    @test calls[] > 1
    @test residual_within_tolerance(buffers)
    solved = Microclimate.unpack_canopy_state!(buffers, buffers.convergence.state)
    @test isapprox(solved.ground_temperature, 290.0u"K"; atol = 0.2u"K")
    # A soil that never settles is an error, not an accepted hour
    @test_throws r"did not converge" canopy_energy_balance!(allocate(model), model, boundary_layer_model,
        canopy_inputs(model); site, ground_temperature_after_hour = (ground_temperature, _, _, _) -> ground_temperature + 1.0u"K")
end

# An air profile that gives NaN on one call, the way a pass that overshoots does
struct OvershootingAirProfile{M} <: Microclimate.AbstractCanopyAirProfileModel
    air_profile_model::M
    overshooting_call::Int
    calls::Base.RefValue{Int}
end
OvershootingAirProfile(overshooting_call) = OvershootingAirProfile(KTheoryAirProfile(), overshooting_call, Ref(0))
Microclimate.allocate_air_profile(model::OvershootingAirProfile, arguments...) =
    Microclimate.allocate_air_profile(model.air_profile_model, arguments...)
Microclimate.canopy_ground_resistance(buffers, model::OvershootingAirProfile, arguments...; keywords...) =
    Microclimate.canopy_ground_resistance(buffers, model.air_profile_model, arguments...; keywords...)
function Microclimate.canopy_heat_vapor_transport!(buffers, model::OvershootingAirProfile, arguments...; keywords...)
    Microclimate.canopy_heat_vapor_transport!(buffers, model.air_profile_model, arguments...; keywords...)
    model.calls[] += 1
    model.calls[] == model.overshooting_call && fill!(buffers.air_temperature, NaN * u"K")
    return nothing
end

@testset "AndersonCanopyConvergence recovers from a pass that overshoots" begin
    overshooting_canopy(overshooting_call) = Microclimate.example_multilayer_canopy(; canopy_height, plant_area_index,
        air_profile_model = OvershootingAirProfile(overshooting_call))
    # Later passes go back towards the last good state and carry on
    model = overshooting_canopy(2)
    result, buffers = solve_canopy(model, canopy_inputs(model))
    @test model.air_profile_model.calls[] > 2
    @test residual_within_tolerance(buffers)
    clean_model = Microclimate.example_multilayer_canopy(; canopy_height, plant_area_index,
        air_profile_model = KTheoryAirProfile())
    _, clean_buffers = solve_canopy(clean_model, canopy_inputs(clean_model))
    @test all(isapprox.(buffers.leaf.leaf_temperature, clean_buffers.leaf.leaf_temperature; atol = 1.0u"K"))
    # The first pass has no good state to go back to
    first_pass_model = overshooting_canopy(1)
    @test_throws r"did not converge" solve_canopy(first_pass_model, canopy_inputs(first_pass_model))
end

@testset "AndersonCanopyConvergence needs a history" begin
    model = Microclimate.example_multilayer_canopy(; canopy_height, plant_area_index,
        convergence_model = AndersonCanopyConvergence(; history = 0))
    @test_throws ArgumentError allocate(model)
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
