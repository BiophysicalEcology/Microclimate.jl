using Microclimate
using Microclimate: AboveCanopyColumn, CanopyAirProfileForcing, GroundExchangeForcing, allocate_air_profile, canopy_heat_vapor_transport!,
    wet_air_properties
using Unitful
using Test

include("helpers.jl")

canopy_height = 1.0u"m"
plant_area_index = 3.0
heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"
boundary_layer_model = MoninObukhov()
site = Microclimate.example_site()
displacement_height = 0.65u"m"
n_layers = 10
T = 290.0u"K"

no_heat = zeros(typeof(0.0u"W/m^2"), n_layers)
no_evaporation = zeros(typeof(0.0u"g/s"), n_layers)
mid_layer = n_layers ÷ 2
mid_heat = [i == mid_layer ? 50.0u"W/m^2" : 0.0u"W/m^2" for i in 1:n_layers]
mid_evaporation = [i == mid_layer ? 1.0e-3u"g/s" : 0.0u"g/s" for i in 1:n_layers]

forcing(; sensible_heat_source = no_heat, evaporation_mass_flow = no_evaporation,
    friction_velocity = 0.3u"m/s", ground_temperature = T) = CanopyAirProfileForcing(;
    above_canopy_column = AboveCanopyColumn(; canopy_height, displacement_height,
        friction_velocity, inverse_obukhov_length = 0.0u"m^-1", enhancement = 1.0),
    wind_speed = fill(1.0u"m/s", n_layers),
    canopy_top_air_temperature = T, canopy_top_relative_humidity = 0.5,
    ground_temperature, ground_relative_humidity = 0.5,
    sensible_heat_source, evaporation_mass_flow,
    leaf_temperature = fill(T, n_layers), atmospheric_pressure = 101325.0u"Pa",
)

top_vapor_density = wet_air_properties(T, 0.5, 101325.0u"Pa").vapour_density

# Fresh buffers with the current air state isothermal at the boundary humidity
function initial_buffers(model)
    buffers = allocate_air_profile(model, canopy_height, displacement_height, plant_area_index, heights, n_layers, boundary_layer_model)
    fill!(buffers.air_temperature_prev, T)
    fill!(buffers.vapor_density_prev, top_vapor_density)
    return buffers
end

function profile(model; kw...)
    buffers = initial_buffers(model)
    canopy_heat_vapor_transport!(buffers, model, boundary_layer_model, forcing(; kw...); site)
    return buffers
end

# A single peak at (or next to) the source layer, everything else below it
is_local_bump(values) = abs(argmax(values) - mid_layer) <= 1 && count(==(maximum(values)), values) == 1

@testset "$(nameof(typeof(model)))" for model in (KTheoryAirProfile(), RaupachLTheoryAirProfile())
    @testset "no forcing: uniform boundaries give a uniform profile" begin
        result = profile(model)
        @test all(t -> isapprox(t, T; atol = 1e-6u"K"), result.air_temperature)
        @test all(v -> isapprox(v, top_vapor_density; rtol = 1e-6), result.vapor_density)
    end

    @testset "a mid-canopy heat source produces a local warm bump" begin
        result = profile(model; sensible_heat_source = mid_heat)
        @test maximum(result.air_temperature) > T
        @test is_local_bump(result.air_temperature)
    end

    @testset "a mid-canopy evaporative source produces a local humidity bump" begin
        result = profile(model; evaporation_mass_flow = mid_evaporation)
        @test maximum(result.vapor_density) > top_vapor_density
        @test is_local_bump(result.vapor_density)
    end

    @testset "stronger mixing flattens the same forcing" begin
        weak = maximum(profile(model; sensible_heat_source = mid_heat, friction_velocity = 0.1u"m/s").air_temperature)
        strong = maximum(profile(model; sensible_heat_source = mid_heat, friction_velocity = 1.0u"m/s").air_temperature)
        @test strong < weak
    end

    @testset "zero allocations: $(nameof(typeof(ground_coupling)))" for ground_coupling in
            (ProfileGroundResistance(), FloorTurbulence(), SakaguchiZeng2009())
        buffers = initial_buffers(model)
        hour_forcing = forcing(; sensible_heat_source = fill(10.0u"W/m^2", n_layers),
            evaporation_mass_flow = fill(1.0e-4u"g/s", n_layers), ground_temperature = 288.0u"K")
        @test allocations(canopy_heat_vapor_transport!, buffers, model, boundary_layer_model, hour_forcing;
            site, ground_coupling) == 0
    end
end

@testset "ground_resistance" begin
    ground_forcing(; friction_velocity = 0.3u"m/s", plant_area_index = 3.0, ground_temperature = T) = GroundExchangeForcing(;
        profile_heat_resistance = 100.0u"s/m", profile_vapor_resistance = 100.0u"s/m", layer_height = 0.1u"m",
        wind_speed = 0.1u"m/s", layer_air_temperature = T, ground_temperature, roughness_height = 0.01u"m",
        friction_velocity, plant_area_index, atmospheric_pressure = 101325.0u"Pa")
    ground_resistance(coupling, forcing) = Microclimate.ground_resistance(coupling, boundary_layer_model, forcing)
    @test ground_resistance(ProfileGroundResistance(), ground_forcing()) == (; heat = 100.0u"s/m", vapor = 100.0u"s/m")
    for coupling in (FloorTurbulence(), SakaguchiZeng2009())
        @test ground_resistance(coupling, ground_forcing(; friction_velocity = 1.0u"m/s")).heat <
              ground_resistance(coupling, ground_forcing()).heat
    end
    # Free convection keeps the floor coupled in calm air
    @test isfinite(ground_resistance(FloorTurbulence(), ground_forcing(; friction_velocity = 0.0u"m/s")).heat)
    # A warm floor under cooler air is more strongly coupled; a cold one less
    @test ground_resistance(FloorTurbulence(), ground_forcing(; ground_temperature = T + 5.0u"K")).heat <
          ground_resistance(FloorTurbulence(), ground_forcing()).heat <
          ground_resistance(FloorTurbulence(), ground_forcing(; ground_temperature = T - 5.0u"K")).heat
    # Dense canopy: Sakaguchi & Zeng's 1/(0.004 u★)
    @test ground_resistance(SakaguchiZeng2009(), ground_forcing(; plant_area_index = 50.0)).heat ≈ 1 / (0.004 * 0.3u"m/s")
end

@testset "KTheoryAirProfile: diffusion decays monotonically away from a source" begin
    result = profile(KTheoryAirProfile(); sensible_heat_source = mid_heat)
    peak = argmax(result.air_temperature)
    @test issorted(result.air_temperature[1:peak])
    @test issorted(result.air_temperature[peak:end]; rev = true)
end
