using Microclimate
using Microclimate: CanopyAirProfileForcing, allocate_air_profile, canopy_air_profile!
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
    friction_velocity, obukhov_length = Inf * u"m", wind_speed = fill(1.0u"m/s", n_layers),
    canopy_top_air_temperature = T, canopy_top_relative_humidity = 0.5,
    ground_temperature, ground_relative_humidity = 0.5,
    sensible_heat_source, evaporation_mass_flow,
    leaf_temperature = fill(T, n_layers), atmospheric_pressure = 101325.0u"Pa",
)

# Fresh buffers in the isothermal, vapour-free initial state (the Raupach
# solve starts from the previous air state, and reseeds humidity from the
# boundary conditions when vapour density is zero)
function initial_buffers(model)
    buffers = allocate_air_profile(model, canopy_height, displacement_height, plant_area_index, heights, n_layers, boundary_layer_model)
    fill!(buffers.air_temperature, T)
    fill!(buffers.vapor_density, 0.0u"kg/m^3")
    return buffers
end

profile(model; kw...) = canopy_air_profile!(initial_buffers(model), model, boundary_layer_model, forcing(; kw...); site)

# A single peak at (or next to) the source layer, everything else below it
is_local_bump(values) = abs(argmax(values) - mid_layer) <= 1 && count(==(maximum(values)), values) == 1

@testset "$(nameof(typeof(model)))" for model in (KTheoryAirProfile(), RaupachLTheoryAirProfile())
    @testset "no forcing: uniform boundaries give a uniform profile" begin
        result = profile(model)
        @test all(t -> isapprox(t, T; atol = 1e-6u"K"), result.air_temperature)
        @test all(rh -> isapprox(rh, 0.5; atol = 1e-6), result.relative_humidity)
    end

    @testset "a mid-canopy heat source produces a local warm bump" begin
        result = profile(model; sensible_heat_source = mid_heat)
        @test maximum(result.air_temperature) > T
        @test is_local_bump(result.air_temperature)
    end

    @testset "a mid-canopy evaporative source produces a local humidity bump" begin
        result = profile(model; evaporation_mass_flow = mid_evaporation)
        @test maximum(result.relative_humidity) > 0.5
        @test is_local_bump(result.relative_humidity)
    end

    @testset "stronger mixing flattens the same forcing" begin
        weak = maximum(profile(model; sensible_heat_source = mid_heat, friction_velocity = 0.1u"m/s").air_temperature)
        strong = maximum(profile(model; sensible_heat_source = mid_heat, friction_velocity = 1.0u"m/s").air_temperature)
        @test strong < weak
    end

    @testset "zero allocations" begin
        buffers = initial_buffers(model)
        hour_forcing = forcing(; sensible_heat_source = fill(10.0u"W/m^2", n_layers),
            evaporation_mass_flow = fill(1.0e-4u"g/s", n_layers), ground_temperature = 288.0u"K")
        @test allocations(canopy_air_profile!, buffers, model, boundary_layer_model, hour_forcing; site) == 0
    end
end

@testset "KTheoryAirProfile: diffusion decays monotonically away from a source" begin
    result = profile(KTheoryAirProfile(); sensible_heat_source = mid_heat)
    peak = argmax(result.air_temperature)
    @test issorted(result.air_temperature[1:peak])
    @test issorted(result.air_temperature[peak:end]; rev = true)
end
