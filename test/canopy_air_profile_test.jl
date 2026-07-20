using Microclimate
using Unitful
using Test

canopy_height = 1.0u"m"
plant_area_index = 3.0
heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"
boundary_layer_model = MoninObukhov()
model = KTheoryAirProfile()
buffers = Microclimate.allocate_air_profile(model, canopy_height, plant_area_index, heights, 10, boundary_layer_model)
n_layers = length(buffers.air_temperature)

displacement_height = 0.65u"m"
friction_velocity = 0.3u"m/s"
atmospheric_pressure = 101325.0u"Pa"
zero_evaporation = zeros(typeof(0.0u"g/s"), n_layers)

@testset "no forcing: uniform boundary temperatures and RH give a uniform profile" begin
    T = 290.0u"K"
    result = Microclimate.canopy_air_profile!(buffers, model, boundary_layer_model;
        canopy_height, displacement_height, friction_velocity,
        canopy_top_air_temperature = T, canopy_top_relative_humidity = 0.5, ground_temperature = T,
        ground_relative_humidity = 0.5, sensible_heat_source = zeros(typeof(0.0u"W/m^2"), n_layers),
        evaporation_mass_flow = zero_evaporation, atmospheric_pressure,
    )
    @test all(t -> isapprox(ustrip(u"K", t), ustrip(u"K", T); atol=1e-8), result.air_temperature)
    @test all(rh -> isapprox(rh, 0.5; atol=1e-6), result.relative_humidity)
end

@testset "a mid-canopy heat source produces a local warm bump" begin
    source = zeros(typeof(0.0u"W/m^2"), n_layers)
    source[n_layers ÷ 2] = 50.0u"W/m^2"
    result = Microclimate.canopy_air_profile!(buffers, model, boundary_layer_model;
        canopy_height, displacement_height, friction_velocity,
        canopy_top_air_temperature = 290.0u"K", canopy_top_relative_humidity = 0.5,
        ground_temperature = 290.0u"K", ground_relative_humidity = 0.5,
        sensible_heat_source = source, evaporation_mass_flow = zero_evaporation, atmospheric_pressure,
    )
    peak = argmax(result.air_temperature)
    @test peak == n_layers ÷ 2
    @test result.air_temperature[peak] > 290.0u"K"
    # monotone decay away from the peak toward both boundaries
    @test issorted(result.air_temperature[1:peak])
    @test issorted(result.air_temperature[peak:end], rev=true)
end

@testset "a mid-canopy evaporative source produces a local humidity bump" begin
    source = zeros(typeof(0.0u"g/s"), n_layers)
    source[n_layers ÷ 2] = 1.0e-3u"g/s"
    result = Microclimate.canopy_air_profile!(buffers, model, boundary_layer_model;
        canopy_height, displacement_height, friction_velocity,
        canopy_top_air_temperature = 290.0u"K", canopy_top_relative_humidity = 0.5,
        ground_temperature = 290.0u"K", ground_relative_humidity = 0.5,
        sensible_heat_source = zeros(typeof(0.0u"W/m^2"), n_layers),
        evaporation_mass_flow = source, atmospheric_pressure,
    )
    peak = argmax(result.relative_humidity)
    @test peak == n_layers ÷ 2
    @test result.relative_humidity[peak] > 0.5
    # monotone decay away from the peak toward both boundaries
    @test issorted(result.relative_humidity[1:peak])
    @test issorted(result.relative_humidity[peak:end], rev=true)
end

@testset "stronger mixing (higher friction velocity) flattens the same forcing" begin
    source = zeros(typeof(0.0u"W/m^2"), n_layers)
    source[n_layers ÷ 2] = 50.0u"W/m^2"
    weak_mixing = Microclimate.canopy_air_profile!(buffers, model, boundary_layer_model;
        canopy_height, displacement_height, friction_velocity = 0.1u"m/s",
        canopy_top_air_temperature = 290.0u"K", canopy_top_relative_humidity = 0.5,
        ground_temperature = 290.0u"K", ground_relative_humidity = 0.5,
        sensible_heat_source = source, evaporation_mass_flow = zero_evaporation, atmospheric_pressure,
    )
    peak_weak = maximum(weak_mixing.air_temperature)
    strong_mixing = Microclimate.canopy_air_profile!(buffers, model, boundary_layer_model;
        canopy_height, displacement_height, friction_velocity = 1.0u"m/s",
        canopy_top_air_temperature = 290.0u"K", canopy_top_relative_humidity = 0.5,
        ground_temperature = 290.0u"K", ground_relative_humidity = 0.5,
        sensible_heat_source = source, evaporation_mass_flow = zero_evaporation, atmospheric_pressure,
    )
    peak_strong = maximum(strong_mixing.air_temperature)
    @test peak_strong < peak_weak
end

@testset "canopy_air_profile! allocation stays bounded" begin
    source = fill(10.0u"W/m^2", n_layers)
    evap_source = fill(1.0e-4u"g/s", n_layers)
    f() = Microclimate.canopy_air_profile!(buffers, model, boundary_layer_model;
        canopy_height, displacement_height, friction_velocity,
        canopy_top_air_temperature = 290.0u"K", canopy_top_relative_humidity = 0.5,
        ground_temperature = 288.0u"K", ground_relative_humidity = 0.5,
        sensible_heat_source = source, evaporation_mass_flow = evap_source, atmospheric_pressure,
    )
    f() # warm up
    @test (@allocated f()) < 5_000
end
