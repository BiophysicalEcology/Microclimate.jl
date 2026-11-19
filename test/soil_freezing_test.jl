using Microclimate, Unitful, Test

model = IsothermalPhaseTransition()
depths = [0.0, 2.5, 5.0]u"cm"
soil_moisture = [0.3, 0.3, 0.3]
buffers = Microclimate.allocate_phase_transition(model, length(depths))
water_mass = 2.5u"cm" * 1000.0u"kg/m" * 0.3
heat_per_kelvin = uconvert(u"J/K", water_mass * 4184.0u"J/kg/K")
max_latent_heat = uconvert(u"J", water_mass * Microclimate.LATENT_HEAT_FUSION)
celsius(t) = u"K"(t * u"°C")

# The first layer's new temperature and latent heat, with the layers below at 5°C
function first_layer(temperature, latent_heat)
    accumulated_latent_heat = [latent_heat, 0.0u"J", 0.0u"J"]
    temperatures = [celsius(temperature), celsius(temperature), celsius(5.0)]
    new_temperatures = Microclimate.phase_transition!(model, buffers; temperatures,
        temperatures_past=temperatures, accumulated_latent_heat, soil_moisture, depths)
    return (u"°C"(new_temperatures[1]), accumulated_latent_heat[1])
end

@testset "Freezing holds the layer at 0°C and stores the heat lost" begin
    (temperature, latent_heat) = first_layer(-2.0, 0.0u"J")
    @test temperature ≈ 0.0u"°C" atol=1e-9u"K"
    @test latent_heat ≈ 2.0u"K" * heat_per_kelvin
    (temperature, latent_heat) = first_layer(-2.0, latent_heat)
    @test temperature ≈ 0.0u"°C" atol=1e-9u"K"
    @test latent_heat ≈ 4.0u"K" * heat_per_kelvin
end

@testset "Once frozen the layer cools freely" begin
    (temperature, latent_heat) = first_layer(-2.0, max_latent_heat)
    @test temperature ≈ -2.0u"°C"
    @test latent_heat == max_latent_heat
end

@testset "Heat beyond the budget changes the temperature" begin
    (temperature, latent_heat) = first_layer(-2.0, max_latent_heat - 0.5u"K" * heat_per_kelvin)
    @test temperature ≈ -1.5u"°C"
    @test latent_heat ≈ max_latent_heat
    (temperature, latent_heat) = first_layer(2.0, 0.5u"K" * heat_per_kelvin)
    @test temperature ≈ 1.5u"°C"
    @test latent_heat == 0.0u"J"
end

@testset "Thawing holds the layer at 0°C until its ice melts" begin
    (temperature, latent_heat) = first_layer(1.0, max_latent_heat)
    @test temperature ≈ 0.0u"°C" atol=1e-9u"K"
    @test latent_heat ≈ max_latent_heat - 1.0u"K" * heat_per_kelvin
end

@testset "Unfrozen soil above 0°C is unchanged" begin
    (temperature, latent_heat) = first_layer(3.0, 0.0u"J")
    @test temperature ≈ 3.0u"°C"
    @test latent_heat == 0.0u"J"
end
