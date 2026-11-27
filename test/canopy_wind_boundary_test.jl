using Microclimate
using Microclimate: canopy_top_flux_boundary, calc_air_volumetric_heat_capacity
using Unitful
using Test

boundary_layer_model = MoninObukhov()
canopy_height = 20.0u"m"
displacement_height = 14.0u"m"
reference_height = 45.0u"m"
reference_temperature = 300.0u"K"
reference_vapor_density = 0.017u"kg/m^3"
previous_top_temperature = 299.0u"K"
previous_top_vapor_density = 0.0175u"kg/m^3"
air_volumetric_heat_capacity = calc_air_volumetric_heat_capacity(reference_temperature)

flux_boundary(friction_velocity; displacement_height = displacement_height, sensible_heat_flux = 50.0u"W/m^2",
    vapor_flux = 1.0e-5u"kg/m^2/s", previous_temperature = previous_top_temperature,
    previous_vapor_density = previous_top_vapor_density) = canopy_top_flux_boundary(
    boundary_layer_model, canopy_height, displacement_height, reference_height, friction_velocity, 0.0u"m^-1",
    sensible_heat_flux, vapor_flux, air_volumetric_heat_capacity,
    reference_temperature, reference_vapor_density, previous_temperature, previous_vapor_density)

@testset "canopy_top_flux_boundary: windy air carries the whole flux to the reference height" begin
    result = flux_boundary(1.0u"m/s")
    @test result.sensible_flux ≈ 50.0u"W/m^2" rtol = 1e-6
    @test result.top_temperature > reference_temperature
    @test result.top_vapor_density > reference_vapor_density
end

@testset "canopy_top_flux_boundary: calm air stays nearer the previous value under a deeper canopy top" begin
    shallow = flux_boundary(1.0e-4u"m/s"; displacement_height = 14.0u"m")
    deep = flux_boundary(1.0e-4u"m/s"; displacement_height = 2.0u"m")
    @test abs(deep.top_temperature - previous_top_temperature) < abs(shallow.top_temperature - previous_top_temperature)
end

@testset "canopy_top_flux_boundary: no flux at the reference values stays there" begin
    result = flux_boundary(1.0e-4u"m/s"; sensible_heat_flux = 0.0u"W/m^2", vapor_flux = 0.0u"kg/m^2/s",
        previous_temperature = reference_temperature, previous_vapor_density = reference_vapor_density)
    @test result.top_temperature ≈ reference_temperature atol = 1e-6u"K"
    @test result.top_vapor_density ≈ reference_vapor_density
end
