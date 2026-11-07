using Microclimate
using Microclimate: canopy_top_flux_boundary
using Unitful
using Test

boundary_layer_model = MoninObukhov()
canopy_height = 20.0u"m"
displacement_height = 14.0u"m"
reference_height = 45.0u"m"
ground_temperature = 298.0u"K"
reference_temperature = 300.0u"K"
ground_vapor_density = 0.018u"kg/m^3"
reference_vapor_density = 0.017u"kg/m^3"
previous_top_temperature = 299.0u"K"
previous_top_vapor_density = 0.0175u"kg/m^3"

flux_boundary(friction_velocity; displacement_height = displacement_height) = canopy_top_flux_boundary(
    boundary_layer_model, canopy_height, displacement_height, reference_height, friction_velocity, 0.0u"m^-1",
    50.0u"W/m^2", 1.0e-5u"kg/m^2/s", ground_temperature, ground_vapor_density,
    reference_temperature, reference_vapor_density, previous_top_temperature, previous_top_vapor_density)

@testset "canopy_top_flux_boundary: windy conditions barely damped" begin
    result = flux_boundary(1.0u"m/s")

    @test isfinite(result.top_temperature)
    # fast equilibration (small residence time): close to reference, not stuck near previous_top_temperature
    @test abs(result.top_temperature - reference_temperature) < 2.0u"K"
end

@testset "canopy_top_flux_boundary: calm conditions damped toward previous value" begin
    # Near-zero friction velocity pins resistances at max_resistance, so
    # residence time scales with z_top alone (only displacement_height
    # varies below) -- bigger z_top should damp closer to previous_top.
    small_top_depth_result = flux_boundary(1.0e-4u"m/s"; displacement_height = 14.0u"m")  # 6 m
    big_top_depth_result = flux_boundary(1.0e-4u"m/s"; displacement_height = 2.0u"m")     # 18 m

    @test isfinite(small_top_depth_result.top_temperature)
    @test isfinite(big_top_depth_result.top_temperature)
    small_deviation = abs(small_top_depth_result.top_temperature - previous_top_temperature)
    big_deviation = abs(big_top_depth_result.top_temperature - previous_top_temperature)
    @test big_deviation < small_deviation
end

@testset "canopy_top_flux_boundary: no-op when already at the steady state" begin
    # A steady state with zero flux and previous == reference/ground-consistent point should stay put
    # regardless of damping strength (target == anchor).
    result = canopy_top_flux_boundary(boundary_layer_model, canopy_height, displacement_height, reference_height,
        1.0e-4u"m/s", 0.0u"m^-1", 0.0u"W/m^2", 0.0u"kg/m^2/s", reference_temperature, reference_vapor_density,
        reference_temperature, reference_vapor_density, reference_temperature, reference_vapor_density)

    @test result.top_temperature ≈ reference_temperature atol=1e-6u"K"
end
