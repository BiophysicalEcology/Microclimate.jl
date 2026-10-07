using Microclimate
using Microclimate: canopy_top_flux_boundary
using Unitful
using Test

boundary_layer_model = MoninObukhov()
canopy_height = 20.0u"m"
displacement_height = 14.0u"m"
roughness_length = 2.0u"m"
reference_height = 45.0u"m"
ground_temperature = 298.0u"K"
reference_temperature = 300.0u"K"
ground_vapour_density = 0.018u"kg/m^3"
reference_vapour_density = 0.017u"kg/m^3"
previous_top_temperature = 299.0u"K"
previous_top_vapour_density = 0.0175u"kg/m^3"

@testset "canopy_top_flux_boundary: windy conditions barely damped" begin
    result = canopy_top_flux_boundary(boundary_layer_model, canopy_height, displacement_height, roughness_length,
        reference_height, 1.0u"m/s", Inf * u"m",
        50.0u"W/m^2", 1.0e-5u"kg/m^2/s", ground_temperature, ground_vapour_density,
        reference_temperature, reference_vapour_density, previous_top_temperature, previous_top_vapour_density)

    @test isfinite(ustrip(u"K", result.top_temperature))
    # fast equilibration (small residence time): close to reference, not stuck near previous_top_temperature
    @test abs(ustrip(u"K", result.top_temperature - reference_temperature)) < 2.0
end

@testset "canopy_top_flux_boundary: calm conditions damped toward previous value" begin
    # Near-zero friction velocity pins resistances at max_resistance, so
    # residence time scales with z_top alone (only displacement_height
    # varies below) -- bigger z_top should damp closer to previous_top.
    run(disp) = canopy_top_flux_boundary(boundary_layer_model, canopy_height, disp, roughness_length,
        reference_height, 1.0e-4u"m/s", Inf * u"m",
        50.0u"W/m^2", 1.0e-5u"kg/m^2/s", ground_temperature, ground_vapour_density,
        reference_temperature, reference_vapour_density, previous_top_temperature, previous_top_vapour_density)

    small_ztop_result = run(14.0u"m")  # z_top = 6m
    big_ztop_result = run(2.0u"m")     # z_top = 18m

    @test isfinite(ustrip(u"K", small_ztop_result.top_temperature))
    @test isfinite(ustrip(u"K", big_ztop_result.top_temperature))
    small_deviation = abs(ustrip(u"K", small_ztop_result.top_temperature - previous_top_temperature))
    big_deviation = abs(ustrip(u"K", big_ztop_result.top_temperature - previous_top_temperature))
    @test big_deviation < small_deviation
end

@testset "canopy_top_flux_boundary: no-op when already at the steady state" begin
    # A steady state with zero flux and previous == reference/ground-consistent point should stay put
    # regardless of damping strength (target == anchor).
    result = canopy_top_flux_boundary(boundary_layer_model, canopy_height, displacement_height, roughness_length,
        reference_height, 1.0e-4u"m/s", Inf * u"m",
        0.0u"W/m^2", 0.0u"kg/m^2/s", reference_temperature, reference_vapour_density,
        reference_temperature, reference_vapour_density, reference_temperature, reference_vapour_density)

    @test ustrip(u"K", result.top_temperature) ≈ ustrip(u"K", reference_temperature) atol=1e-6
end
