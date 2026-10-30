using Microclimate
using Microclimate: AboveCanopyColumn, canopy_top_enhancement, canopy_top_flux_boundary, column_heat_resistance
using Unitful
using Test

include("helpers.jl")

boundary_layer_model = MoninObukhov()
canopy_height = 20.0u"m"
displacement_height = 14.0u"m"
reference_height = 45.0u"m"
ground_temperature = 298.0u"K"
reference_temperature = 300.0u"K"
ground_vapour_density = 0.018u"kg/m^3"
reference_vapour_density = 0.017u"kg/m^3"
previous_top_temperature = 299.0u"K"
previous_top_vapour_density = 0.0175u"kg/m^3"

column(; friction_velocity, displacement = displacement_height) = AboveCanopyColumn(; canopy_height,
    displacement_height = displacement, friction_velocity, inverse_obukhov_length = 0.0u"m^-1", enhancement = 1.0)

@testset "canopy_top_flux_boundary: windy conditions barely damped" begin
    result = canopy_top_flux_boundary(boundary_layer_model, column(; friction_velocity = 1.0u"m/s"), reference_height,
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
    run(disp) = canopy_top_flux_boundary(boundary_layer_model,
        column(; friction_velocity = 1.0e-4u"m/s", displacement = disp), reference_height,
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
    result = canopy_top_flux_boundary(boundary_layer_model, column(; friction_velocity = 1.0e-4u"m/s"), reference_height,
        0.0u"W/m^2", 0.0u"kg/m^2/s", reference_temperature, reference_vapour_density,
        reference_temperature, reference_vapour_density, reference_temperature, reference_vapour_density)

    @test ustrip(u"K", result.top_temperature) ≈ ustrip(u"K", reference_temperature) atol=1e-6
end

@testset "canopy_top_flux_boundary: vapour takes the floor's vapour resistance" begin
    boundary(floor_resistance) = canopy_top_flux_boundary(boundary_layer_model, column(; friction_velocity = 0.3u"m/s"),
        reference_height, 0.0u"W/m^2", 0.0u"kg/m^2/s", ground_temperature, ground_vapour_density,
        reference_temperature, reference_vapour_density, previous_top_temperature, previous_top_vapour_density;
        floor_resistance)
    heat_only = boundary((; heat = 50.0u"s/m", vapor = 0.0u"s/m"))
    both = boundary((; heat = 50.0u"s/m", vapor = 50.0u"s/m"))
    @test heat_only.top_temperature == both.top_temperature
    @test heat_only.top_vapour_density != both.top_vapour_density
end

@testset "roughness sublayer" begin
    enhancement(model) = canopy_top_enhancement(model, canopy_height, 5.0, 0.004u"m", boundary_layer_model.karman_constant)
    @test enhancement(NoRoughnessSublayer()) == 1.0
    @test enhancement(CanopyTopMixingEnhancement()) ≈ 3.0
    # Fades to 1 in slow flow over a short canopy
    fade(u) = Microclimate.turbulent_enhancement(CanopyTopMixingEnhancement(), 3.0, 0.05u"m", u, 1.5e-5u"m^2/s")
    @test fade(0.01u"m/s") < 1.1 && fade(5.0u"m/s") > 2.9
    # A canopy barely above the ground roughness closes no gap
    @test canopy_top_enhancement(CanopyTopMixingEnhancement(), 0.05u"m", 5.0, 0.02u"m",
        boundary_layer_model.karman_constant) == 1.0
    # No fade, even in still air
    @test Microclimate.turbulent_enhancement(CanopyTopMixingEnhancement(; transition_reynolds_number = 0.0), 3.0,
        0.05u"m", 0.0u"m/s", 1.5e-5u"m^2/s") == 3.0

    column(e) = AboveCanopyColumn(; canopy_height, displacement_height, friction_velocity = 0.5u"m/s",
        inverse_obukhov_length = -0.01u"m^-1", enhancement = e)
    R(e, z₁, z₂) = column_heat_resistance(column(e), boundary_layer_model, z₁, z₂)
    @test R(3.0, canopy_height, 30.0u"m") + R(3.0, 30.0u"m", reference_height) ≈ R(3.0, canopy_height, reference_height)
    @test R(3.0, canopy_height, reference_height) < R(1.0, canopy_height, reference_height)
    @test allocations(column_heat_resistance, column(3.0), boundary_layer_model, canopy_height, reference_height) == 0
end
