using Microclimate
using Unitful
using Test

include("helpers.jl")

body = Microclimate.leaf_body(0.05u"m", 0.02u"m", 1.0)
stomatal_conductance = LeafEvaporationParameters(;
    abaxial_vapour_conductance = 0.3u"mol/m^2/s", adaxial_vapour_conductance = 0.0u"mol/m^2/s",
    cuticular_conductance = 0.01u"mol/m^2/s",
)
air_temperature = 293.0u"K"
relative_humidity = 0.5
air_vapor_density = Microclimate.wet_air_properties(air_temperature, relative_humidity, 101325.0u"Pa").vapour_density
wind_speed = 1.0u"m/s"
atmospheric_pressure = 101325.0u"Pa"
leaf_emissivity = 0.97
leaf_water_potential = 0.0u"J/kg"

@testset "leaf_body sizing" begin
    b = Microclimate.leaf_body(0.05u"m", 0.02u"m", 1.0)
    expected_radius = sqrt(0.05 * 0.02) * Microclimate.leaf_angle_width_factor(1.0)
    @test ustrip(u"m", b.geometry.length.b_semi_minor_skin) ≈ expected_radius rtol=1e-10
    @test b.geometry.length.b_semi_minor_skin == b.geometry.length.c_semi_minor_skin
end

@testset "leaf_heat_balance: net residual plus fluxes closes the balance away from equilibrium" begin
    # Far from the equilibrium leaf temperature, so out.net is a real nonzero residual.
    Rabs = 300.0u"W/m^2"
    Tleaf = air_temperature + 40.0u"K"
    out = Microclimate.leaf_heat_balance(Tleaf, Rabs, air_temperature, air_vapor_density,
        wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body, 1.0u"m^2")
    @test abs(out.net) > 1.0u"W"

    reallocated_sensible = out.conv.convection_flow + out.net
    closed_balance = reallocated_sensible + out.evap.evaporation_heat_flow + out.emitted_longwave
    @test ustrip(u"W", closed_balance) ≈ ustrip(u"W", Rabs * 1.0u"m^2") atol=1e-9
end

@testset "leaf_heat_balance is smooth through T_leaf = T_air" begin
    net(Tleaf, convection_model) = Microclimate.leaf_heat_balance(Tleaf, 300.0u"W/m^2", air_temperature, air_vapor_density,
        0.1u"m/s", atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body, 1.0u"m^2",
        convection_model).net
    # Slopes just below and just above T_air
    δ = 1.0e-3u"K"
    slopes(convection_model) = (net(air_temperature - δ, convection_model) - net(air_temperature - 2δ, convection_model),
        net(air_temperature + 2δ, convection_model) - net(air_temperature + δ, convection_model)) ./ δ
    below, above = slopes(ElaborateLeafConvection())
    @test below ≈ above rtol = 0.01
    # The exact |ΔT| has a cusp there
    exact_below, exact_above = slopes(ElaborateLeafConvection(; temperature_difference_model = ExactTemperatureDifference()))
    @test abs(exact_above - exact_below) > 100 * abs(above - below)
end
