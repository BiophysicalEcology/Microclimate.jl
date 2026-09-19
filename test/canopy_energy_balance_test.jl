using Microclimate
using Unitful
using Test

body = Microclimate.leaf_body(0.05u"m", 0.02u"m", 1.0)
stomatal_conductance = LeafEvaporationParameters(;
    abaxial_vapour_conductance = 0.3u"mol/m^2/s", adaxial_vapour_conductance = 0.0u"mol/m^2/s",
    cuticular_conductance = 0.01u"mol/m^2/s",
)
air_temperature = 293.0u"K"
relative_humidity = 0.5
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

@testset "leaf temperature increases monotonically with absorbed radiation" begin
    for solver in (LinearizedLeafTemperature(), RootFindLeafTemperature())
        temperatures = [
            ustrip(u"K", Microclimate.leaf_temperature(solver, Rabs, air_temperature, relative_humidity,
                wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body))
            for Rabs in (50.0, 200.0, 400.0, 700.0) .* u"W/m^2"
        ]
        @test issorted(temperatures)
    end
end

@testset "iterated LinearizedLeafTemperature converges to RootFindLeafTemperature" begin
    Rabs = 300.0u"W/m^2"
    guess = air_temperature
    for _ in 1:8
        guess = Microclimate.leaf_temperature(LinearizedLeafTemperature(), Rabs, air_temperature, relative_humidity,
            wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body;
            leaf_temperature_guess=guess)
    end
    exact = Microclimate.leaf_temperature(RootFindLeafTemperature(), Rabs, air_temperature, relative_humidity,
        wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body)
    @test ustrip(u"K", guess) ≈ ustrip(u"K", exact) atol=1e-2
end

@testset "water stress warms the leaf (less evaporative cooling)" begin
    Rabs = 300.0u"W/m^2"
    wet = Microclimate.leaf_temperature(RootFindLeafTemperature(), Rabs, air_temperature, relative_humidity,
        wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, 0.0u"J/kg", body)
    dry = Microclimate.leaf_temperature(RootFindLeafTemperature(), Rabs, air_temperature, relative_humidity,
        wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, -2.0e6u"J/kg", body)
    @test dry > wet
end

@testset "leaf_temperature is allocation-light" begin
    solver = LinearizedLeafTemperature()
    f() = Microclimate.leaf_temperature(solver, 300.0u"W/m^2", air_temperature, relative_humidity,
        wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body)
    f() # warm up
    @test (@allocated f()) < 100
end

@testset "energy balance closes at the solved temperature" begin
    for solver in (LinearizedLeafTemperature(), RootFindLeafTemperature())
        Rabs = 300.0u"W/m^2"
        Tleaf = Microclimate.leaf_temperature(solver, Rabs, air_temperature, relative_humidity,
            wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body)
        out = Microclimate.leaf_heat_balance(Tleaf, Rabs, air_temperature, relative_humidity,
            wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body, 1.0u"m^2")
        # LinearizedLeafTemperature is a single step from air_temperature, so its
        # residual is only small when the true solution is close to air_temperature;
        # allow a looser tolerance for it than for the exact root-find.
        tol = solver isa RootFindLeafTemperature ? 2e-2u"W" : 10.0u"W"
        @test abs(out.net) < tol
    end
end

@testset "clamped temperature: net residual reroutes to sensible, balance still closes" begin
    # Mimics _canopy_picard_pass!'s clamp: evaluate away from the true
    # equilibrium, so out.net is a real nonzero residual, not solver noise.
    Rabs = 300.0u"W/m^2"
    Tclamped = air_temperature + 40.0u"K"
    out = Microclimate.leaf_heat_balance(Tclamped, Rabs, air_temperature, relative_humidity,
        wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body, 1.0u"m^2")
    @test abs(out.net) > 1.0u"W"  # confirms this scenario actually has a residual to reroute

    reallocated_sensible = out.conv.convection_flow + out.net
    closed_balance = reallocated_sensible + out.evap.evaporation_heat_flow + out.emitted_longwave
    @test ustrip(u"W", closed_balance) ≈ ustrip(u"W", Rabs * 1.0u"m^2") atol=1e-9
end
