using Microclimate
using Microclimate: litter_shortwave, litter_longwave, litter_transfer, litter_wind_speed, litter_temperature,
    litter_net_heat_flux, soil_surface_heat_flux, soil_vapor_exchange, ground_vapor_source, ground_surface, ground_reflectance,
    litter_element_humidity, litter_maximum_water_content, litter_water_exchange, litter_water_step, LitterVapor,
    litter_element_evaporation, pore_vapor_density,
    canopy_energy_balance!, GroundConditions, GroundAirExchange, ReferenceAirExchange, σ
using Unitful
using Test

include("helpers.jl")

litter = LitterLayer(; dry_mass = 1.5u"kg/m^2", thickness = 2.5u"cm", cover_fraction = 0.9, albedo = 0.3)
soil(; emissivity = 0.95) = ground(; albedo = 0.2, emissivity, wetness = 0.1)

# Under a canopy: its lowest layer's air and conductances
canopy_air = GroundAirExchange(; air_temperature = 293.0u"K", air_relative_humidity = 0.5, wind_speed = 0.3u"m/s",
    height = 0.5u"m", heat_conductance = 10.0u"W/m^2/K", vapor_conductance = 0.01u"m/s")
reference_air = ReferenceAirExchange(; boundary_layer_model = MoninObukhov(), air_temperature = 293.0u"K",
    air_relative_humidity = 0.5, wind_speed = 2.0u"m/s", height = 2.0u"m", roughness_length = 0.004u"m",
    zenith_angle = 30.0u"°", obukhov_length_prev = Ref(-0.3u"m"))
ground_conditions(air; downward_longwave = 350.0u"W/m^2", longwave_return_fraction = 0.0, exposed_fraction = 1.0) =
    GroundConditions(; direct_shortwave = 400.0u"W/m^2", diffuse_shortwave = 100.0u"W/m^2",
        sun_elevation_sine = cosd(30.0), downward_longwave, longwave_return_fraction, exposed_fraction,
        atmospheric_pressure = 101325.0u"Pa", vapour_pressure_equation = GoffGratch(), air)
conditions = ground_conditions(canopy_air)
heat_flux(layer, conditions, T) = soil_surface_heat_flux(layer, soil(), conditions, BulkTransferEvaporation(), T)

@testset "shortwave: absorbed by litter and soil plus reflected is the incident" begin
    for cover_fraction in (0.0, 0.5, 0.95), sun_elevation_sine in (0.2, 0.9), radiation in (GapLitterRadiation(), ShawLitterRadiation())
        layer = LitterLayer(; dry_mass = 1.5u"kg/m^2", thickness = 2.5u"cm", cover_fraction, radiation)
        result = litter_shortwave(layer, 400.0u"W/m^2", 100.0u"W/m^2", sun_elevation_sine, 0.2)
        @test result.litter_absorbed + result.soil_absorbed + result.reflected ≈ 500.0u"W/m^2"
        @test result.litter_absorbed >= 0.0u"W/m^2" && result.soil_absorbed >= 0.0u"W/m^2"
    end
    sparse = litter_shortwave(LitterLayer(; dry_mass = 1.5u"kg/m^2", thickness = 2.5u"cm", cover_fraction = 0.2), 400.0u"W/m^2", 100.0u"W/m^2", 0.9, 0.2)
    dense = litter_shortwave(LitterLayer(; dry_mass = 1.5u"kg/m^2", thickness = 2.5u"cm", cover_fraction = 0.9), 400.0u"W/m^2", 100.0u"W/m^2", 0.9, 0.2)
    @test dense.soil_absorbed < sparse.soil_absorbed
end

@testset "longwave: no net exchange in a black isothermal enclosure, and energy closes" begin
    T = 290.0u"K"
    isothermal = litter_longwave(litter, soil(; emissivity = 1.0), ground_conditions(canopy_air; downward_longwave = σ * T^4), T, T)
    @test abs(isothermal.litter_net) < 1.0e-9u"W/m^2"
    @test abs(isothermal.soil_net) < 1.0e-9u"W/m^2"
    result = litter_longwave(litter, soil(), conditions, 300.0u"K", 290.0u"K")
    @test 350.0u"W/m^2" - result.upward ≈ result.litter_net + result.soil_net
    # Part of what leaves the ground comes back down
    returning = litter_longwave(litter, soil(), ground_conditions(canopy_air; longwave_return_fraction = 0.3), 300.0u"K", 290.0u"K")
    downward = 350.0u"W/m^2" + 0.3 * returning.upward
    @test downward - returning.upward ≈ returning.litter_net + returning.soil_net
end

@testset "the ground surface: emission and reflectance give the upward longwave" begin
    return_conditions = ground_conditions(canopy_air; longwave_return_fraction = 0.3)
    T_soil = 290.0u"K"
    soil_with_litter = ground(; albedo = 0.2, wetness = 0.1, litter)
    surface = ground_surface(soil_with_litter, return_conditions, T_soil)
    T_litter = litter_temperature(litter, soil_with_litter, return_conditions,
        litter_shortwave(litter, 400.0u"W/m^2", 100.0u"W/m^2", cosd(30.0), 0.2).litter_absorbed, T_soil).temperature
    (; upward) = litter_longwave(litter, soil_with_litter, return_conditions, T_litter, T_soil)
    downward = 350.0u"W/m^2" + 0.3 * upward
    @test upward ≈ surface.emitted_longwave + surface.longwave_reflectance * downward
    (; direct, diffuse) = ground_reflectance(soil_with_litter, cosd(30.0))
    @test 0.0 < direct < 1.0 && 0.0 < diffuse < 1.0
    # Almost no litter is bare soil
    almost_none = ground(; albedo = 0.2, wetness = 0.1,
        litter = LitterLayer(; dry_mass = 1.5u"kg/m^2", thickness = 2.5u"cm", cover_fraction = 1.0e-6))
    bare = ground_surface(soil(), return_conditions, T_soil)
    nearly_bare = ground_surface(almost_none, return_conditions, T_soil)
    @test nearly_bare.temperature ≈ bare.temperature rtol = 1e-5
    @test nearly_bare.emitted_longwave ≈ bare.emitted_longwave rtol = 1e-4
    @test nearly_bare.reflectance.direct ≈ bare.reflectance rtol = 1e-4
    @test nearly_bare.reflectance.diffuse ≈ bare.reflectance rtol = 1e-4
    @test nearly_bare.longwave_reflectance ≈ bare.longwave_reflectance atol = 1e-6
end

@testset "conductance: still-air convection plus material conduction" begin
    # Bulk density 60 kg/m³ fills 60/170 of the layer with material; the rest is air
    material = 60.0 / 170.0
    wind = 0.3 * exp(-0.1 * (60.0 - 20.0) / 4)
    @test ustrip(u"m/s", litter_wind_speed(litter, canopy_air)) ≈ wind
    (; heat_conductance, vapor_resistance) = litter_transfer(litter, 0.0, litter_wind_speed(litter, canopy_air), 293.15u"K")
    conductivity = 0.025 * (1 + 0.007 * 20.0) * (1 + 4.0 * wind) * (1 - material) + 0.05 * material
    @test ustrip(u"W/m^2/K", heat_conductance) ≈ conductivity / 0.025 rtol = 1e-3
    @test vapor_resistance > 0.0u"s/m"
    h_air = 0.01u"m/s"
    gap = 1 - litter.cover_fraction
    h_soil = (1 - gap) / (1 / h_air + vapor_resistance) + gap * h_air
    # Through dry litter the soil's vapor crosses the whole layer in series with the air
    dry_vapor = LitterVapor(; covered_fraction = 1 - gap, half_resistance = vapor_resistance / 2,
        element_conductance = 0.0u"m/s", element_vapor_density = 0.0u"kg/m^3", soil_wetness = 1.0)
    exchange = soil_vapor_exchange(dry_vapor, h_air, 293.0u"K", 0.5, 101325.0u"Pa", GoffGratch())
    @test exchange.mass_transfer_coefficient ≈ h_soil
    @test exchange.relative_humidity ≈ 0.5
    # The ground's source gives the soil's flux across the air's resistance alone
    air_density = 0.008u"kg/m^3"
    source = ground_vapor_source(dry_vapor, h_air, 0.012u"kg/m^3", air_density)
    @test h_air * (source - air_density) ≈ h_soil * (0.012u"kg/m^3" - air_density)
end

@testset "litter temperature balances its fluxes, and sits above the soil in sun" begin
    T_soil = 295.0u"K"
    for air in (canopy_air, reference_air)
        air_conditions = ground_conditions(air)
        absorbed = litter_shortwave(litter, air_conditions.direct_shortwave, air_conditions.diffuse_shortwave,
            air_conditions.sun_elevation_sine, 0.2).litter_absorbed
        (; temperature) = litter_temperature(litter, soil(), air_conditions, absorbed, T_soil)
        wind_speed = litter_wind_speed(litter, air)
        flux(T) = litter_net_heat_flux(litter, soil(), air_conditions, absorbed, wind_speed, T, T_soil)
        @test flux(temperature - litter.tolerance) > 0.0u"W/m^2" > flux(temperature + litter.tolerance)
        @test temperature > T_soil
    end
end

@testset "litter insulates the soil surface from midday sun" begin
    T_soil = 295.0u"K"
    @test heat_flux(litter, conditions, T_soil) < heat_flux(NoLitter(), conditions, T_soil)
end

@testset "soil between litter elements acts as bare ground" begin
    for air in (canopy_air, reference_air), exposed_fraction in (1.0, 0.6)
        air_conditions = ground_conditions(air; exposed_fraction)
        T_soil = 295.0u"K"
        bare = heat_flux(NoLitter(), air_conditions, T_soil)
        almost_none = LitterLayer(; dry_mass = 1.5u"kg/m^2", thickness = 2.5u"cm", cover_fraction = 1.0e-6)
        @test heat_flux(almost_none, air_conditions, T_soil) ≈ bare rtol = 1e-4
        full = LitterLayer(; dry_mass = 1.5u"kg/m^2", thickness = 2.5u"cm", cover_fraction = 0.999)
        for T in (295.0u"K", 285.0u"K")
            covered = heat_flux(full, air_conditions, T)
            open = heat_flux(NoLitter(), air_conditions, T)
            half = heat_flux(LitterLayer(; dry_mass = 1.5u"kg/m^2", thickness = 2.5u"cm", cover_fraction = 0.5), air_conditions, T)
            @test min(covered, open) < half < max(covered, open)
        end
    end
end

@testset "the litter solves without allocating" begin
    for air in (canopy_air, reference_air)
        @test allocations(soil_surface_heat_flux, litter, soil(), ground_conditions(air), BulkTransferEvaporation(), 295.0u"K") == 0
    end
    soil_with_litter = ground(; litter)
    @test allocations(ground_surface, soil_with_litter, ground_conditions(canopy_air; longwave_return_fraction = 0.3), 295.0u"K") == 0
end

@testset "a canopy over litter solves" begin
    model = example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 2.0)
    site = example_site()
    buffers = Microclimate.allocate_canopy(model, vcat(0.1:0.1:1.0, [2.0]) .* u"m", MoninObukhov(); ground_roughness_length = 0.004u"m",
        vapour_pressure_equation = GoffGratch(), soil_hydraulic_model = nothing)
    result = canopy_energy_balance!(buffers, model, MoninObukhov(), canopy_inputs(model; litter); site)
    @test all(isfinite, ustrip.(u"K", buffers.leaf.leaf_temperature))
    @test result.ground_downward_longwave > 0.0u"W/m^2"
end

with_litter(problem, litter) = MicroProblem(Microclimate.setproperties(problem.model,
    (; surface = SurfaceModel(; litter_model = litter))), problem.inputs; days = problem.days, time_mode = problem.time_mode)
surface_range(result) = maximum(result.soil_temperature[:, 1]) - minimum(result.soil_temperature[:, 1])

@testset "a litter-covered simulation runs and differs from bare soil" begin
    bare_problem = example_microclimate_problem(; heights = [0.01, 2.0]u"m")
    covered = solve(with_litter(bare_problem, litter))
    bare = solve(bare_problem)
    @test all(isfinite, ustrip.(covered.soil_temperature))
    # Litter damps the soil surface's daily range
    @test surface_range(covered) < surface_range(bare)
end

@testset "litter under a canopy, with the soil solved or predicted" begin
    heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"
    for ground_temperature_model in (IntegratedGroundTemperature(), PredictedGroundTemperature())
        canopy_model = example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 2.0,
            convergence_model = AndersonCanopyConvergence(; ground_temperature_model))
        bare_problem = example_microclimate_problem(; heights, canopy_model)
        covered = solve(with_litter(bare_problem, litter))
        @test all(t -> 200.0u"K" < t < 340.0u"K", covered.soil_temperature)
        @test surface_range(covered) < surface_range(solve(bare_problem))
    end
end

wet = WetLitter(; element_vapor_resistance = 5000.0u"s/m")
wet_litter = LitterLayer(; dry_mass = 1.5u"kg/m^2", thickness = 2.5u"cm", cover_fraction = 0.9, albedo = 0.3, moisture = wet)
wet_soil(water_content) = ground(; albedo = 0.2, wetness = 0.1,
    litter = wet_litter, litter_water_content = water_content)

@testset "wet litter: its elements' humidity rises with water content to the maximum it holds" begin
    T = 293.0u"K"
    maximum_content = litter_maximum_water_content(wet, T)
    @test litter_element_humidity(wet, maximum_content, T) ≈ wet.maximum_relative_humidity
    @test litter_element_humidity(wet, 0.0, T) == 0.0
    contents = (0.01, 0.1, 0.5, maximum_content)
    @test issorted(map(w -> litter_element_humidity(wet, w, T), contents))
    # Myrold et al. (1981): wheat straw holds a few times its dry mass at saturation
    @test 1.0 < maximum_content < 5.0
end

@testset "wet litter: evaporation cools it, and its balance closes" begin
    T_soil = 295.0u"K"
    for air in (canopy_air, reference_air)
        air_conditions = ground_conditions(air)
        absorbed = litter_shortwave(wet_litter, air_conditions.direct_shortwave, air_conditions.diffuse_shortwave,
            air_conditions.sun_elevation_sine, 0.2).litter_absorbed
        dry = litter_temperature(wet_litter, wet_soil(0.0), air_conditions, absorbed, T_soil).temperature
        damp = litter_temperature(wet_litter, wet_soil(1.0), air_conditions, absorbed, T_soil).temperature
        @test damp < dry
        wind_speed = litter_wind_speed(wet_litter, air)
        flux(T) = litter_net_heat_flux(wet_litter, wet_soil(1.0), air_conditions, absorbed, wind_speed, T, T_soil)
        @test flux(damp - wet_litter.tolerance) > 0.0u"W/m^2" > flux(damp + wet_litter.tolerance)
        @test litter_water_exchange(wet_soil(1.0), air_conditions, T_soil).vapor_flux > 0.0u"kg/m^2/s"
        # Dry litter takes up water vapor from the air
        @test litter_water_exchange(wet_soil(0.0), air_conditions, T_soil).vapor_flux < 0.0u"kg/m^2/s"
    end
    # Dry litter exchanges no water
    @test litter_water_exchange(soil(), conditions, T_soil).vapor_flux == 0.0u"kg/m^2/s"
end

@testset "wet litter: water is conserved as it catches rain and evaporates" begin
    T = 293.0u"K"
    capacity = litter_maximum_water_content(wet, T) * wet_litter.dry_mass
    for water_content in (0.0, 0.5, 2.0), rain in (0.0u"kg/m^2", 0.5u"kg/m^2", 20.0u"kg/m^2"),
            vapor_flux in (-1.0e-4u"kg/m^2/s", 0.0u"kg/m^2/s", 1.0e-5u"kg/m^2/s")
        duration = 3600.0u"s"
        step = litter_water_step(wet_litter, water_content; vapor_flux, temperature = T, rain, duration, buried = false)
        before = water_content * wet_litter.dry_mass
        after = step.water_content * wet_litter.dry_mass
        evaporated = min(vapor_flux * duration, before)
        @test after + step.throughfall ≈ before + rain - evaporated atol = 1e-12u"kg/m^2"
        @test 0.0u"kg/m^2" <= after <= capacity * (1 + 1e-12)
        @test step.throughfall >= (1.0 - wet_litter.cover_fraction) * rain - 1e-12u"kg/m^2"
    end
    # Under snow the litter keeps its water and lets the rain pass
    buried = litter_water_step(wet_litter, 0.5; vapor_flux = 1.0e-5u"kg/m^2/s", temperature = T, rain = 1.0u"kg/m^2",
        duration = 3600.0u"s", buried = true)
    @test buried == (; water_content = 0.5, throughfall = 1.0u"kg/m^2")
    @test litter_water_step(litter, 0.0; vapor_flux = 0.0u"kg/m^2/s", temperature = T, rain = 1.0u"kg/m^2",
        duration = 3600.0u"s", buried = false) == (; water_content = 0.0, throughfall = 1.0u"kg/m^2")
end

@testset "wet litter: vapor balances in the litter's air spaces, and wet litter holds the soil's back" begin
    vapor = LitterVapor(; covered_fraction = 0.9, half_resistance = 200.0u"s/m", element_conductance = 1 / 500.0u"s/m",
        element_vapor_density = 0.02u"kg/m^3", soil_wetness = 0.4)
    h_air = 0.01u"m/s"
    soil_density, air_density = 0.015u"kg/m^3", 0.008u"kg/m^3"
    pore = pore_vapor_density(vapor, h_air, soil_density, air_density)
    covered, gap, wetness = 0.9, 0.1, 0.4
    from_soil = wetness * (covered * (soil_density - pore) / 200.0u"s/m" + gap * h_air * (soil_density - air_density))
    from_elements = litter_element_evaporation(vapor, h_air, soil_density, air_density)
    to_air = covered * (pore - air_density) / (200.0u"s/m" + 1 / h_air) + gap * wetness * h_air * (soil_density - air_density)
    @test from_soil + from_elements ≈ to_air
    # The ground's single source gives the same flux to the air
    source = ground_vapor_source(vapor, h_air, soil_density, air_density)
    @test h_air * (source - air_density) ≈ to_air
    # What the soil evaporates into is moister under wet litter than under dry
    T_air, P = 293.0u"K", 101325.0u"Pa"
    saturated = Microclimate.wet_air_properties(T_air, 1.0, P; vapour_pressure_equation = GoffGratch()).vapour_density
    wet_exchange = soil_vapor_exchange(vapor, h_air, T_air, air_density / saturated, P, GoffGratch())
    dry_exchange = soil_vapor_exchange(LitterVapor(; covered_fraction = 0.9, half_resistance = 200.0u"s/m",
        element_conductance = 0.0u"m/s", element_vapor_density = 0.0u"kg/m^3", soil_wetness = 0.4), h_air, T_air,
        air_density / saturated, P, GoffGratch())
    @test wet_exchange.relative_humidity > dry_exchange.relative_humidity
    # and the soil's wet fraction's flux into it is its flux into the litter's spaces, plus between the elements
    soil_flux = wetness * wet_exchange.mass_transfer_coefficient * (soil_density - wet_exchange.relative_humidity * saturated)
    @test soil_flux ≈ from_soil
    T_soil = 295.0u"K"
    @test allocations(soil_surface_heat_flux, wet_litter, wet_soil(1.0), conditions, BulkTransferEvaporation(), T_soil) == 0
    @test allocations(ground_surface, wet_soil(1.0), conditions, T_soil) == 0
end

@testset "a simulation with wet litter runs and stays cooler than with dry litter" begin
    bare_problem = example_microclimate_problem(; heights = [0.01, 2.0]u"m")
    damp_litter = LitterLayer(; dry_mass = 1.5u"kg/m^2", thickness = 2.5u"cm", cover_fraction = 0.9, albedo = 0.3,
        moisture = WetLitter(; element_vapor_resistance = 5000.0u"s/m", initial_water_content = 1.0))
    damp = solve(with_litter(bare_problem, damp_litter))
    dry = solve(with_litter(bare_problem, litter))
    @test all(isfinite, ustrip.(damp.soil_temperature))
    @test all(>=(0.0), damp.litter_water_content)
    @test all(iszero, dry.litter_water_content)
    @test sum(damp.soil_temperature[:, 1]) < sum(dry.soil_temperature[:, 1])
    # Litter that starts dry takes up water from the air
    drying = LitterLayer(; dry_mass = 1.5u"kg/m^2", thickness = 2.5u"cm", cover_fraction = 0.9, albedo = 0.3,
        moisture = WetLitter(; element_vapor_resistance = 5000.0u"s/m"))
    started_dry = solve(with_litter(bare_problem, drying))
    @test all(isfinite, ustrip.(started_dry.soil_temperature))
    @test maximum(started_dry.litter_water_content) > 0.0
end

@testset "litter slows the soil's drying, as a mulch" begin
    n = length(Microclimate.DEFAULT_DEPTHS)
    bare_problem = example_microclimate_problem(; heights = [0.01, 2.0]u"m", soil_moisture_strategy = DynamicSoilMoisture(),
        initial_soil_moisture = fill(0.3, n))
    bare = solve(bare_problem)
    covered = solve(with_litter(bare_problem, litter))
    surface_water(result) = sum(result.soil_moisture[end, 1:4])
    @test surface_water(covered) > surface_water(bare)
end
