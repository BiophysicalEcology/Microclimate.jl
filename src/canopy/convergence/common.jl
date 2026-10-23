# Placeholder hourly forcing that fixes the solver cache's types
initial_residual_forcing(model, soil_hydraulic_model) = CanopyResidualForcing(;
    atmospheric_pressure=101325.0u"Pa", reference_temperature=288.0u"K", reference_wind_speed=1.0u"m/s",
    reference_humidity=0.5, incoming_longwave=300.0u"W/m^2", ground_temperature=288.0u"K",
    ground_emissivity=0.95, ground_relative_humidity=0.5, leaf_water_potential=0.0u"J/kg",
    dry_conductance=stomatal_conductance(model.stomatal_model, 0.0u"°", 0.0u"J/kg", soil_hydraulic_model),
    canopy_top_air_temperature=288.0u"K", canopy_top_relative_humidity=0.5)

"""
    canopy_hour_setup!(state, forcing::Ref, buffers, model::MultilayerCanopy, inputs::CanopyEnergyBalanceInputs,
        initial_wind; site, vapour_pressure_equation, soil_hydraulic_model)

Write the hour's [`CanopyResidualForcing`](@ref) into `forcing` and warm-start `state` from
the previous hour's solution.
"""
function canopy_hour_setup!(state, forcing::Ref, buffers, model::MultilayerCanopy, inputs::CanopyEnergyBalanceInputs,
    initial_wind; site, vapour_pressure_equation, soil_hydraulic_model,
)
    (; environment_instant, zenith_angle, ground_temperature, ground_emissivity, ground_relative_humidity,
       leaf_water_potential) = inputs
    (; atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity) = environment_instant
    leaf_temperature = buffers.leaf.leaf_temperature
    air = buffers.air_profile
    n_layers = length(leaf_temperature)

    forcing[] = CanopyResidualForcing(; atmospheric_pressure, reference_temperature, reference_wind_speed,
        reference_humidity, ground_temperature, ground_emissivity, ground_relative_humidity, leaf_water_potential,
        incoming_longwave=canopy_incoming_longwave(model; site, environment_instant, vapour_pressure_equation),
        dry_conductance=stomatal_conductance(model.stomatal_model, zenith_angle, leaf_water_potential, soil_hydraulic_model),
        canopy_top_air_temperature=initial_wind.canopy_top_air_temperature,
        canopy_top_relative_humidity=initial_wind.canopy_top_relative_humidity)

    # Warm start from last hour's solution; seed vapor density on the first call.
    if all(iszero, air.vapor_density)
        fill!(air.vapor_density, wet_air_properties(initial_wind.canopy_top_air_temperature,
            initial_wind.canopy_top_relative_humidity, atmospheric_pressure; vapour_pressure_equation).vapour_density)
    end
    for i in 1:n_layers
        state[i] = ustrip(u"K", leaf_temperature[i])
        state[n_layers + i] = ustrip(u"K", air.air_temperature[i])
        state[2n_layers + i] = ustrip(u"g/m^3", air.vapor_density[i])
    end
    state[3n_layers + 1] = ustrip(u"m^-1", initial_wind.inverse_obukhov_length)
    return state
end

"""
    canopy_solution!(residual, state, buffers, model::MultilayerCanopy, boundary_layer_model,
        forcing::CanopyResidualForcing, inputs::CanopyEnergyBalanceInputs, iterations;
        site, vapour_pressure_equation, soil_hydraulic_model)

Leave `buffers` at the solved `state` and return [`converge_canopy!`](@ref)'s result.
"""
function canopy_solution!(residual, state, buffers, model::MultilayerCanopy, boundary_layer_model,
    forcing::CanopyResidualForcing, inputs::CanopyEnergyBalanceInputs, iterations;
    site, vapour_pressure_equation, soil_hydraulic_model,
)
    (; zenith_angle) = inputs
    (; atmospheric_pressure) = inputs.environment_instant
    (; leaf_temperature, potential_evaporation_mass_flow, leaf_body) = buffers.leaf
    air = buffers.air_profile
    (; layer_transmission) = buffers.longwave
    (; layer_plant_area_index_boundaries) = buffers.shortwave
    n_layers = length(leaf_temperature)

    # Leave buffers at the solution: the residual's own intermediates, with
    # the air profile at the solved state rather than its transport.
    result = evaluate_canopy!(residual, state, buffers, model, boundary_layer_model, forcing;
        site, vapour_pressure_equation)
    unstressed_dry_conductance = stomatal_conductance(model.stomatal_model, zenith_angle, 0.0u"J/kg", soil_hydraulic_model)
    for i in 1:n_layers
        air.air_temperature[i] = air.air_temperature_prev[i]
        air.vapor_density[i] = air.vapor_density_prev[i]
        saturation = wet_air_properties(air.air_temperature[i], 1.0, atmospheric_pressure; vapour_pressure_equation).vapour_density
        air.relative_humidity[i] = clamp(air.vapor_density[i] / saturation, 0.0, 1.0)

        # Unstressed (water_potential=0) reference transpiration, for a
        # soil-hydraulics demand term — not the leaf's own energy balance.
        exchange_fraction = 1.0 - layer_transmission[i]
        if iszero(exchange_fraction)
            potential_evaporation_mass_flow[i] = zero(potential_evaporation_mass_flow[i])
            continue
        end
        layer_pai = layer_plant_area_index_boundaries[i + 1] - layer_plant_area_index_boundaries[i]
        aerodynamic_area = 2.0 * layer_pai * 1.0u"m^2"
        conv = leaf_convection(model.leaf_convection_model, leaf_body, aerodynamic_area, air.air_temperature[i],
            leaf_temperature[i], buffers.wind.wind_speed[i], atmospheric_pressure)
        atmos = AtmosphericConditions(air.relative_humidity[i], buffers.wind.wind_speed[i], atmospheric_pressure)
        potential_evaporation_mass_flow[i] = HeatExchange.evaporation(unstressed_dry_conductance,
            conv.mass_transfer_coefficient, atmos, aerodynamic_area, leaf_temperature[i], air.air_temperature[i];
            water_potential=0.0u"J/kg").transpiration_mass_flow
    end

    return result.longwave_result, iterations, result.canopy_top_air_temperature,
        result.canopy_top_relative_humidity, result.friction_velocity, result.inverse_obukhov_length
end
