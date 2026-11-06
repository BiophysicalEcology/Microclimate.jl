# The tolerance for each kind of residual row, in the state's units (see canopy_state_length)
canopy_tolerances(convergence_model) = (; leaf_balance=ustrip(u"W/m^2", convergence_model.leaf_balance_tolerance),
    air_temperature=ustrip(u"K", convergence_model.air_temperature_tolerance),
    vapor_density=ustrip(u"g/m^3", convergence_model.vapor_density_tolerance),
    inverse_obukhov_length=ustrip(u"m^-1", convergence_model.inverse_obukhov_length_tolerance))

# The tolerance of residual row `row`, in the state's units
function canopy_row_tolerance(tolerances, buffers, row)
    rows = canopy_state_rows(length(buffers.leaf.leaf_temperature))
    if row in rows.leaf_temperature
        # A layer without leaves gives a temperature difference. Leaf rows come first, so the row is the layer.
        iszero(layer_leaf_areas(buffers, row).exchange_fraction) && return tolerances.air_temperature
        return tolerances.leaf_balance
    elseif row in rows.air_temperature
        return tolerances.air_temperature
    elseif row in rows.vapor_density
        return tolerances.vapor_density
    else
        return tolerances.inverse_obukhov_length
    end
end

"""
    AbstractCanopyFailure

What to do with an hour the canopy solve can't finish: [`ThrowError`](@ref) or
[`SmallestResidual`](@ref).
"""
abstract type AbstractCanopyFailure end

"""
    ThrowError()

Stop with an error that lists the hour's conditions.
"""
struct ThrowError <: AbstractCanopyFailure end

"""
    SmallestResidual()

Keep the best guess the solve reached and carry on.
"""
struct SmallestResidual <: AbstractCanopyFailure end

# The largest residual, measured in units of its own tolerance. At or below 1, the hour is solved.
function canopy_scaled_residual(residual, tolerances, buffers)
    largest = 0.0
    for row in eachindex(residual)
        largest = max(largest, abs(residual[row]) / canopy_row_tolerance(tolerances, buffers, row))
    end
    return largest
end

# Handle an hour that didn't finish, as the convergence model's `on_failure` asks
canopy_failure!(::ThrowError, convergence_model, return_code, best_state, buffers, model, boundary_layer_model,
    forcing, inputs; kw...) = canopy_convergence_error(convergence_model, return_code, inputs)
function canopy_failure!(::SmallestResidual, convergence_model, return_code, best_state, buffers, model, boundary_layer_model,
    forcing, inputs; site, vapour_pressure_equation, soil_hydraulic_model,
)
    (; state, residual) = buffers.convergence
    copyto!(state, best_state)
    return canopy_solution!(residual, state, buffers, model, boundary_layer_model, forcing, inputs, convergence_model.maxiters;
        site, vapour_pressure_equation, soil_hydraulic_model)
end

canopy_convergence_error(convergence_model, return_code, inputs) =
    error(nameof(typeof(convergence_model)), " did not converge: ", return_code,
        " (reference temperature ", inputs.environment_instant.reference_temperature,
        ", wind ", inputs.environment_instant.reference_wind_speed,
        ", ground temperature ", inputs.ground_temperature, ", zenith angle ", inputs.zenith_angle, ")")

# A placeholder hour, so the solver's buffers can be built with the right types
initial_residual_forcing(model, soil_hydraulic_model) = CanopyResidualForcing(;
    atmospheric_pressure=101325.0u"Pa", reference_temperature=288.0u"K", reference_wind_speed=1.0u"m/s",
    reference_humidity=0.5, incoming_longwave=300.0u"W/m^2", ground_temperature=288.0u"K",
    ground_emissivity=0.95, ground_relative_humidity=0.5, leaf_water_potential=0.0u"J/kg",
    dry_conductance=stomatal_conductance(model.stomatal_model, 0.0u"°", 0.0u"J/kg", soil_hydraulic_model),
    canopy_top_air_temperature=288.0u"K", canopy_top_relative_humidity=0.5,
    stability_limits=(; unstable=-Inf * u"m^-1", stable=Inf * u"m^-1"))

"""
    canopy_hour_setup!(state, forcing::Ref, buffers, model::MultilayerCanopy, boundary_layer_model,
        inputs::CanopyEnergyBalanceInputs, initial_wind; site, vapour_pressure_equation, soil_hydraulic_model)

Prepares the solve for one hour: writes the hour's [`CanopyResidualForcing`](@ref) into
`forcing`, and starts `state` from the previous hour's answer.
"""
function canopy_hour_setup!(state, forcing::Ref, buffers, model::MultilayerCanopy, boundary_layer_model,
    inputs::CanopyEnergyBalanceInputs, initial_wind; site, vapour_pressure_equation, soil_hydraulic_model,
)
    (; environment_instant, zenith_angle, ground_temperature, ground_emissivity, ground_relative_humidity,
       leaf_water_potential) = inputs
    (; atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity) = environment_instant
    leaf_temperature = buffers.leaf.leaf_temperature
    air = buffers.air_profile
    rows = canopy_state_rows(length(leaf_temperature))

    forcing[] = CanopyResidualForcing(; atmospheric_pressure, reference_temperature, reference_wind_speed,
        reference_humidity, ground_temperature, ground_emissivity, ground_relative_humidity, leaf_water_potential,
        incoming_longwave=canopy_incoming_longwave(model; site, environment_instant, vapour_pressure_equation),
        dry_conductance=stomatal_conductance(model.stomatal_model, zenith_angle, leaf_water_potential, soil_hydraulic_model),
        canopy_top_air_temperature=initial_wind.canopy_top_air_temperature,
        canopy_top_relative_humidity=initial_wind.canopy_top_relative_humidity,
        stability_limits=canopy_stability_limits(buffers, boundary_layer_model))

    # The first hour has no previous answer for vapor density
    if all(iszero, air.vapor_density)
        fill!(air.vapor_density, wet_air_properties(initial_wind.canopy_top_air_temperature,
            initial_wind.canopy_top_relative_humidity, atmospheric_pressure; vapour_pressure_equation).vapour_density)
    end
    for i in eachindex(leaf_temperature)
        state[rows.leaf_temperature[i]] = ustrip(u"K", leaf_temperature[i])
        state[rows.air_temperature[i]] = ustrip(u"K", air.air_temperature[i])
        state[rows.vapor_density[i]] = ustrip(u"g/m^3", air.vapor_density[i])
    end
    state[rows.inverse_obukhov_length] = ustrip(u"m^-1",
        limit_inverse_obukhov_length(initial_wind.inverse_obukhov_length, forcing[].stability_limits))
    return state
end

"""
    canopy_solution!(residual, state, buffers, model::MultilayerCanopy, boundary_layer_model,
        forcing::CanopyResidualForcing, inputs::CanopyEnergyBalanceInputs, iterations;
        site, vapour_pressure_equation, soil_hydraulic_model)

Fills `buffers` with the solved `state` and returns the result of [`converge_canopy!`](@ref).
"""
function canopy_solution!(residual, state, buffers, model::MultilayerCanopy, boundary_layer_model,
    forcing::CanopyResidualForcing, inputs::CanopyEnergyBalanceInputs, iterations;
    site, vapour_pressure_equation, soil_hydraulic_model,
)
    (; zenith_angle) = inputs
    (; atmospheric_pressure) = inputs.environment_instant
    (; leaf_temperature, potential_evaporation_mass_flow, leaf_body) = buffers.leaf
    air = buffers.air_profile
    n_layers = length(leaf_temperature)

    # Evaluating fills the buffers, but leaves the air profile at what the leaves
    # produce, so it is then set back to the solved values
    result = evaluate_canopy!(residual, state, buffers, model, boundary_layer_model, forcing;
        site, vapour_pressure_equation)
    unstressed_dry_conductance = stomatal_conductance(model.stomatal_model, zenith_angle, 0.0u"J/kg", soil_hydraulic_model)
    for i in 1:n_layers
        air.air_temperature[i] = air.current_air_temperature[i]
        air.vapor_density[i] = air.current_vapor_density[i]
        saturation_vapor_density = wet_air_properties(air.air_temperature[i], 1.0, atmospheric_pressure;
            vapour_pressure_equation).vapour_density
        air.relative_humidity[i] = clamp(air.vapor_density[i] / saturation_vapor_density, 0.0, 1.0)

        # What the leaves would transpire with no water stress, which the soil hydraulics
        # model uses as demand
        (; exchange_fraction, aerodynamic_area) = layer_leaf_areas(buffers, i)
        if iszero(exchange_fraction)
            potential_evaporation_mass_flow[i] = zero(potential_evaporation_mass_flow[i])
            continue
        end
        convection = leaf_convection(model.leaf_convection_model, leaf_body, aerodynamic_area, air.air_temperature[i],
            leaf_temperature[i], buffers.wind.wind_speed[i], atmospheric_pressure)
        atmosphere = AtmosphericConditions(; vapour_density=air.vapor_density[i], wind_speed=buffers.wind.wind_speed[i],
            atmospheric_pressure)
        potential_evaporation_mass_flow[i] = HeatExchange.evaporation(unstressed_dry_conductance,
            convection.mass_transfer_coefficient, atmosphere, aerodynamic_area, leaf_temperature[i], air.air_temperature[i];
            water_potential=0.0u"J/kg").transpiration_mass_flow
    end

    return result.longwave_result, iterations, result.canopy_top_air_temperature,
        result.canopy_top_relative_humidity, result.friction_velocity, result.inverse_obukhov_length
end
