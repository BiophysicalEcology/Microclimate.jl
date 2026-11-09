# The tolerance for each kind of residual, as plain numbers in the state's units
canopy_tolerances(convergence_model) = (; leaf_balance=in_state_units(convergence_model.leaf_balance_tolerance),
    air_temperature=in_state_units(convergence_model.air_temperature_tolerance),
    vapor_density=in_state_units(convergence_model.vapor_density_tolerance),
    ground_temperature=in_state_units(convergence_model.ground_temperature_tolerance),
    inverse_obukhov_length=in_state_units(convergence_model.inverse_obukhov_length_tolerance))

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

Warn, keep the best guess the solve reached, and carry on.
"""
struct SmallestResidual <: AbstractCanopyFailure end

# The largest residual from `evaluate_canopy!`, measured in units of its own tolerance.
# The ground temperature residual needs the soil, so the solver adds it.
function canopy_scaled_residual(residual, tolerances)
    largest = abs(in_state_units(residual.inverse_obukhov_length[])) / tolerances.inverse_obukhov_length
    for layer in eachindex(residual.leaf_balance)
        largest = max(largest,
            abs(in_state_units(residual.leaf_balance[layer])) / tolerances.leaf_balance,
            abs(in_state_units(residual.leafless_temperature_difference[layer])) / tolerances.air_temperature,
            abs(in_state_units(residual.air_temperature[layer])) / tolerances.air_temperature,
            abs(in_state_units(residual.vapor_density[layer])) / tolerances.vapor_density)
    end
    return largest
end

# Handle an hour that didn't finish, as the convergence model's `on_failure` asks
canopy_failure!(::ThrowError, convergence_model, return_code, best_state, smallest_residual, iterations, buffers, model,
    boundary_layer_model, inputs, hour_conditions; kw...) = canopy_convergence_error(convergence_model, return_code, inputs)
function canopy_failure!(::SmallestResidual, convergence_model, return_code, best_state, smallest_residual, iterations, buffers,
    model, boundary_layer_model, inputs, hour_conditions; site, vapour_pressure_equation, soil_hydraulic_model,
)
    (; state, residual) = buffers.convergence
    copyto!(state, best_state)
    result = canopy_solution!(residual, state, buffers, model, boundary_layer_model, inputs, hour_conditions, iterations;
        site, vapour_pressure_equation, soil_hydraulic_model)
    @warn "Canopy solve did not converge, keeping the state with the smallest residual" return_code iterations residual_over_tolerance =
        smallest_residual maxlog = 5
    return result
end

canopy_convergence_error(convergence_model, return_code, inputs) =
    error(nameof(typeof(convergence_model)), " did not converge: ", return_code,
        " (reference temperature ", inputs.environment_instant.reference_temperature,
        ", wind ", inputs.environment_instant.reference_wind_speed,
        ", ground temperature ", inputs.ground_temperature, ", zenith angle ", inputs.zenith_angle, ")")

"""
    canopy_hour_conditions(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs,
        initial_wind; site, vapour_pressure_equation, soil_hydraulic_model)

What stays fixed while the canopy is solved for one hour, as `(; incoming_longwave, dry_conductance,
stability_limits, canopy_top_air_temperature, canopy_top_relative_humidity)`. The canopy-top values are
those when the hour starts. `stability_limits` are the most unstable and most stable inverse Obukhov
lengths allowed, `(; unstable, stable)`.
"""
canopy_hour_conditions(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs,
    initial_wind; site, vapour_pressure_equation, soil_hydraulic_model) = (;
    incoming_longwave=canopy_incoming_longwave(model; site, inputs.environment_instant, vapour_pressure_equation),
    dry_conductance=stomatal_conductance(model.stomatal_model, inputs.zenith_angle, inputs.leaf_water_potential,
        soil_hydraulic_model),
    stability_limits=canopy_stability_limits(buffers, boundary_layer_model),
    canopy_top_air_temperature=initial_wind.canopy_top_air_temperature,
    canopy_top_relative_humidity=initial_wind.canopy_top_relative_humidity)

"""
    canopy_starting_state!(state, buffers, inputs::CanopyEnergyBalanceInputs, initial_wind, stability_limits;
        vapour_pressure_equation)

Starts `state` from the previous hour's answer.
"""
function canopy_starting_state!(state, buffers, inputs::CanopyEnergyBalanceInputs, initial_wind, stability_limits;
    vapour_pressure_equation,
)
    air = buffers.air_profile
    # The first hour has no previous answer for vapor density
    if all(iszero, air.vapor_density)
        fill!(air.vapor_density, wet_air_properties(initial_wind.canopy_top_air_temperature,
            initial_wind.canopy_top_relative_humidity, inputs.environment_instant.atmospheric_pressure;
            vapour_pressure_equation).vapour_density)
    end
    return pack_canopy_state!(state, buffers.leaf.leaf_temperature, air.air_temperature, air.vapor_density,
        inputs.ground_temperature, limit_inverse_obukhov_length(initial_wind.inverse_obukhov_length, stability_limits))
end

"""
    canopy_solution!(residual, state, buffers, model::MultilayerCanopy, boundary_layer_model,
        inputs::CanopyEnergyBalanceInputs, hour_conditions, iterations;
        site, vapour_pressure_equation, soil_hydraulic_model)

Fills `buffers` with the solved `state` and returns the result of [`converge_canopy!`](@ref).
"""
function canopy_solution!(residual, state, buffers, model::MultilayerCanopy, boundary_layer_model,
    inputs::CanopyEnergyBalanceInputs, hour_conditions, iterations;
    site, vapour_pressure_equation, soil_hydraulic_model,
)
    (; zenith_angle) = inputs
    (; atmospheric_pressure) = inputs.environment_instant
    (; leaf_temperature, potential_evaporation_mass_flux, leaf_body) = buffers.leaf
    air = buffers.air_profile
    n_layers = length(leaf_temperature)

    # Evaluating fills the buffers, but leaves the air profile at what the leaves
    # produce, so it is then set back to the solved values
    (; ground_temperature, inverse_obukhov_length) = unpack_canopy_state!(buffers, state)
    result = evaluate_canopy!(residual, inverse_obukhov_length, ground_temperature, buffers, model, boundary_layer_model,
        inputs, hour_conditions; site, vapour_pressure_equation)
    unstressed_dry_conductance = stomatal_conductance(model.stomatal_model, zenith_angle, 0.0u"J/kg", soil_hydraulic_model)
    for i in 1:n_layers
        air.air_temperature[i] = air.current_air_temperature[i]
        air.vapor_density[i] = air.current_vapor_density[i]
        air.relative_humidity[i] = layer_relative_humidity(air, i, atmospheric_pressure, vapour_pressure_equation)

        # What the leaves would transpire with no water stress, which the soil hydraulics
        # model uses as demand
        (; exchange_fraction, aerodynamic_area) = layer_leaf_areas(buffers, i)
        if iszero(exchange_fraction)
            potential_evaporation_mass_flux[i] = zero(potential_evaporation_mass_flux[i])
            continue
        end
        convection = leaf_convection(model.leaf_convection_model, leaf_body, aerodynamic_area, air.air_temperature[i],
            leaf_temperature[i], buffers.wind.wind_speed[i], atmospheric_pressure)
        atmosphere = AtmosphericConditions(; vapour_density=air.vapor_density[i], wind_speed=buffers.wind.wind_speed[i],
            atmospheric_pressure)
        potential_evaporation_mass_flux[i] = HeatExchange.evaporation(unstressed_dry_conductance,
            convection.mass_transfer_coefficient, atmosphere, aerodynamic_area, leaf_temperature[i], air.air_temperature[i];
            water_potential=0.0u"J/kg").transpiration_mass_flow / GROUND_AREA
    end

    return (; result..., iterations)
end
