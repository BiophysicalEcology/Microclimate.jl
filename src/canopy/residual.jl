"""
    CanopyResidualForcing(; atmospheric_pressure, reference_temperature, reference_wind_speed,
        reference_humidity, incoming_longwave, ground_emissivity,
        ground_relative_humidity, leaf_water_potential, dry_conductance,
        canopy_top_air_temperature, canopy_top_relative_humidity, stability_limits)

Everything [`evaluate_canopy!`](@ref) needs for the hour that doesn't change as the canopy
is solved. `canopy_top_air_temperature` and `canopy_top_relative_humidity` are the values at
the canopy top when the hour starts. `stability_limits` are the most unstable and most stable
inverse Obukhov lengths allowed, `(; unstable, stable)`.
"""
@kwdef struct CanopyResidualForcing{C}
    atmospheric_pressure::typeof(1.0u"Pa")
    reference_temperature::typeof(1.0u"K")
    reference_wind_speed::typeof(1.0u"m/s")
    reference_humidity::Float64
    incoming_longwave::typeof(1.0u"W/m^2")
    ground_emissivity::Float64
    ground_relative_humidity::Float64
    leaf_water_potential::typeof(1.0u"J/kg")
    dry_conductance::C
    canopy_top_air_temperature::typeof(1.0u"K")
    canopy_top_relative_humidity::Float64
    stability_limits::@NamedTuple{unstable::typeof(1.0u"m^-1"), stable::typeof(1.0u"m^-1")}
end
# Convert each field to its fixed unit, inferring only the conductance type
CanopyResidualForcing(atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity,
    incoming_longwave, ground_emissivity, ground_relative_humidity, leaf_water_potential,
    dry_conductance::C, canopy_top_air_temperature, canopy_top_relative_humidity, stability_limits) where C =
    CanopyResidualForcing{C}(atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity,
        incoming_longwave, ground_emissivity, ground_relative_humidity, leaf_water_potential,
        dry_conductance, canopy_top_air_temperature, canopy_top_relative_humidity, stability_limits)

Unitful.@derived_dimension AreaPower Unitful.𝐌 * Unitful.𝐓^-3

# The solver works on one vector of plain numbers, in these units
const CANOPY_STATE_UNITS = (; temperature=u"K", vapor_density=u"g/m^3", inverse_obukhov_length=u"m^-1", leaf_balance=u"W/m^2")

in_state_units(value::Unitful.Temperature) = ustrip(CANOPY_STATE_UNITS.temperature, value)
in_state_units(value::Unitful.Density) = ustrip(CANOPY_STATE_UNITS.vapor_density, value)
in_state_units(value::Unitful.Wavenumber) = ustrip(CANOPY_STATE_UNITS.inverse_obukhov_length, value)
in_state_units(value::AreaPower) = ustrip(CANOPY_STATE_UNITS.leaf_balance, value)

# The solver's vector holds a leaf temperature, air temperature and vapor density for each
# layer, then one ground temperature and one inverse Obukhov length. Only the three functions
# below know that layout.
canopy_state_length(layer_count) = 3layer_count + 2

"""
    pack_canopy_state!(state, leaf_temperature, air_temperature, vapor_density, ground_temperature, inverse_obukhov_length)

Writes one value of each kind per layer, the ground temperature and the inverse Obukhov length
into the solver's vector.
"""
function pack_canopy_state!(state, leaf_temperature, air_temperature, vapor_density, ground_temperature, inverse_obukhov_length)
    layer_count = length(leaf_temperature)
    for layer in 1:layer_count
        state[layer] = in_state_units(leaf_temperature[layer])
        state[layer_count + layer] = in_state_units(air_temperature[layer])
        state[2layer_count + layer] = in_state_units(vapor_density[layer])
    end
    state[end - 1] = in_state_units(ground_temperature)
    state[end] = in_state_units(inverse_obukhov_length)
    return state
end

"""
    unpack_canopy_state!(buffers, state)

Writes the solver's vector into the leaf temperature, air temperature and vapor density in
`buffers`, and returns `(; ground_temperature, inverse_obukhov_length)`.
"""
function unpack_canopy_state!(buffers, state)
    (; leaf_temperature) = buffers.leaf
    (; current_air_temperature, current_vapor_density) = buffers.air_profile
    layer_count = length(leaf_temperature)
    for layer in 1:layer_count
        leaf_temperature[layer] = state[layer] * CANOPY_STATE_UNITS.temperature
        current_air_temperature[layer] = state[layer_count + layer] * CANOPY_STATE_UNITS.temperature
        current_vapor_density[layer] = state[2layer_count + layer] * CANOPY_STATE_UNITS.vapor_density
    end
    return (; ground_temperature=state[end - 1] * CANOPY_STATE_UNITS.temperature,
        inverse_obukhov_length=state[end] * CANOPY_STATE_UNITS.inverse_obukhov_length)
end

# A solver's vector with the same value for every layer, for the limits and scales that depend
# only on the kind of value
canopy_state_filled(layer_count, leaf_temperature, air_temperature, vapor_density, ground_temperature, inverse_obukhov_length) =
    pack_canopy_state!(zeros(canopy_state_length(layer_count)), fill(leaf_temperature, layer_count),
        fill(air_temperature, layer_count), fill(vapor_density, layer_count), ground_temperature, inverse_obukhov_length)

"""
    allocate_canopy_residual(layer_count)

Vectors for how far each layer is from solved, written by [`evaluate_canopy!`](@ref):

- `leaf_balance`: each leaf's net energy balance per unit leaf area. Zero in a layer without leaves.
- `leafless_temperature_difference`: leaf minus air temperature in a layer without leaves, which
  has no energy balance of its own. Zero everywhere else.
- `air_temperature`: air temperature minus what the air-profile model gets from the leaves.
- `vapor_density`: vapor density minus what the air-profile model gets from the leaves.
- `ground_temperature`: one value, minus what the soil reaches by the end of the hour.
- `inverse_obukhov_length`: one value for the canopy, minus what the canopy's heat flux implies.
"""
allocate_canopy_residual(layer_count) = (;
    leaf_balance=zeros(typeof(0.0u"W/m^2"), layer_count),
    leafless_temperature_difference=zeros(typeof(0.0u"K"), layer_count),
    air_temperature=zeros(typeof(0.0u"K"), layer_count),
    vapor_density=zeros(typeof(0.0u"g/m^3"), layer_count),
    ground_temperature=Ref(0.0u"K"),
    inverse_obukhov_length=Ref(0.0u"m^-1"))

# Marks an evaluation as unusable, so that the solver treats it as an overshoot
function invalidate_canopy_residual!(residual)
    fill!(residual.leaf_balance, NaN * u"W/m^2")
    fill!(residual.leafless_temperature_difference, NaN * u"K")
    fill!(residual.air_temperature, NaN * u"K")
    fill!(residual.vapor_density, NaN * u"g/m^3")
    residual.ground_temperature[] = NaN * u"K"
    residual.inverse_obukhov_length[] = NaN * u"m^-1"
    return residual
end

"""
    layer_leaf_areas(buffers, layer)

The areas of a layer's leaves, as `(; exchange_fraction, leaf_area, aerodynamic_area)`.

`exchange_fraction` is the share of longwave radiation the layer stops. `leaf_area` is the
two-sided area that exchanges radiation, from that share. `aerodynamic_area` is the two-sided
area that exchanges heat and vapor with the air, from the layer's plant area index.
"""
function layer_leaf_areas(buffers, layer)
    exchange_fraction = 1.0 - buffers.longwave.layer_transmission[layer]
    boundaries = buffers.shortwave.layer_plant_area_index_boundaries
    layer_plant_area_index = boundaries[layer + 1] - boundaries[layer]
    return (; exchange_fraction, leaf_area=2.0 * exchange_fraction * 1.0u"m^2",
        aerodynamic_area=2.0 * layer_plant_area_index * 1.0u"m^2")
end

"""
    CANOPY_TEMPERATURE_BOUNDS

The lowest and highest temperatures a canopy state can have. Outside these,
[`evaluate_canopy!`](@ref) gives a NaN residual.
"""
const CANOPY_TEMPERATURE_BOUNDS = (100.0u"K", 1000.0u"K")

"""
    layer_leaf_balance(buffers, model::MultilayerCanopy, forcing::CanopyResidualForcing, layer, leaf_temperature)

A layer's leaf energy balance at `leaf_temperature`, as `(; balance, leaf_area)`.
"""
function layer_leaf_balance(buffers, model::MultilayerCanopy, forcing::CanopyResidualForcing, layer, leaf_temperature)
    (; atmospheric_pressure, leaf_water_potential, dry_conductance) = forcing
    air = buffers.air_profile
    air_temperature = air.current_air_temperature[layer]
    (; leaf_area, aerodynamic_area) = layer_leaf_areas(buffers, layer)
    conductance = blend_stomatal_conductance(dry_conductance, wet_canopy_fraction(model, buffers, layer))
    balance = leaf_heat_balance(leaf_temperature, buffers.leaf.absorbed_radiation[layer], air_temperature, air.current_vapor_density[layer],
        buffers.wind.wind_speed[layer], atmospheric_pressure, model.leaf_parameters.leaf_emissivity, conductance,
        leaf_water_potential, buffers.leaf.leaf_body, leaf_area, model.leaf_convection_model; aerodynamic_area)
    return (; balance, leaf_area)
end

# The relative humidity of the air the solver has in `layer`
function layer_relative_humidity(air_profile, layer, atmospheric_pressure, vapour_pressure_equation)
    saturation_vapor_density = wet_air_properties(air_profile.current_air_temperature[layer], 1.0, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    return clamp(air_profile.current_vapor_density[layer] / saturation_vapor_density, 0.0, 1.0)
end

# A value inside its bounds. NaN, from a failed step, goes to the lower bound.
_within_bounds(value, lower_bound, upper_bound=typemax(value)) =
    isnan(value) ? lower_bound : clamp(value, lower_bound, upper_bound)

"""
    evaluate_canopy!(residual, inverse_obukhov_length, ground_temperature, buffers, model,
        boundary_layer_model, forcing::CanopyResidualForcing; site, vapour_pressure_equation)

Measures how far the canopy in `buffers`, with this `inverse_obukhov_length` and
`ground_temperature`, is from solved. The hour is solved when every value written to `residual`
(see [`allocate_canopy_residual`](@ref)) is close to zero. `residual.ground_temperature` needs the
soil, so the solver writes that one.

The wind, radiation, leaf fluxes and air profile are left in `buffers`.

Returns `(; longwave_result, canopy_top_air_temperature,
canopy_top_relative_humidity, friction_velocity, inverse_obukhov_length)`.
"""
function evaluate_canopy!(residual, inverse_obukhov_length, ground_temperature, buffers, model::MultilayerCanopy, boundary_layer_model,
    forcing::CanopyResidualForcing; site, vapour_pressure_equation,
)
    (; atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity, incoming_longwave,
       ground_emissivity, ground_relative_humidity) = forcing
    (; leaf_temperature, sensible_heat_source, evaporation_mass_flow, absorbed_radiation, net_balance) = buffers.leaf
    air = buffers.air_profile
    (; absorbed_shortwave) = buffers.shortwave
    (; absorbed_longwave) = buffers.longwave
    n_layers = length(leaf_temperature)

    # A canopy out of bounds is evaluated at the bounds, and its residual is NaN
    lowest_temperature, highest_temperature = CANOPY_TEMPERATURE_BOUNDS
    in_bounds = !isnan(inverse_obukhov_length) && lowest_temperature <= ground_temperature <= highest_temperature
    ground_temperature = _within_bounds(ground_temperature, lowest_temperature, highest_temperature)
    for i in 1:n_layers
        in_bounds &= lowest_temperature <= leaf_temperature[i] <= highest_temperature &&
            lowest_temperature <= air.current_air_temperature[i] <= highest_temperature &&
            !isnan(air.current_vapor_density[i])
        leaf_temperature[i] = _within_bounds(leaf_temperature[i], lowest_temperature, highest_temperature)
        air.current_air_temperature[i] = _within_bounds(air.current_air_temperature[i], lowest_temperature, highest_temperature)
        air.current_vapor_density[i] = _within_bounds(air.current_vapor_density[i], zero(air.current_vapor_density[i]))
    end

    environment_instant = (; atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity)
    wind = canopy_wind_profile!(buffers, model, boundary_layer_model;
        site, environment_instant, canopy_source_temperature=leaf_temperature[1], vapour_pressure_equation,
        inverse_obukhov_length=isnan(inverse_obukhov_length) ? zero(inverse_obukhov_length) : inverse_obukhov_length)
    wind_speed = buffers.wind.wind_speed
    longwave_result = canopy_longwave!(buffers, model;
        leaf_temperature, ground_temperature, ground_emissivity, incoming_longwave)

    for i in 1:n_layers
        air_temperature = air.current_air_temperature[i]
        (; exchange_fraction) = layer_leaf_areas(buffers, i)
        if iszero(exchange_fraction)
            sensible_heat_source[i] = zero(sensible_heat_source[i])
            evaporation_mass_flow[i] = zero(evaporation_mass_flow[i])
            absorbed_radiation[i] = zero(absorbed_radiation[i])
            net_balance[i] = zero(net_balance[i])
            residual.leaf_balance[i] = zero(residual.leaf_balance[i])
            residual.leafless_temperature_difference[i] = leaf_temperature[i] - air_temperature
            continue
        end
        # absorbed_shortwave/absorbed_longwave are per unit ground area
        absorbed_radiation[i] = (absorbed_shortwave[i] + absorbed_longwave[i]) / (2.0 * exchange_fraction)
        (; balance, leaf_area) = layer_leaf_balance(buffers, model, forcing, i, leaf_temperature[i])
        sensible_heat_source[i] = balance.convection.convection_flow / 1.0u"m^2"
        evaporation_mass_flow[i] = balance.evaporation.transpiration_mass_flow
        net_balance[i] = balance.net
        residual.leaf_balance[i] = balance.net / leaf_area
        residual.leafless_temperature_difference[i] = zero(residual.leafless_temperature_difference[i])
    end

    reference_vapor_density = wet_air_properties(reference_temperature, reference_humidity, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    ground_vapor_density = wet_air_properties(ground_temperature, ground_relative_humidity, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    initial_top_vapor_density = wet_air_properties(forcing.canopy_top_air_temperature, forcing.canopy_top_relative_humidity,
        atmospheric_pressure; vapour_pressure_equation).vapour_density
    boundary = canopy_top_flux_boundary(boundary_layer_model, model.canopy_height, buffers.wind.displacement_height,
        buffers.wind.above_canopy_profile.heights[2], wind.friction_velocity, wind.inverse_obukhov_length,
        sum(sensible_heat_source), sum(evaporation_mass_flow) / 1.0u"m^2",
        ground_temperature, ground_vapor_density, reference_temperature, reference_vapor_density,
        forcing.canopy_top_air_temperature, initial_top_vapor_density)
    in_bounds &= lowest_temperature <= boundary.top_temperature <= highest_temperature &&
        boundary.top_vapor_density >= zero(boundary.top_vapor_density)
    canopy_top_air_temperature = clamp(boundary.top_temperature, lowest_temperature, highest_temperature)
    # The air profile takes vapor density, not relative humidity: humidity capped at 100%
    # would lose any extra vapor from the water balance
    canopy_top_vapor_density = max(boundary.top_vapor_density, zero(boundary.top_vapor_density))
    top_saturation = wet_air_properties(canopy_top_air_temperature, 1.0, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    canopy_top_relative_humidity = clamp(canopy_top_vapor_density / top_saturation, 0.0, 1.0)

    air_forcing = CanopyAirProfileForcing(; friction_velocity=wind.friction_velocity,
        inverse_obukhov_length=wind.inverse_obukhov_length, wind_speed, canopy_top_air_temperature,
        canopy_top_vapor_density, ground_temperature, ground_relative_humidity, sensible_heat_source,
        evaporation_mass_flow, atmospheric_pressure)
    canopy_heat_vapor_transport!(air, model.air_profile_model, boundary_layer_model, air_forcing;
        site, vapour_pressure_equation)

    for i in 1:n_layers
        residual.air_temperature[i] = air.current_air_temperature[i] - air.air_temperature[i]
        residual.vapor_density[i] = air.current_vapor_density[i] - air.vapor_density[i]
    end
    # The stability the whole canopy's sensible flux through its top implies; the Obukhov
    # relation takes the downward flux
    obukhov_coefficient = calc_ρ_cp(reference_temperature) * reference_temperature / (boundary_layer_model.karman_constant * g_n)
    inverse_obukhov_length_update = limit_inverse_obukhov_length(stability_from_fluxes(InverseObukhovLength(), obukhov_coefficient,
        wind.friction_velocity, -boundary.sensible_flux), forcing.stability_limits)
    residual.inverse_obukhov_length[] = wind.inverse_obukhov_length - inverse_obukhov_length_update
    in_bounds || invalidate_canopy_residual!(residual)

    return (; longwave_result, canopy_top_air_temperature, canopy_top_relative_humidity,
        friction_velocity=wind.friction_velocity, inverse_obukhov_length=wind.inverse_obukhov_length)
end

# The canopy-top inverse Obukhov lengths the boundary layer model's `stability_limit` allows
canopy_stability_limits(buffers, boundary_layer_model) = inverse_obukhov_length_limits(boundary_layer_model.stability_limit,
    boundary_layer_model.stability_functions,
    last(buffers.wind.above_canopy_profile.heights) - buffers.wind.displacement_height, buffers.wind.roughness_length)
