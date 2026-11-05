"""
    CanopyResidualForcing(; atmospheric_pressure, reference_temperature, reference_wind_speed,
        reference_humidity, incoming_longwave, ground_temperature, ground_emissivity,
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
    ground_temperature::typeof(1.0u"K")
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
    incoming_longwave, ground_temperature, ground_emissivity, ground_relative_humidity, leaf_water_potential,
    dry_conductance::C, canopy_top_air_temperature, canopy_top_relative_humidity, stability_limits) where C =
    CanopyResidualForcing{C}(atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity,
        incoming_longwave, ground_temperature, ground_emissivity, ground_relative_humidity, leaf_water_potential,
        dry_conductance, canopy_top_air_temperature, canopy_top_relative_humidity, stability_limits)

"""
    canopy_state_length(buffers)

How many values the canopy solve finds: a leaf temperature (K), air temperature (K) and
vapor density (g/m³) for each layer, then one inverse Obukhov length (1/m).
"""
canopy_state_length(buffers) = 3 * length(buffers.leaf.leaf_temperature) + 1

"""
    canopy_state_rows(layer_count)

Where each kind of value sits in the canopy state, as
`(; leaf_temperature, air_temperature, vapor_density, inverse_obukhov_length)`.
"""
canopy_state_rows(layer_count) = (; leaf_temperature=1:layer_count, air_temperature=(layer_count + 1):2layer_count,
    vapor_density=(2layer_count + 1):3layer_count, inverse_obukhov_length=3layer_count + 1)

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
    air_temperature = air.air_temperature_prev[layer]
    (; leaf_area, aerodynamic_area) = layer_leaf_areas(buffers, layer)
    conductance = blend_stomatal_conductance(dry_conductance, wet_canopy_fraction(model, buffers, layer))
    balance = leaf_heat_balance(leaf_temperature, buffers.leaf.absorbed_radiation[layer], air_temperature, air.vapor_density_prev[layer],
        buffers.wind.wind_speed[layer], atmospheric_pressure, model.leaf_parameters.leaf_emissivity, conductance,
        leaf_water_potential, buffers.leaf.leaf_body, leaf_area, model.leaf_convection_model; aerodynamic_area)
    return (; balance, leaf_area)
end

# A value inside its bounds. NaN, from a failed step, goes to the lower bound.
_within_bounds(value, lower_bound, upper_bound=typemax(value)) =
    isnan(value) ? lower_bound : clamp(value, lower_bound, upper_bound)

"""
    evaluate_canopy!(residual, state, buffers, model, boundary_layer_model, forcing::CanopyResidualForcing;
        site, vapour_pressure_equation)

Measures how far `state` (see [`canopy_state_length`](@ref)) is from solving the canopy. The
hour is solved when every value written to `residual` is close to zero:

- each leaf's net energy balance per unit leaf area (W/m²). A layer without leaves
  has no temperature of its own, so it gives leaf minus air temperature (K)
- each layer's air temperature, minus what the air-profile model gets from the leaves (K)
- each layer's vapor density, minus what the air-profile model gets from the leaves (g/m³)
- the inverse Obukhov length, minus what the canopy's heat flux implies (1/m)

The wind, radiation, leaf fluxes and air profile are left in `buffers`.

Returns `(; longwave_result, canopy_top_air_temperature,
canopy_top_relative_humidity, friction_velocity, inverse_obukhov_length)`.
"""
function evaluate_canopy!(residual, state, buffers, model::MultilayerCanopy, boundary_layer_model,
    forcing::CanopyResidualForcing; site, vapour_pressure_equation,
)
    (; atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity, incoming_longwave,
       ground_temperature, ground_emissivity, ground_relative_humidity) = forcing
    (; leaf_temperature, sensible_heat_source, evaporation_mass_flow, absorbed_radiation, net_balance) = buffers.leaf
    air = buffers.air_profile
    (; absorbed_shortwave) = buffers.shortwave
    (; absorbed_longwave) = buffers.longwave
    n_layers = length(leaf_temperature)
    rows = canopy_state_rows(n_layers)

    # A state out of bounds is evaluated at the bounds, and its residual is NaN
    lowest_temperature, highest_temperature = CANOPY_TEMPERATURE_BOUNDS
    in_bounds = true
    for i in 1:n_layers
        leaf = state[rows.leaf_temperature[i]] * u"K"
        air_temperature = state[rows.air_temperature[i]] * u"K"
        vapor_density = state[rows.vapor_density[i]] * u"g/m^3"
        # A negative vapor density is evaluated at zero. Its residual still points back to positive.
        in_bounds &= lowest_temperature <= leaf <= highest_temperature &&
            lowest_temperature <= air_temperature <= highest_temperature && !isnan(vapor_density)
        leaf_temperature[i] = _within_bounds(leaf, lowest_temperature, highest_temperature)
        air.air_temperature_prev[i] = _within_bounds(air_temperature, lowest_temperature, highest_temperature)
        air.vapor_density_prev[i] = _within_bounds(vapor_density, zero(vapor_density))
    end
    inverse_obukhov_length = state[rows.inverse_obukhov_length] * u"m^-1"
    in_bounds &= !isnan(inverse_obukhov_length)

    environment_instant = (; atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity)
    wind = canopy_wind_profile!(buffers, model, boundary_layer_model;
        site, environment_instant, canopy_source_temperature=leaf_temperature[1], vapour_pressure_equation,
        inverse_obukhov_length=isnan(inverse_obukhov_length) ? zero(inverse_obukhov_length) : inverse_obukhov_length)
    wind_speed = buffers.wind.wind_speed
    longwave_result = canopy_longwave!(buffers, model;
        leaf_temperature, ground_temperature, ground_emissivity, incoming_longwave)

    for i in 1:n_layers
        air_temperature = air.air_temperature_prev[i]
        (; exchange_fraction) = layer_leaf_areas(buffers, i)
        if iszero(exchange_fraction)
            sensible_heat_source[i] = zero(sensible_heat_source[i])
            evaporation_mass_flow[i] = zero(evaporation_mass_flow[i])
            absorbed_radiation[i] = zero(absorbed_radiation[i])
            net_balance[i] = zero(net_balance[i])
            residual[rows.leaf_temperature[i]] = ustrip(u"K", leaf_temperature[i] - air_temperature)
            continue
        end
        # absorbed_shortwave/absorbed_longwave are per unit ground area
        absorbed_radiation[i] = (absorbed_shortwave[i] + absorbed_longwave[i]) / (2.0 * exchange_fraction)
        (; balance, leaf_area) = layer_leaf_balance(buffers, model, forcing, i, leaf_temperature[i])
        sensible_heat_source[i] = balance.conv.convection_flow / 1.0u"m^2"
        evaporation_mass_flow[i] = balance.evap.transpiration_mass_flow
        net_balance[i] = balance.net
        residual[rows.leaf_temperature[i]] = ustrip(u"W/m^2", balance.net / leaf_area)
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
        boundary.top_vapour_density >= zero(boundary.top_vapour_density)
    canopy_top_air_temperature = clamp(boundary.top_temperature, lowest_temperature, highest_temperature)
    # The air profile takes vapor density, not relative humidity: humidity capped at 100%
    # would lose any extra vapor from the water balance
    canopy_top_vapor_density = max(boundary.top_vapour_density, zero(boundary.top_vapour_density))
    top_saturation = wet_air_properties(canopy_top_air_temperature, 1.0, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    canopy_top_relative_humidity = clamp(canopy_top_vapor_density / top_saturation, 0.0, 1.0)

    air_forcing = CanopyAirProfileForcing(; friction_velocity=wind.friction_velocity,
        inverse_obukhov_length=wind.inverse_obukhov_length, wind_speed, canopy_top_air_temperature,
        canopy_top_vapor_density, ground_temperature, ground_relative_humidity, sensible_heat_source,
        evaporation_mass_flow, leaf_temperature, atmospheric_pressure)
    canopy_heat_vapor_transport!(air, model.air_profile_model, boundary_layer_model, air_forcing;
        site, vapour_pressure_equation)

    for i in 1:n_layers
        residual[rows.air_temperature[i]] = ustrip(u"K", air.air_temperature_prev[i] - air.air_temperature[i])
        residual[rows.vapor_density[i]] = state[rows.vapor_density[i]] - ustrip(u"g/m^3", air.vapor_density[i])
    end
    # The stability the whole canopy's sensible flux through its top implies; the Obukhov
    # relation takes the downward flux
    ρcpTκg = calc_ρ_cp(reference_temperature) * reference_temperature / (boundary_layer_model.karman_constant * g_n)
    inverse_obukhov_length_update = limit_inverse_obukhov_length(stability_from_fluxes(InverseObukhovLength(), ρcpTκg,
        wind.friction_velocity, -boundary.sensible_flux), forcing.stability_limits)
    residual[rows.inverse_obukhov_length] = ustrip(u"m^-1", wind.inverse_obukhov_length - inverse_obukhov_length_update)
    in_bounds || fill!(residual, NaN)

    return (; longwave_result, canopy_top_air_temperature, canopy_top_relative_humidity,
        friction_velocity=wind.friction_velocity, inverse_obukhov_length=wind.inverse_obukhov_length)
end

# The canopy-top inverse Obukhov lengths the boundary layer model's `stability_limit` allows
canopy_stability_limits(buffers, boundary_layer_model) = inverse_obukhov_length_limits(boundary_layer_model.stability_limit,
    boundary_layer_model.stability_functions,
    last(buffers.wind.above_canopy_profile.heights) - buffers.wind.displacement_height, buffers.wind.roughness_length)
