"""
    CanopyResidualForcing(; atmospheric_pressure, reference_temperature, reference_wind_speed,
        reference_humidity, incoming_longwave, ground_temperature, ground_emissivity,
        ground_relative_humidity, leaf_water_potential, dry_conductance,
        canopy_top_air_temperature, canopy_top_relative_humidity)

The hour's forcing for [`canopy_residual!`](@ref): everything the residual
reads that doesn't depend on canopy state. Concretely typed so the solver
cache built by `allocate_canopy` can hold it, and overwritten each hour.
`canopy_top_air_temperature`/`canopy_top_relative_humidity` are the hour's
initial canopy-top state, the anchor of [`canopy_top_flux_boundary`](@ref).
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
end
# Convert each field to its fixed unit, inferring only the conductance type
CanopyResidualForcing(atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity,
    incoming_longwave, ground_temperature, ground_emissivity, ground_relative_humidity, leaf_water_potential,
    dry_conductance::C, canopy_top_air_temperature, canopy_top_relative_humidity) where C =
    CanopyResidualForcing{C}(atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity,
        incoming_longwave, ground_temperature, ground_emissivity, ground_relative_humidity, leaf_water_potential,
        dry_conductance, canopy_top_air_temperature, canopy_top_relative_humidity)

"""
    canopy_state_length(buffers)

Length of the canopy state vector: leaf temperature (K), air temperature (K)
and vapor density (g/m³) per layer, then the inverse Obukhov length (1/m).
"""
canopy_state_length(buffers) = 3 * length(buffers.leaf.leaf_temperature) + 1

"""
    evaluate_canopy!(residual, state, buffers, model, boundary_layer_model, forcing::CanopyResidualForcing;
        site, vapour_pressure_equation)

Write the canopy residual for `state` (see [`canopy_state_length`](@ref))
into `residual`, and every intermediate (wind, longwave, leaf sources, air
profile) into `buffers`:

- leaf net energy balance per unit leaf area (W/m²); a layer without leaf
  area has no leaf temperature of its own, so it gives `T_leaf - T_air` (K)
- air temperature minus the air-profile model's transport of the leaf
  sources (K)
- vapor density minus its transport (g/m³)
- inverse Obukhov length minus what the Obukhov relation gives back for it (1/m)

Returns `(; longwave_result, canopy_top_air_temperature,
canopy_top_relative_humidity, friction_velocity, inverse_obukhov_length)`.
"""
function evaluate_canopy!(residual, state, buffers, model::MultilayerCanopy, boundary_layer_model,
    forcing::CanopyResidualForcing; site, vapour_pressure_equation,
)
    (; atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity, incoming_longwave,
       ground_temperature, ground_emissivity, ground_relative_humidity, leaf_water_potential, dry_conductance) = forcing
    (; leaf_temperature, sensible_heat_source, evaporation_mass_flow, absorbed_radiation, net_balance, leaf_body) = buffers.leaf
    air = buffers.air_profile
    (; absorbed_shortwave, layer_plant_area_index_boundaries) = buffers.shortwave
    (; absorbed_longwave, layer_transmission) = buffers.longwave
    n_layers = length(leaf_temperature)

    for i in 1:n_layers
        leaf_temperature[i] = state[i] * u"K"
        air.air_temperature_prev[i] = state[n_layers + i] * u"K"
        air.vapor_density_prev[i] = state[2n_layers + i] * u"g/m^3"
    end

    environment_instant = (; atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity)
    wind = canopy_wind_profile!(buffers, model, boundary_layer_model;
        site, environment_instant, canopy_source_temperature=leaf_temperature[1], vapour_pressure_equation,
        inverse_obukhov_length=state[3n_layers + 1] * u"m^-1")
    wind_speed = buffers.wind.wind_speed
    longwave_result = canopy_longwave!(buffers, model;
        leaf_temperature, ground_temperature, ground_emissivity, incoming_longwave)

    for i in 1:n_layers
        air_temperature = air.air_temperature_prev[i]
        # Radiative exchange area is the layer's interception fraction (1-τ);
        # aerodynamic exchange scales with leaf surface (PAI). 2.0 = both sides.
        exchange_fraction = 1.0 - layer_transmission[i]
        if iszero(exchange_fraction)
            sensible_heat_source[i] = zero(sensible_heat_source[i])
            evaporation_mass_flow[i] = zero(evaporation_mass_flow[i])
            absorbed_radiation[i] = zero(absorbed_radiation[i])
            net_balance[i] = zero(net_balance[i])
            residual[i] = ustrip(u"K", leaf_temperature[i] - air_temperature)
            continue
        end
        saturation = wet_air_properties(air_temperature, 1.0, atmospheric_pressure; vapour_pressure_equation).vapour_density
        relative_humidity = clamp(air.vapor_density_prev[i] / saturation, 0.0, 1.0)
        leaf_area = 2.0 * exchange_fraction * 1.0u"m^2"
        layer_pai = layer_plant_area_index_boundaries[i + 1] - layer_plant_area_index_boundaries[i]
        aerodynamic_area = 2.0 * layer_pai * 1.0u"m^2"
        # absorbed_shortwave/absorbed_longwave are per unit ground area
        absorbed_radiation[i] = (absorbed_shortwave[i] + absorbed_longwave[i]) / (2.0 * exchange_fraction)
        conductance = blend_stomatal_conductance(dry_conductance, wet_canopy_fraction(model, buffers, i))
        balance = leaf_heat_balance(leaf_temperature[i], absorbed_radiation[i], air_temperature, relative_humidity,
            wind_speed[i], atmospheric_pressure, model.leaf_parameters.leaf_emissivity, conductance,
            leaf_water_potential, leaf_body, leaf_area, model.leaf_convection_model; aerodynamic_area)
        sensible_heat_source[i] = balance.conv.convection_flow / 1.0u"m^2"
        evaporation_mass_flow[i] = balance.evap.transpiration_mass_flow
        net_balance[i] = balance.net
        residual[i] = ustrip(u"W/m^2", balance.net / leaf_area)
    end

    reference_vapor_density = wet_air_properties(reference_temperature, reference_humidity, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    ground_vapor_density = wet_air_properties(ground_temperature, ground_relative_humidity, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    initial_top_vapor_density = wet_air_properties(forcing.canopy_top_air_temperature, forcing.canopy_top_relative_humidity,
        atmospheric_pressure; vapour_pressure_equation).vapour_density
    boundary = canopy_top_flux_boundary(boundary_layer_model, model.canopy_height, buffers.wind.displacement_height,
        buffers.wind.roughness_length, buffers.wind.above_canopy_profile.heights[2], wind.friction_velocity,
        wind.inverse_obukhov_length, sum(sensible_heat_source), sum(evaporation_mass_flow) / 1.0u"m^2",
        ground_temperature, ground_vapor_density, reference_temperature, reference_vapor_density,
        forcing.canopy_top_air_temperature, initial_top_vapor_density)
    canopy_top_air_temperature = boundary.top_temperature
    top_saturation = wet_air_properties(canopy_top_air_temperature, 1.0, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    canopy_top_relative_humidity = clamp(boundary.top_vapour_density / top_saturation, 0.0, 1.0)

    air_forcing = CanopyAirProfileForcing(; friction_velocity=wind.friction_velocity,
        inverse_obukhov_length=wind.inverse_obukhov_length, wind_speed, canopy_top_air_temperature,
        canopy_top_relative_humidity, ground_temperature, ground_relative_humidity, sensible_heat_source,
        evaporation_mass_flow, leaf_temperature, atmospheric_pressure)
    canopy_heat_vapor_transport!(air, model.air_profile_model, boundary_layer_model, air_forcing; site, vapour_pressure_equation)

    for i in 1:n_layers
        residual[n_layers + i] = ustrip(u"K", air.air_temperature_prev[i] - air.air_temperature[i])
        residual[2n_layers + i] = ustrip(u"g/m^3", air.vapor_density_prev[i] - air.vapor_density[i])
    end
    residual[3n_layers + 1] = ustrip(u"m^-1", wind.inverse_obukhov_length - wind.inverse_obukhov_length_update)

    return (; longwave_result, canopy_top_air_temperature, canopy_top_relative_humidity,
        friction_velocity=wind.friction_velocity, inverse_obukhov_length=wind.inverse_obukhov_length)
end

"""
    canopy_residual!(residual, state, buffers, (model, boundary_layer_model, forcing, site, vapour_pressure_equation))

[`evaluate_canopy!`](@ref) in the `f!(y, x, contexts...)` form
DifferentiationInterface differentiates: `buffers` is the `Cache` context,
the tuple the `Constant` one, with `forcing` a `Ref` overwritten each hour.
"""
function canopy_residual!(residual, state, buffers, (model, boundary_layer_model, forcing, site, vapour_pressure_equation))
    evaluate_canopy!(residual, state, buffers, model, boundary_layer_model, forcing[]; site, vapour_pressure_equation)
    return nothing
end
