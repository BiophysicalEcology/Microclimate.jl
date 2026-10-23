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
    CANOPY_TEMPERATURE_BOUNDS

Leaf, air and canopy-top temperatures outside which [`evaluate_canopy!`](@ref)
returns a NaN residual, as it does for negative vapor density.
"""
const CANOPY_TEMPERATURE_BOUNDS = (100.0u"K", 1000.0u"K")

# Floor resistance on the canopy-top boundary's ground path. The profile's own
# ground resistance is only known inside its transport, so ProfileGroundResistance
# adds none.
_canopy_top_floor_resistance(::ProfileGroundResistance, boundary_layer_model, forcing::GroundExchangeForcing) =
    (; heat=0.0u"s/m", vapor=0.0u"s/m")
_canopy_top_floor_resistance(coupling::AbstractGroundCoupling, boundary_layer_model, forcing::GroundExchangeForcing) =
    ground_resistance(coupling, boundary_layer_model, forcing)

"""
    layer_leaf_balance(buffers, model::MultilayerCanopy, forcing::CanopyResidualForcing, layer, leaf_temperature;
        vapour_pressure_equation)

Layer `layer`'s `leaf_heat_balance` at `leaf_temperature`, as `(; balance, leaf_area)`.
"""
function layer_leaf_balance(buffers, model::MultilayerCanopy, forcing::CanopyResidualForcing, layer, leaf_temperature;
    vapour_pressure_equation,
)
    (; atmospheric_pressure, leaf_water_potential, dry_conductance) = forcing
    air = buffers.air_profile
    air_temperature = air.air_temperature_prev[layer]
    exchange_fraction = 1.0 - buffers.longwave.layer_transmission[layer]
    layer_plant_area_index_boundaries = buffers.shortwave.layer_plant_area_index_boundaries
    saturation = wet_air_properties(air_temperature, 1.0, atmospheric_pressure; vapour_pressure_equation).vapour_density
    relative_humidity = clamp(air.vapor_density_prev[layer] / saturation, 0.0, 1.0)
    leaf_area = 2.0 * exchange_fraction * 1.0u"m^2"
    layer_pai = layer_plant_area_index_boundaries[layer + 1] - layer_plant_area_index_boundaries[layer]
    aerodynamic_area = 2.0 * layer_pai * 1.0u"m^2"
    conductance = blend_stomatal_conductance(dry_conductance, wet_canopy_fraction(model, buffers, layer))
    balance = leaf_heat_balance(leaf_temperature, buffers.leaf.absorbed_radiation[layer], air_temperature, relative_humidity,
        buffers.wind.wind_speed[layer], atmospheric_pressure, model.leaf_parameters.leaf_emissivity, conductance,
        leaf_water_potential, buffers.leaf.leaf_body, leaf_area, model.leaf_convection_model; aerodynamic_area)
    return (; balance, leaf_area)
end

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
       ground_temperature, ground_emissivity, ground_relative_humidity) = forcing
    (; leaf_temperature, sensible_heat_source, evaporation_mass_flow, absorbed_radiation, net_balance) = buffers.leaf
    air = buffers.air_profile
    (; absorbed_shortwave) = buffers.shortwave
    (; absorbed_longwave, layer_transmission) = buffers.longwave
    n_layers = length(leaf_temperature)

    # Evaluate out-of-bounds states at the bounds; their residual is NaN
    lo, hi = CANOPY_TEMPERATURE_BOUNDS
    in_bounds = true
    for i in 1:n_layers
        leaf = state[i] * u"K"
        air_temperature = state[n_layers + i] * u"K"
        vapor_density = state[2n_layers + i] * u"g/m^3"
        in_bounds &= lo <= leaf <= hi && lo <= air_temperature <= hi && vapor_density >= zero(vapor_density)
        leaf_temperature[i] = clamp(leaf, lo, hi)
        air.air_temperature_prev[i] = clamp(air_temperature, lo, hi)
        air.vapor_density_prev[i] = max(vapor_density, zero(vapor_density))
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
        # absorbed_shortwave/absorbed_longwave are per unit ground area
        absorbed_radiation[i] = (absorbed_shortwave[i] + absorbed_longwave[i]) / (2.0 * exchange_fraction)
        (; balance, leaf_area) = layer_leaf_balance(buffers, model, forcing, i, leaf_temperature[i]; vapour_pressure_equation)
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
    kinematic_viscosity = dry_air_properties(reference_temperature, atmospheric_pressure).kinematic_viscosity
    enhancement = turbulent_enhancement(model.roughness_sublayer, buffers.roughness_sublayer.enhancement,
        model.canopy_height, wind.canopy_top_wind_speed, kinematic_viscosity)
    column = AboveCanopyColumn(; canopy_height=model.canopy_height, displacement_height=buffers.wind.displacement_height,
        friction_velocity=wind.friction_velocity, inverse_obukhov_length=wind.inverse_obukhov_length, enhancement)
    floor = _canopy_top_floor_resistance(model.ground_coupling, boundary_layer_model, GroundExchangeForcing(;
        layer_height=buffers.wind.ground_layer_height, wind_speed=wind_speed[n_layers],
        layer_air_temperature=air.air_temperature_prev[n_layers], ground_temperature,
        roughness_height=site.roughness_height, friction_velocity=wind.friction_velocity,
        plant_area_index=buffers.air_profile.plant_area_index, atmospheric_pressure))
    boundary = canopy_top_flux_boundary(boundary_layer_model, column, buffers.wind.above_canopy_profile.heights[2],
        sum(sensible_heat_source), sum(evaporation_mass_flow) / 1.0u"m^2",
        ground_temperature, ground_vapor_density, reference_temperature, reference_vapor_density,
        forcing.canopy_top_air_temperature, initial_top_vapor_density; floor_resistance=floor)
    in_bounds &= lo <= boundary.top_temperature <= hi
    canopy_top_air_temperature = clamp(boundary.top_temperature, lo, hi)
    top_saturation = wet_air_properties(canopy_top_air_temperature, 1.0, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    canopy_top_relative_humidity = clamp(boundary.top_vapour_density / top_saturation, 0.0, 1.0)

    air_forcing = CanopyAirProfileForcing(; above_canopy_column=column, wind_speed, canopy_top_air_temperature,
        canopy_top_relative_humidity, ground_temperature, ground_relative_humidity, sensible_heat_source,
        evaporation_mass_flow, leaf_temperature, atmospheric_pressure)
    canopy_heat_vapor_transport!(air, model.air_profile_model, boundary_layer_model, air_forcing;
        site, vapour_pressure_equation, ground_coupling=model.ground_coupling)

    for i in 1:n_layers
        residual[n_layers + i] = ustrip(u"K", air.air_temperature_prev[i] - air.air_temperature[i])
        residual[2n_layers + i] = state[2n_layers + i] - ustrip(u"g/m^3", air.vapor_density[i])
    end
    residual[3n_layers + 1] = ustrip(u"m^-1", wind.inverse_obukhov_length - wind.inverse_obukhov_length_update)
    in_bounds || fill!(residual, NaN)

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
