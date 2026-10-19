"""
    CANOPY_TIMESTEP

Fixed one-hour step for rain-interception storage depletion, matching
`MicroModel`'s hourly forcing cadence.
"""
const CANOPY_TIMESTEP = 1.0u"hr"

"""
    CanopyEnergyBalanceInputs(model::MultilayerCanopy; environment_instant, zenith_angle,
                               direct_horizontal_irradiance, diffuse_horizontal_irradiance,
                               ground_reflectance, ground_temperature, ground_emissivity,
                               ground_relative_humidity, canopy_source_temperature, rainfall=0.0u"kg/m^2",
                               leaf_water_potential=model.leaf_parameters.leaf_water_potential)

Hourly forcing and lagged state for [`canopy_energy_balance!`](@ref).
"""
struct CanopyEnergyBalanceInputs{EI,ZA,DHI,DFI,GR,GT,GE,GRH,CST,RF,LWP}
    environment_instant::EI
    zenith_angle::ZA
    direct_horizontal_irradiance::DHI
    diffuse_horizontal_irradiance::DFI
    ground_reflectance::GR
    ground_temperature::GT
    ground_emissivity::GE
    ground_relative_humidity::GRH
    canopy_source_temperature::CST
    rainfall::RF
    leaf_water_potential::LWP
end

function CanopyEnergyBalanceInputs(model::MultilayerCanopy;
    environment_instant, zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance,
    ground_reflectance, ground_temperature, ground_emissivity, ground_relative_humidity, canopy_source_temperature,
    rainfall=0.0u"kg/m^2", leaf_water_potential=model.leaf_parameters.leaf_water_potential,
)
    CanopyEnergyBalanceInputs(environment_instant, zenith_angle, direct_horizontal_irradiance,
        diffuse_horizontal_irradiance, ground_reflectance, ground_temperature, ground_emissivity,
        ground_relative_humidity, canopy_source_temperature, rainfall, leaf_water_potential)
end

# Bounds the wet-surface share of evaporation_mass_flow to what storage can
# supply (positive demand) or absorb (negative demand, i.e. condensation);
# reroutes any unrealized latent flux to sensible so mass/energy stay exact
# either direction. Dew/frost classification of the result happens below,
# in the per-layer loop's condensation_model call.
function _cap_wet_surface_evaporation(evaporation_mass_flow, sensible_heat_flow, wet_fraction, available_mass,
    layer_capacity, air_temperature,
)
    demanded_mass = uconvert(u"kg", evaporation_mass_flow * wet_fraction * CANOPY_TIMESTEP) / 1.0u"m^2"
    actual_mass = clamp(demanded_mass, available_mass - layer_capacity, available_mass)
    shortfall_mass = demanded_mass - actual_mass
    iszero(shortfall_mass) && return (evaporation_mass_flow, sensible_heat_flow, actual_mass)

    shortfall_flow = uconvert(unit(evaporation_mass_flow), shortfall_mass * 1.0u"m^2" / CANOPY_TIMESTEP)
    shortfall_power = uconvert(unit(sensible_heat_flow), shortfall_flow * enthalpy_of_vaporisation(air_temperature) / 1.0u"m^2")
    return (evaporation_mass_flow - shortfall_flow, sensible_heat_flow + shortfall_power, actual_mass)
end

"""
    canopy_energy_balance!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs;
        site, vapour_pressure_equation=GoffGratch(), soil_hydraulic_model=nothing)

Hourly solve tying every canopy sub-model together. Shortwave, the
above-canopy wind/temperature profile, and rain interception are computed
once per call (none depends on leaf temperature). Longwave, the in-canopy
air profile, and per-layer leaf temperature are then solved together by
[`converge_canopy!`](@ref), dispatching on `model.convergence_model`
([`NewtonCanopyConvergence`](@ref)).

A wetted leaf surface ([`wet_canopy_fraction`](@ref)) blends stomatal
conductance toward [`WET_SURFACE_CONDUCTANCE`](@ref)
([`blend_stomatal_conductance`](@ref)); storage depletion (or, for
condensation, refill) from that evaporation is applied once, after
convergence, bounded to [`available_canopy_water`](@ref)/
[`canopy_water_capacity`](@ref) with any unrealized latent flux rerouted to
sensible -- a conservative closure, not a re-solved leaf temperature.

`inputs.ground_temperature`/`canopy_source_temperature`/`leaf_water_potential`
are lagged (previous-hour) boundary conditions; only canopy-internal state is
iterated. The `leaf_water_potential` lag is bounded by how fast it itself
moves (a soil-moisture function, day-timescale) and self-corrects next hour.

Per-layer relative humidity is resolved by `model.air_profile_model` from
each layer's own evaporation source, the same way per-layer air temperature
is.

Returns `(; ground_absorbed_shortwave, net_absorbed_shortwave,
ground_absorbed_longwave, net_absorbed_longwave, ground_throughfall,
canopy_potential_transpiration, ground_heat_conductance, ground_vapor_conductance,
iterations, canopy_top_air_temperature, canopy_top_relative_humidity,
friction_velocity, inverse_obukhov_length, displacement_height)`. `canopy_potential_transpiration`
is the canopy-summed unstressed (water_potential=0) transpiration rate, for
feeding a soil-hydraulics demand term (e.g. `CampbellSoilHydraulics`).
`ground_heat_conductance`/`ground_vapor_conductance` are `model.air_profile_model`'s
own ground-to-lowest-layer conductances, for the soil model's surface exchange to
reuse instead of inverting the log law at a reference height too close to the
roughness length (see `ground_convection_conditions`). The last five fields are
the solved above-canopy boundary state —
diagnostic-oriented (e.g. comparing against an independent implementation's own
boundary solve), not consumed elsewhere in this package.
"""
function canopy_energy_balance!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs;
    site, vapour_pressure_equation=GoffGratch(), soil_hydraulic_model=nothing,
)
    (; environment_instant, zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance,
       ground_reflectance, ground_temperature, ground_emissivity, ground_relative_humidity, canopy_source_temperature,
       rainfall, leaf_water_potential) = inputs

    leaf_temperature = buffers.leaf.leaf_temperature
    air_temperature = buffers.air_profile.air_temperature
    evaporation_mass_flow = buffers.leaf.evaporation_mass_flow
    sensible_heat_source = buffers.leaf.sensible_heat_source
    potential_evaporation_mass_flow = buffers.leaf.potential_evaporation_mass_flow
    leaf_dew_formed = buffers.leaf.leaf_dew_formed
    leaf_frost_formed = buffers.leaf.leaf_frost_formed
    leaf_body = buffers.leaf.leaf_body
    (; wind_speed) = buffers.wind
    (; layer_plant_area_index_boundaries) = buffers.shortwave
    (; layer_transmission) = buffers.longwave
    atmospheric_pressure = environment_instant.atmospheric_pressure
    n_layers = length(leaf_temperature)

    shortwave_result = canopy_shortwave!(buffers, model;
        zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance, ground_reflectance)

    wind_result = canopy_wind_profile!(buffers, model, boundary_layer_model;
        site, environment_instant, canopy_source_temperature, vapour_pressure_equation)

    interception_result = canopy_interception!(buffers, model; rainfall, wind_speed)

    # Warm-start from last hour's solution; only seed on the very first call
    # (buffer still zero, and the guess must not be 0 K).
    iszero(air_temperature[1]) && fill!(air_temperature, wind_result.canopy_top_air_temperature)
    iszero(leaf_temperature[1]) && fill!(leaf_temperature, wind_result.canopy_top_air_temperature)

    longwave_result, iter, canopy_top_air_temperature, canopy_top_relative_humidity, friction_velocity, inverse_obukhov_length =
        converge_canopy!(buffers, model, boundary_layer_model, inputs, wind_result, model.convergence_model;
            site, vapour_pressure_equation, soil_hydraulic_model)

    melt_drip_total = 0.0u"kg/m^2"
    for layer in 1:n_layers
        # evaporation_mass_flow blends transpiration and wet-surface
        # evaporation; only the wet fraction should draw down intercepted
        # rain — transpiration draws from soil/plant water instead.
        wet_fraction = wet_canopy_fraction(model, buffers, layer)
        available_mass = available_canopy_water(model, buffers, layer)
        layer_capacity = canopy_water_capacity(model, buffers, layer)
        evaporation_mass_flow[layer], sensible_heat_source[layer], actual_mass = _cap_wet_surface_evaporation(
            evaporation_mass_flow[layer], sensible_heat_source[layer], wet_fraction, available_mass, layer_capacity,
            air_temperature[layer])

        leaf_temp = leaf_temperature[layer]
        leaf_dew_formed[layer] = 0.0u"kg/m^2"
        leaf_frost_formed[layer] = 0.0u"kg/m^2"
        # Standing frost melts into standing dew (and leaf_surface_water) once the leaf
        # warms above 0°C; any overflow past capacity drips to the ground (stemflow).
        standing_frost_layer = leaf_standing_frost(model, buffers, layer)
        if standing_frost_layer > 0.0u"kg/m^2" && leaf_temp > u"K"(0.0u"°C")
            add_leaf_standing_dew!(model, buffers, layer, standing_frost_layer)
            add_leaf_standing_frost!(model, buffers, layer, -standing_frost_layer)
            deplete_canopy_water!(model, buffers, layer, -standing_frost_layer)
            overflow = max(available_canopy_water(model, buffers, layer) - layer_capacity, 0.0u"kg/m^2")
            if overflow > 0.0u"kg/m^2"
                deplete_canopy_water!(model, buffers, layer, overflow)
                melt_drip_total += overflow
            end
        end

        # condensation_model classifies actual_mass<0 as dew or frost (or
        # neither); it does not supply the mass itself -- unlike ground,
        # leaf_surface_water/leaf_standing_dew/leaf_standing_frost aren't
        # independent reservoirs, so the only physically-consistent mass
        # is actual_mass, already fixed by _cap_wet_surface_evaporation and
        # tied to evaporation_mass_flow/sensible_heat_source above.
        exchange_fraction = 1.0 - layer_transmission[layer]
        leaf_area = 2.0 * exchange_fraction * 1.0u"m^2"
        layer_pai = layer_plant_area_index_boundaries[layer + 1] - layer_plant_area_index_boundaries[layer]
        aerodynamic_area = 2.0 * layer_pai * 1.0u"m^2"
        conv = leaf_convection(model.leaf_convection_model, leaf_body, aerodynamic_area, air_temperature[layer],
            leaf_temp, wind_speed[layer], atmospheric_pressure)
        dew_frost_energy_flux = condensation_energy_flux(model.condensation_model;
            leaf_temperature=leaf_temp, air_temperature=air_temperature[layer],
            relative_humidity=buffers.air_profile.relative_humidity[layer],
            absorbed_radiation=buffers.leaf.absorbed_radiation[layer],
            leaf_emissivity=model.leaf_parameters.leaf_emissivity,
            heat_transfer_coefficient=conv.heat_transfer_coefficient.combined,
            leaf_area, aerodynamic_area, vapour_pressure_equation,
        )
        condensing = dew_frost_energy_flux < 0.0u"W/m^2"

        # actual_mass < 0 is condensation; frost isn't part of leaf_surface_water
        # (tracked, uncapped, in leaf_standing_frost instead, like ground).
        if actual_mass < 0.0u"kg/m^2" && condensing && leaf_temp <= u"K"(0.0u"°C")
            leaf_frost_formed[layer] = -actual_mass
            add_leaf_standing_frost!(model, buffers, layer, leaf_frost_formed[layer])
        elseif actual_mass < 0.0u"kg/m^2" && condensing
            leaf_dew_formed[layer] = -actual_mass
            add_leaf_standing_dew!(model, buffers, layer, leaf_dew_formed[layer])
            deplete_canopy_water!(model, buffers, layer, actual_mass)
        elseif actual_mass < 0.0u"kg/m^2"
            # Bulk-transfer condenses but condensation_model disagrees (e.g. a
            # momentary vapour-density inversion despite net radiative gain) --
            # ordinary interception gain, not dew/frost.
            deplete_canopy_water!(model, buffers, layer, actual_mass)
        else
            frost_sublimated = min(actual_mass, leaf_standing_frost(model, buffers, layer))
            add_leaf_standing_frost!(model, buffers, layer, -frost_sublimated)
            deplete_canopy_water!(model, buffers, layer, actual_mass - frost_sublimated)
        end
        # Floor rather than exact molecular attribution: drip/evaporation may have
        # removed leaf_surface_water that was standing dew without us tracking which.
        clamp_leaf_standing_dew!(model, buffers, layer, available_canopy_water(model, buffers, layer))
    end

    canopy_potential_transpiration = sum(potential_evaporation_mass_flow) / 1.0u"m^2"
    ground_heat_conductance = buffers.air_profile.ground_heat_conductance[]
    ground_vapor_conductance = buffers.air_profile.ground_vapor_conductance[]

    return (; shortwave_result..., longwave_result..., interception_result...,
        ground_throughfall = interception_result.ground_throughfall + melt_drip_total,
        canopy_potential_transpiration,
        ground_heat_conductance, ground_vapor_conductance, iterations=iter,
        canopy_top_air_temperature, canopy_top_relative_humidity,
        friction_velocity, obukhov_length=inv(inverse_obukhov_length),
        displacement_height=buffers.wind.displacement_height)
end

"""
    converge_canopy!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs,
        initial_wind, convergence_model; site, vapour_pressure_equation, soil_hydraulic_model)

Solve the hour's leaf/air-temperature state per `convergence_model`
([`NewtonCanopyConvergence`](@ref)). `initial_wind` is the hour's
[`canopy_wind_profile!`](@ref) result, whose canopy-top air temperature and
humidity anchor the canopy-top boundary. Returns `(longwave_result,
iterations, canopy_top_air_temperature, canopy_top_relative_humidity,
friction_velocity, inverse_obukhov_length)` at the solution.
"""
function converge_canopy! end
