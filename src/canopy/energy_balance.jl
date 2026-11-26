"""
    CANOPY_TIMESTEP

The canopy's one-hour step.
"""
const CANOPY_TIMESTEP = 1.0u"hr"

"""
    CanopyEnergyBalanceInputs(model::MultilayerCanopy; environment_instant, zenith_angle,
                               direct_horizontal_irradiance, diffuse_horizontal_irradiance,
                               ground_reflectance, ground_temperature, ground_emissivity,
                               ground_wetness, canopy_source_temperature, rainfall=0.0u"kg/m^2",
                               leaf_water_potential=model.leaf_parameters.leaf_water_potential)

Hourly forcing and lagged state for [`canopy_energy_balance!`](@ref). `ground_wetness` is the
wet fraction of the soil surface, which evaporates as the soil's own evaporation does.
"""
struct CanopyEnergyBalanceInputs{EI,ZA,DHI,DFI,GR,GT,GE,GW,CST,RF,LWP}
    environment_instant::EI
    zenith_angle::ZA
    direct_horizontal_irradiance::DHI
    diffuse_horizontal_irradiance::DFI
    ground_reflectance::GR
    ground_temperature::GT
    ground_emissivity::GE
    ground_wetness::GW
    canopy_source_temperature::CST
    rainfall::RF
    leaf_water_potential::LWP
end

function CanopyEnergyBalanceInputs(model::MultilayerCanopy;
    environment_instant, zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance,
    ground_reflectance, ground_temperature, ground_emissivity, ground_wetness, canopy_source_temperature,
    rainfall=0.0u"kg/m^2", leaf_water_potential=model.leaf_parameters.leaf_water_potential,
)
    CanopyEnergyBalanceInputs(environment_instant, zenith_angle, direct_horizontal_irradiance,
        diffuse_horizontal_irradiance, ground_reflectance, ground_temperature, ground_emissivity,
        ground_wetness, canopy_source_temperature, rainfall, leaf_water_potential)
end

# Bounds wet-surface evaporation, or condensation, to what storage can supply or hold, and turns any latent heat that can't happen into sensible heat
function _cap_wet_surface_evaporation(evaporation_mass_flux, sensible_heat_flux, wet_surface_evaporation_mass_flux,
    available_mass, layer_capacity, air_temperature,
)
    demanded_mass = uconvert(u"kg/m^2", wet_surface_evaporation_mass_flux * CANOPY_TIMESTEP)
    actual_mass = clamp(demanded_mass, available_mass - layer_capacity, available_mass)
    shortfall_mass = demanded_mass - actual_mass
    iszero(shortfall_mass) && return (evaporation_mass_flux, sensible_heat_flux, actual_mass)

    shortfall_mass_flux = uconvert(unit(evaporation_mass_flux), shortfall_mass / CANOPY_TIMESTEP)
    shortfall_heat_flux = uconvert(unit(sensible_heat_flux), shortfall_mass_flux * enthalpy_of_vaporisation(air_temperature))
    return (evaporation_mass_flux - shortfall_mass_flux, sensible_heat_flux + shortfall_heat_flux, actual_mass)
end

# The ground stays at the temperature it was given
keep_ground_temperature(ground_temperature, ground_shortwave, ground_absorbed_longwave) = ground_temperature

"""
    canopy_energy_balance!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs;
        site, vapour_pressure_equation=GoffGratch(), soil_hydraulic_model=nothing,
        ground_temperature_after_hour=keep_ground_temperature)

Solves the canopy for one hour. Shortwave and rain interception don't depend on leaf
temperature, so they are computed once. Leaf temperature, air temperature, vapor density,
ground temperature, longwave and stability are then solved together by [`converge_canopy!`](@ref).

`ground_temperature_after_hour(ground_temperature, ground_shortwave, ground_absorbed_longwave)`
returns the ground temperature the soil reaches by the end of the hour, with `ground_shortwave` the
[`canopy_shortwave!`](@ref) result. Any soil solve can be used.

The wet part of a leaf ([`wet_canopy_fraction`](@ref)) evaporates as free water
([`WET_SURFACE_EVAPORATION`](@ref)), and the rest transpires. The water that
evaporates or condenses is taken from storage once, after the solve. It is limited by
[`available_canopy_water`](@ref) and [`canopy_water_capacity`](@ref), and any latent heat that
can't happen is added to the sensible heat, without solving leaf temperature again.

`inputs.ground_temperature` is where the ground starts. `inputs.canopy_source_temperature` and
`inputs.leaf_water_potential` come from the caller and stay fixed during the solve.

Returns `(; ground_absorbed_shortwave, net_absorbed_shortwave, ground_absorbed_longwave,
net_absorbed_longwave, ground_throughfall, canopy_potential_transpiration,
ground_heat_conductance, ground_vapor_conductance, iterations, canopy_top_air_temperature,
canopy_top_relative_humidity, friction_velocity, obukhov_length, displacement_height)`.

- `canopy_potential_transpiration` is what the canopy would transpire with no water stress,
  which the soil hydraulics model (e.g. `CampbellSoilHydraulics`) uses as demand.
- `ground_heat_conductance` and `ground_vapor_conductance` connect the ground to the lowest
  layer. The soil's surface exchange uses them (see `ground_convection_conditions`).
"""
function canopy_energy_balance!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs;
    site, vapour_pressure_equation=GoffGratch(), soil_hydraulic_model=nothing,
    ground_temperature_after_hour::G=keep_ground_temperature,
) where G
    (; environment_instant, zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance,
       ground_reflectance, ground_temperature, ground_emissivity, canopy_source_temperature,
       rainfall, leaf_water_potential) = inputs

    leaf_temperature = buffers.leaf.leaf_temperature
    air_temperature = buffers.air_profile.air_temperature
    evaporation_mass_flux = buffers.leaf.evaporation_mass_flux
    sensible_heat_flux = buffers.leaf.sensible_heat_flux
    potential_evaporation_mass_flux = buffers.leaf.potential_evaporation_mass_flux
    leaf_dew_formed = buffers.leaf.leaf_dew_formed
    leaf_frost_formed = buffers.leaf.leaf_frost_formed
    leaf_body = buffers.leaf.leaf_body
    (; wind_speed) = buffers.wind
    atmospheric_pressure = environment_instant.atmospheric_pressure
    n_layers = length(leaf_temperature)

    shortwave_result = canopy_shortwave!(buffers, model;
        zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance, ground_reflectance)

    wind_result = canopy_wind_profile!(buffers, model, boundary_layer_model;
        site, environment_instant, canopy_source_temperature, vapour_pressure_equation)

    interception_result = canopy_interception!(buffers, model; rainfall, wind_speed)

    # The first hour starts from the canopy-top air, later hours from the previous hour's answer
    if !buffers.has_previous_hour[]
        fill!(air_temperature, wind_result.canopy_top_air_temperature)
        fill!(leaf_temperature, wind_result.canopy_top_air_temperature)
        fill!(buffers.air_profile.vapor_density, wet_air_properties(wind_result.canopy_top_air_temperature,
            wind_result.canopy_top_relative_humidity, atmospheric_pressure; vapour_pressure_equation).vapour_density)
    end

    longwave_result, iterations, canopy_top_air_temperature, canopy_top_relative_humidity, friction_velocity, inverse_obukhov_length =
        converge_canopy!(buffers, model, boundary_layer_model, inputs, wind_result, model.convergence_model;
            site, vapour_pressure_equation, soil_hydraulic_model,
            ground_temperature_after_hour=(ground_temperature, ground_absorbed_longwave) -> ground_temperature_after_hour(
                ground_temperature, shortwave_result, ground_absorbed_longwave))
    buffers.has_previous_hour[] = true

    melt_drip_total = 0.0u"kg/m^2"
    for layer in 1:n_layers
        # Only wet-surface evaporation draws down intercepted rain; transpiration draws from the soil
        available_mass = available_canopy_water(model, buffers, layer)
        layer_capacity = canopy_water_capacity(model, buffers, layer)
        evaporation_mass_flux[layer], sensible_heat_flux[layer], actual_mass = _cap_wet_surface_evaporation(
            evaporation_mass_flux[layer], sensible_heat_flux[layer], buffers.leaf.wet_surface_evaporation_mass_flux[layer],
            available_mass, layer_capacity, air_temperature[layer])

        layer_leaf_temperature = leaf_temperature[layer]
        leaf_dew_formed[layer] = 0.0u"kg/m^2"
        leaf_frost_formed[layer] = 0.0u"kg/m^2"
        # Standing frost melts into standing dew (and leaf_surface_water) once the leaf
        # warms above 0°C; any overflow past capacity drips to the ground (stemflow).
        standing_frost_layer = leaf_standing_frost(model, buffers, layer)
        if standing_frost_layer > 0.0u"kg/m^2" && layer_leaf_temperature > u"K"(0.0u"°C")
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
        # tied to evaporation_mass_flux/sensible_heat_flux above.
        (; leaf_area, aerodynamic_area) = layer_leaf_areas(buffers, layer)
        convection = leaf_convection(model.leaf_convection_model, leaf_body, aerodynamic_area, air_temperature[layer],
            layer_leaf_temperature, wind_speed[layer], atmospheric_pressure)
        dew_frost_energy_flux = condensation_energy_flux(model.condensation_model;
            leaf_temperature=layer_leaf_temperature, air_temperature=air_temperature[layer],
            relative_humidity=buffers.air_profile.relative_humidity[layer],
            absorbed_radiation=buffers.leaf.absorbed_radiation[layer],
            leaf_emissivity=model.leaf_parameters.leaf_emissivity,
            heat_transfer_coefficient=convection.heat_transfer_coefficient.combined,
            leaf_area, aerodynamic_area, vapour_pressure_equation,
        )
        condensing = dew_frost_energy_flux < 0.0u"W/m^2"

        # actual_mass < 0 is condensation; frost isn't part of leaf_surface_water
        # (tracked, uncapped, in leaf_standing_frost instead, like ground).
        if actual_mass < 0.0u"kg/m^2" && condensing && layer_leaf_temperature <= u"K"(0.0u"°C")
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

    canopy_potential_transpiration = sum(potential_evaporation_mass_flux)
    ground_heat_conductance = buffers.air_profile.ground_heat_conductance[]
    ground_vapor_conductance = buffers.air_profile.ground_vapor_conductance[]

    return (; shortwave_result..., longwave_result..., interception_result...,
        ground_throughfall = interception_result.ground_throughfall + melt_drip_total,
        canopy_potential_transpiration,
        ground_heat_conductance, ground_vapor_conductance, iterations,
        canopy_top_air_temperature, canopy_top_relative_humidity,
        friction_velocity, obukhov_length=inv(inverse_obukhov_length),
        displacement_height=buffers.wind.displacement_height)
end

"""
    converge_canopy!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs,
        initial_wind, convergence_model; site, vapour_pressure_equation, soil_hydraulic_model,
        ground_temperature_after_hour)

Solves the hour's canopy state with `convergence_model` ([`AndersonCanopyConvergence`](@ref)).
`initial_wind` is the hour's [`canopy_wind_profile!`](@ref) result, which gives the canopy-top
air temperature and humidity the hour starts from. Returns `(longwave_result, iterations,
canopy_top_air_temperature, canopy_top_relative_humidity, friction_velocity,
inverse_obukhov_length)` at the solution.
"""
function converge_canopy! end
