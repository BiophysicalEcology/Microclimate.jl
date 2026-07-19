"""
    CANOPY_TIMESTEP

Fixed one-hour step for rain-interception storage depletion, matching
`MicroModel`'s hourly forcing cadence.
"""
const CANOPY_TIMESTEP = 1.0u"hr"

"""
    CanopyEnergyBalanceInputs(model::MultilayerCanopy; site, environment_instant, zenith_angle,
                               direct_horizontal_irradiance, diffuse_horizontal_irradiance,
                               ground_reflectance, ground_temperature, ground_emissivity,
                               canopy_source_temperature, rainfall=0.0u"kg/m^2",
                               leaf_water_potential=model.leaf_parameters.leaf_water_potential,
                               vapour_pressure_equation=GoffGratch())

Per-hour inputs to [`canopy_energy_balance!`](@ref), mirroring
`SoilEnergyInputs`: built once, refreshed in place via
[`update_canopy_energy_balance_inputs!`](@ref).
"""
mutable struct CanopyEnergyBalanceInputs{S,EI,ZA,DHI,DFI,GR,GT,GE,CST,RF,LWP,VP}
    site::S
    environment_instant::EI
    zenith_angle::ZA
    direct_horizontal_irradiance::DHI
    diffuse_horizontal_irradiance::DFI
    ground_reflectance::GR
    ground_temperature::GT
    ground_emissivity::GE
    canopy_source_temperature::CST
    rainfall::RF
    leaf_water_potential::LWP
    vapour_pressure_equation::VP
end

function CanopyEnergyBalanceInputs(model::MultilayerCanopy;
    site, environment_instant, zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance,
    ground_reflectance, ground_temperature, ground_emissivity, canopy_source_temperature,
    rainfall=0.0u"kg/m^2", leaf_water_potential=model.leaf_parameters.leaf_water_potential,
    vapour_pressure_equation=GoffGratch(),
)
    CanopyEnergyBalanceInputs(site, environment_instant, zenith_angle, direct_horizontal_irradiance,
        diffuse_horizontal_irradiance, ground_reflectance, ground_temperature, ground_emissivity,
        canopy_source_temperature, rainfall, leaf_water_potential, vapour_pressure_equation)
end

"""
    update_canopy_energy_balance_inputs!(inputs::CanopyEnergyBalanceInputs; environment_instant, zenith_angle,
                                          direct_horizontal_irradiance, diffuse_horizontal_irradiance,
                                          ground_reflectance, ground_temperature, ground_emissivity,
                                          canopy_source_temperature, rainfall=0.0u"kg/m^2",
                                          leaf_water_potential=inputs.leaf_water_potential)

In-place per-hour update, no allocation. `site`/`vapour_pressure_equation`
aren't refreshed (static for a run).
"""
function update_canopy_energy_balance_inputs!(inputs::CanopyEnergyBalanceInputs;
    environment_instant, zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance,
    ground_reflectance, ground_temperature, ground_emissivity, canopy_source_temperature,
    rainfall=0.0u"kg/m^2", leaf_water_potential=inputs.leaf_water_potential,
)
    inputs.environment_instant = environment_instant
    inputs.zenith_angle = zenith_angle
    inputs.direct_horizontal_irradiance = direct_horizontal_irradiance
    inputs.diffuse_horizontal_irradiance = diffuse_horizontal_irradiance
    inputs.ground_reflectance = ground_reflectance
    inputs.ground_temperature = ground_temperature
    inputs.ground_emissivity = ground_emissivity
    inputs.canopy_source_temperature = canopy_source_temperature
    inputs.rainfall = rainfall
    inputs.leaf_water_potential = leaf_water_potential
    return inputs
end

"""
    canopy_energy_balance!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs)

Hourly Picard loop tying every canopy sub-model together. Shortwave, the
above-canopy wind/temperature profile, and rain interception are computed
once per call (none depends on leaf temperature). Longwave, the in-canopy
air profile, and per-layer leaf temperature are then iterated until
`model.convergence` is satisfied, reusing
`AbstractSoilTemperatureConvergence`'s dispatch.

A wetted leaf surface ([`wet_canopy_fraction`](@ref)) blends stomatal
conductance toward [`WET_SURFACE_CONDUCTANCE`](@ref)
([`blend_stomatal_conductance`](@ref)); storage depletion from that
evaporation is applied once, after convergence.

`inputs.ground_temperature`/`canopy_source_temperature` are lagged
(previous-hour) boundary conditions; only canopy-internal state is
Picard-iterated.

Relative humidity is uniform through the canopy for now (the
reference/above-canopy value).

Returns `(; ground_absorbed_shortwave, canopy_absorbed_shortwave,
ground_absorbed_longwave, canopy_absorbed_longwave, ground_throughfall, iterations)`.
"""
function canopy_energy_balance!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs)
    (; site, environment_instant, zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance,
       ground_reflectance, ground_temperature, ground_emissivity, canopy_source_temperature, rainfall,
       leaf_water_potential, vapour_pressure_equation) = inputs

    (; leaf_body, leaf_temperature_prev, sensible_heat_source, evaporation_mass_flow) = buffers.leaf
    leaf_temperature_buffer = buffers.leaf.leaf_temperature  # avoids shadowing the leaf_temperature(solver, ...) function
    (; absorbed_shortwave) = buffers.shortwave
    (; absorbed_longwave, layer_transmission) = buffers.longwave
    (; wind_speed) = buffers.wind
    n_layers = length(leaf_temperature_buffer)
    (; atmospheric_pressure) = environment_instant
    relative_humidity = environment_instant.reference_humidity
    leaf_emissivity = model.leaf_parameters.leaf_emissivity

    shortwave_result = canopy_shortwave!(buffers, model;
        zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance, ground_reflectance)

    wind_result = canopy_wind_profile!(buffers, model, boundary_layer_model;
        site, environment_instant, canopy_source_temperature, vapour_pressure_equation)

    interception_result = canopy_interception!(buffers, model; rainfall, wind_speed)

    fill!(buffers.air_profile.air_temperature, wind_result.canopy_top_air_temperature)
    # LinearizedLeafTemperature divides by leaf_temperature_guess, so it must not start at 0 K.
    fill!(leaf_temperature_buffer, wind_result.canopy_top_air_temperature)

    niter = max_iterations(model.convergence)
    local longwave_result
    iter = 0
    for i in 1:niter
        iter = i
        leaf_temperature_prev .= leaf_temperature_buffer

        longwave_result = canopy_longwave!(buffers, model;
            leaf_temperature=leaf_temperature_buffer, ground_temperature, ground_emissivity,
            site, environment_instant, vapour_pressure_equation)

        @inbounds for layer in 1:n_layers
            # Exchange area is the layer's interception fraction (1-τ), not PAI
            # directly (self-shading means (1-τ)<PAI); canopy_longwave!'s own
            # emission scales the same way.
            leaf_area = 2.0 * (1.0 - layer_transmission[layer]) * 1.0u"m^2"
            # absorbed_shortwave/absorbed_longwave are per unit ground area;
            # leaf_heat_balance/leaf_temperature want per unit leaf_area.
            absorbed_radiation = (absorbed_shortwave[layer] + absorbed_longwave[layer]) / (2.0 * (1.0 - layer_transmission[layer]))
            dry_conductance = stomatal_conductance(model.stomatal_model, zenith_angle, leaf_water_potential)
            conductance = blend_stomatal_conductance(dry_conductance, wet_canopy_fraction(model, buffers, layer))
            air_temperature_layer = buffers.air_profile.air_temperature[layer]

            leaf_temperature_buffer[layer] = leaf_temperature(model.leaf_temperature_solver, absorbed_radiation,
                air_temperature_layer, relative_humidity, wind_speed[layer], atmospheric_pressure,
                leaf_emissivity, conductance, leaf_water_potential, leaf_body;
                leaf_area, leaf_temperature_guess=leaf_temperature_buffer[layer])

            balance = leaf_heat_balance(leaf_temperature_buffer[layer], absorbed_radiation, air_temperature_layer,
                relative_humidity, wind_speed[layer], atmospheric_pressure, leaf_emissivity,
                conductance, leaf_water_potential, leaf_body, leaf_area)
            sensible_heat_source[layer] = balance.conv.convection_flow / 1.0u"m^2"
            evaporation_mass_flow[layer] = balance.evap.transpiration_mass_flow
        end

        canopy_air_profile!(buffers, model, boundary_layer_model;
            canopy_height=model.canopy_height, displacement_height=buffers.wind.displacement_height,
            friction_velocity=wind_result.friction_velocity, canopy_top_air_temperature=wind_result.canopy_top_air_temperature,
            ground_temperature, sensible_heat_source)

        is_converged(model.convergence, iter, niter, leaf_temperature_buffer, leaf_temperature_prev) && break
    end

    @inbounds for layer in 1:n_layers
        evaporated_mass = uconvert(u"kg", evaporation_mass_flow[layer] * CANOPY_TIMESTEP) / 1.0u"m^2"
        deplete_canopy_water!(model, buffers, layer, evaporated_mass)
    end

    return (; shortwave_result..., longwave_result..., interception_result..., iterations=iter)
end
