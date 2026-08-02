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
                               ground_relative_humidity, canopy_source_temperature, rainfall=0.0u"kg/m^2",
                               leaf_water_potential=model.leaf_parameters.leaf_water_potential,
                               vapour_pressure_equation=GoffGratch())

Per-hour inputs to [`canopy_energy_balance!`](@ref), mirroring
`SoilEnergyInputs`: built once, refreshed in place via
[`update_canopy_energy_balance_inputs!`](@ref).
"""
mutable struct CanopyEnergyBalanceInputs{S,EI,ZA,DHI,DFI,GR,GT,GE,GRH,CST,RF,LWP,VP}
    site::S
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
    vapour_pressure_equation::VP
end

function CanopyEnergyBalanceInputs(model::MultilayerCanopy;
    site, environment_instant, zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance,
    ground_reflectance, ground_temperature, ground_emissivity, ground_relative_humidity, canopy_source_temperature,
    rainfall=0.0u"kg/m^2", leaf_water_potential=model.leaf_parameters.leaf_water_potential,
    vapour_pressure_equation=GoffGratch(),
)
    CanopyEnergyBalanceInputs(site, environment_instant, zenith_angle, direct_horizontal_irradiance,
        diffuse_horizontal_irradiance, ground_reflectance, ground_temperature, ground_emissivity,
        ground_relative_humidity, canopy_source_temperature, rainfall, leaf_water_potential, vapour_pressure_equation)
end

"""
    update_canopy_energy_balance_inputs!(inputs::CanopyEnergyBalanceInputs; environment_instant, zenith_angle,
                                          direct_horizontal_irradiance, diffuse_horizontal_irradiance,
                                          ground_reflectance, ground_temperature, ground_emissivity,
                                          ground_relative_humidity, canopy_source_temperature, rainfall=0.0u"kg/m^2",
                                          leaf_water_potential=inputs.leaf_water_potential)

In-place per-hour update, no allocation. `site`/`vapour_pressure_equation`
aren't refreshed (static for a run).
"""
function update_canopy_energy_balance_inputs!(inputs::CanopyEnergyBalanceInputs;
    environment_instant, zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance,
    ground_reflectance, ground_temperature, ground_emissivity, ground_relative_humidity, canopy_source_temperature,
    rainfall=0.0u"kg/m^2", leaf_water_potential=inputs.leaf_water_potential,
)
    inputs.environment_instant = environment_instant
    inputs.zenith_angle = zenith_angle
    inputs.direct_horizontal_irradiance = direct_horizontal_irradiance
    inputs.diffuse_horizontal_irradiance = diffuse_horizontal_irradiance
    inputs.ground_reflectance = ground_reflectance
    inputs.ground_temperature = ground_temperature
    inputs.ground_emissivity = ground_emissivity
    inputs.ground_relative_humidity = ground_relative_humidity
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

Per-layer relative humidity is resolved by `model.air_profile_model` from
each layer's own evaporation source, the same way per-layer air temperature
is.

Returns `(; ground_absorbed_shortwave, canopy_absorbed_shortwave,
ground_absorbed_longwave, canopy_absorbed_longwave, ground_throughfall,
canopy_potential_transpiration, iterations)`. `canopy_potential_transpiration`
is the canopy-summed unstressed (water_potential=0) transpiration rate, for
feeding a soil-hydraulics demand term (e.g. `CampbellSoilHydraulics`).
"""
function canopy_energy_balance!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs)
    (; site, environment_instant, zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance,
       ground_reflectance, ground_temperature, ground_emissivity, ground_relative_humidity, canopy_source_temperature,
       rainfall, leaf_water_potential, vapour_pressure_equation) = inputs

    (; leaf_body, leaf_temperature_prev, sensible_heat_source, evaporation_mass_flow,
       potential_evaporation_mass_flow, net_balance, air_temperature_prev, relative_humidity_prev) = buffers.leaf
    leaf_temperature_buffer = buffers.leaf.leaf_temperature  # avoids shadowing the leaf_temperature(solver, ...) function
    absorbed_radiation_buffer = buffers.leaf.absorbed_radiation  # avoids shadowing the per-layer local below
    (; absorbed_shortwave) = buffers.shortwave
    (; absorbed_longwave, layer_transmission) = buffers.longwave
    (; wind_speed) = buffers.wind
    n_layers = length(leaf_temperature_buffer)
    (; atmospheric_pressure) = environment_instant
    leaf_emissivity = model.leaf_parameters.leaf_emissivity

    shortwave_result = canopy_shortwave!(buffers, model;
        zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance, ground_reflectance)

    wind_result = canopy_wind_profile!(buffers, model, boundary_layer_model;
        site, environment_instant, canopy_source_temperature, vapour_pressure_equation)

    interception_result = canopy_interception!(buffers, model; rainfall, wind_speed)

    fill!(buffers.air_profile.air_temperature, wind_result.canopy_top_air_temperature)
    fill!(buffers.air_profile.relative_humidity, wind_result.canopy_top_relative_humidity)
    # Warm-start from last hour's converged leaf temperatures; only reset on
    # the very first call (buffer still zero, and the guess must not be 0 K).
    iszero(leaf_temperature_buffer[1]) && fill!(leaf_temperature_buffer, wind_result.canopy_top_air_temperature)

    niter = max_iterations(model.convergence)
    local longwave_result
    iter = 0
    for i in 1:niter
        iter = i
        leaf_temperature_prev .= leaf_temperature_buffer
        # A binding clamp can hold the value unchanged between iterations,
        # falsely satisfying the convergence check while net stays nonzero.
        any_clamped = false

        longwave_result = canopy_longwave!(buffers, model;
            leaf_temperature=leaf_temperature_buffer, ground_temperature, ground_emissivity,
            site, environment_instant, vapour_pressure_equation)

        @inbounds for layer in 1:n_layers
            # Exchange area is the layer's interception fraction (1-τ), not PAI
            # directly (self-shading means (1-τ)<PAI); canopy_longwave!'s own
            # emission scales the same way. 2.0 = both leaf surfaces.
            exchange_fraction = 1.0 - layer_transmission[layer]
            air_temperature_layer = buffers.air_profile.air_temperature[layer]
            relative_humidity_layer = buffers.air_profile.relative_humidity[layer]

            if iszero(exchange_fraction)
                # No leaf area (zero PAI): avoid 0/0 below, no distinct leaf temperature.
                leaf_temperature_buffer[layer] = air_temperature_layer
                sensible_heat_source[layer] = zero(sensible_heat_source[layer])
                evaporation_mass_flow[layer] = zero(evaporation_mass_flow[layer])
                absorbed_radiation_buffer[layer] = zero(absorbed_radiation_buffer[layer])
                net_balance[layer] = zero(net_balance[layer])
                potential_evaporation_mass_flow[layer] = zero(potential_evaporation_mass_flow[layer])
                continue
            end

            leaf_area = 2.0 * exchange_fraction * 1.0u"m^2"
            # absorbed_shortwave/absorbed_longwave are per unit ground area;
            # leaf_heat_balance/leaf_temperature want per unit leaf_area.
            absorbed_radiation = (absorbed_shortwave[layer] + absorbed_longwave[layer]) / (2.0 * exchange_fraction)
            dry_conductance = stomatal_conductance(model.stomatal_model, zenith_angle, leaf_water_potential)
            conductance = blend_stomatal_conductance(dry_conductance, wet_canopy_fraction(model, buffers, layer))

            new_leaf_temperature = leaf_temperature(model.leaf_temperature_solver, absorbed_radiation,
                air_temperature_layer, relative_humidity_layer, wind_speed[layer], atmospheric_pressure,
                leaf_emissivity, conductance, leaf_water_potential, leaf_body;
                leaf_area, leaf_temperature_guess=leaf_temperature_buffer[layer])
            # Under-relax, then clamp to a physically sane range (same bracket
            # RootFindLeafTemperature uses) as a hard backstop. Free/tunable,
            # uncited; asymmetric since radiative/convective heating under low
            # wind can push a leaf further above air temperature than
            # evaporative cooling pushes it below.
            relaxed_leaf_temperature = model.relaxation * new_leaf_temperature +
                (1.0 - model.relaxation) * leaf_temperature_prev[layer]
            clamp_lo = air_temperature_layer - 30.0u"K"
            clamp_hi = air_temperature_layer + 40.0u"K"
            any_clamped |= relaxed_leaf_temperature < clamp_lo || relaxed_leaf_temperature > clamp_hi
            leaf_temperature_buffer[layer] = clamp(relaxed_leaf_temperature, clamp_lo, clamp_hi)

            balance = leaf_heat_balance(leaf_temperature_buffer[layer], absorbed_radiation, air_temperature_layer,
                relative_humidity_layer, wind_speed[layer], atmospheric_pressure, leaf_emissivity,
                conductance, leaf_water_potential, leaf_body, leaf_area)
            sensible_heat_source[layer] = balance.conv.convection_flow / 1.0u"m^2"
            evaporation_mass_flow[layer] = balance.evap.transpiration_mass_flow
            absorbed_radiation_buffer[layer] = absorbed_radiation
            net_balance[layer] = balance.net

            # Unstressed (water_potential=0) reference transpiration, for a
            # soil-hydraulics demand term — not the leaf's own energy balance.
            atmos = AtmosphericConditions(relative_humidity_layer, wind_speed[layer], atmospheric_pressure)
            potential_evap = HeatExchange.evaporation(dry_conductance, balance.conv.mass_transfer_coefficient,
                atmos, leaf_area, leaf_temperature_buffer[layer], air_temperature_layer; water_potential=0.0u"J/kg")
            potential_evaporation_mass_flow[layer] = potential_evap.transpiration_mass_flow
        end

        air_temperature_prev .= buffers.air_profile.air_temperature
        relative_humidity_prev .= buffers.air_profile.relative_humidity

        canopy_air_profile!(buffers, model, boundary_layer_model;
            canopy_height=model.canopy_height, displacement_height=buffers.wind.displacement_height,
            friction_velocity=wind_result.friction_velocity, canopy_top_air_temperature=wind_result.canopy_top_air_temperature,
            canopy_top_relative_humidity=wind_result.canopy_top_relative_humidity,
            ground_temperature, ground_relative_humidity, sensible_heat_source, evaporation_mass_flow,
            obukhov_length=wind_result.obukhov_length, atmospheric_pressure, vapour_pressure_equation)

        # Under-relax the air profile too: it feeds back into next iteration's
        # leaf temperature via sensible_heat_source, so damping only the leaf
        # side doesn't stop a joint leaf/air runaway.
        buffers.air_profile.air_temperature .= model.relaxation .* buffers.air_profile.air_temperature .+
            (1.0 - model.relaxation) .* air_temperature_prev
        buffers.air_profile.relative_humidity .= model.relaxation .* buffers.air_profile.relative_humidity .+
            (1.0 - model.relaxation) .* relative_humidity_prev

        !any_clamped && is_converged(model.convergence, iter, niter, leaf_temperature_buffer, leaf_temperature_prev) && break
    end

    @inbounds for layer in 1:n_layers
        # evaporation_mass_flow blends transpiration and wet-surface
        # evaporation; only the wet fraction should draw down intercepted
        # rain — transpiration draws from soil/plant water instead.
        wet_fraction = wet_canopy_fraction(model, buffers, layer)
        evaporated_mass = uconvert(u"kg", evaporation_mass_flow[layer] * wet_fraction * CANOPY_TIMESTEP) / 1.0u"m^2"
        deplete_canopy_water!(model, buffers, layer, evaporated_mass)
    end

    canopy_potential_transpiration = sum(potential_evaporation_mass_flow) / 1.0u"m^2"

    return (; shortwave_result..., longwave_result..., interception_result..., canopy_potential_transpiration, iterations=iter)
end
