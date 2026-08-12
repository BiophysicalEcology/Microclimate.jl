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
                               vapour_pressure_equation=GoffGratch(), soil_hydraulic_model=nothing)

Per-hour inputs to [`canopy_energy_balance!`](@ref), mirroring
`SoilEnergyInputs`: built once, refreshed in place via
[`update_canopy_energy_balance_inputs!`](@ref). `soil_hydraulic_model`
(`MicroModel.soil_hydraulic_model`) is static like `site` — passed through
to [`stomatal_conductance`](@ref) so models like
[`MoistureResponsiveStomatalConductance`](@ref) share its parameters
instead of storing their own copy.
"""
mutable struct CanopyEnergyBalanceInputs{S,EI,ZA,DHI,DFI,GR,GT,GE,GRH,CST,RF,LWP,VP,SHM}
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
    soil_hydraulic_model::SHM
end

function CanopyEnergyBalanceInputs(model::MultilayerCanopy;
    site, environment_instant, zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance,
    ground_reflectance, ground_temperature, ground_emissivity, ground_relative_humidity, canopy_source_temperature,
    rainfall=0.0u"kg/m^2", leaf_water_potential=model.leaf_parameters.leaf_water_potential,
    vapour_pressure_equation=GoffGratch(), soil_hydraulic_model=nothing,
)
    CanopyEnergyBalanceInputs(site, environment_instant, zenith_angle, direct_horizontal_irradiance,
        diffuse_horizontal_irradiance, ground_reflectance, ground_temperature, ground_emissivity,
        ground_relative_humidity, canopy_source_temperature, rainfall, leaf_water_potential, vapour_pressure_equation,
        soil_hydraulic_model)
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

Hourly solve tying every canopy sub-model together. Shortwave, the
above-canopy wind/temperature profile, and rain interception are computed
once per call (none depends on leaf temperature). Longwave, the in-canopy
air profile, and per-layer leaf temperature are then driven to a fixed point
by [`converge_canopy!`](@ref), dispatching on `model.convergence_model`
([`PicardCanopyConvergence`](@ref)/[`NonlinearSolveCanopyConvergence`](@ref)).

A wetted leaf surface ([`wet_canopy_fraction`](@ref)) blends stomatal
conductance toward [`WET_SURFACE_CONDUCTANCE`](@ref)
([`blend_stomatal_conductance`](@ref)); storage depletion from that
evaporation is applied once, after convergence.

`inputs.ground_temperature`/`canopy_source_temperature` are lagged
(previous-hour) boundary conditions; only canopy-internal state is iterated.

Per-layer relative humidity is resolved by `model.air_profile_model` from
each layer's own evaporation source, the same way per-layer air temperature
is.

Returns `(; ground_absorbed_shortwave, net_absorbed_shortwave,
ground_absorbed_longwave, net_absorbed_longwave, ground_throughfall,
canopy_potential_transpiration, ground_heat_conductance, ground_vapor_conductance,
iterations, canopy_top_air_temperature, canopy_top_relative_humidity,
friction_velocity, obukhov_length, displacement_height)`. `canopy_potential_transpiration`
is the canopy-summed unstressed (water_potential=0) transpiration rate, for
feeding a soil-hydraulics demand term (e.g. `CampbellSoilHydraulics`).
`ground_heat_conductance`/`ground_vapor_conductance` are `model.air_profile_model`'s
own ground-to-lowest-layer conductances, for the soil model's surface exchange to
reuse instead of inverting the log law at a reference height too close to the
roughness length (see `ground_convection_conditions`). The last five fields are
the converged above-canopy boundary state `canopy_air_profile!` was driven with —
diagnostic-oriented (e.g. comparing against an independent implementation's own
boundary solve), not consumed elsewhere in this package.
"""
function canopy_energy_balance!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs)
    (; site, environment_instant, zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance,
       ground_reflectance, ground_temperature, ground_emissivity, ground_relative_humidity, canopy_source_temperature,
       rainfall, leaf_water_potential, vapour_pressure_equation, soil_hydraulic_model) = inputs

    leaf_temperature_buffer = buffers.leaf.leaf_temperature  # avoids shadowing the leaf_temperature(solver, ...) function
    vapor_density_buffer = buffers.air_profile.vapor_density
    air_temperature = buffers.air_profile.air_temperature
    evaporation_mass_flow = buffers.leaf.evaporation_mass_flow
    potential_evaporation_mass_flow = buffers.leaf.potential_evaporation_mass_flow
    (; wind_speed) = buffers.wind
    n_layers = length(leaf_temperature_buffer)

    shortwave_result = canopy_shortwave!(buffers, model;
        zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance, ground_reflectance)

    wind_result = canopy_wind_profile!(buffers, model, boundary_layer_model;
        site, environment_instant, canopy_source_temperature, vapour_pressure_equation)

    interception_result = canopy_interception!(buffers, model; rainfall, wind_speed)

    # Warm-start from last hour's converged state; only reset on the very
    # first call (buffer still zero, and the guess must not be 0 K/0 kg/m^3).
    # air_temperature must warm-start too -- NonlinearSolveCanopyConvergence
    # reads it directly as its starting guess, unlike PicardCanopyConvergence
    # which re-derives it every pass regardless of the starting value.
    iszero(air_temperature[1]) && fill!(air_temperature, wind_result.canopy_top_air_temperature)
    fill!(buffers.air_profile.relative_humidity, wind_result.canopy_top_relative_humidity)
    iszero(leaf_temperature_buffer[1]) && fill!(leaf_temperature_buffer, wind_result.canopy_top_air_temperature)
    # vapor_density has no analogous "already warm" check -- it's read as
    # NonlinearSolveCanopyConvergence's starting guess every call, so leaving
    # it at its zero-allocation default on the first call is a degenerate
    # (zero-humidity) corner, not just a poor guess.
    iszero(vapor_density_buffer[1]) && fill!(vapor_density_buffer,
        wet_air_properties(wind_result.canopy_top_air_temperature, wind_result.canopy_top_relative_humidity,
            environment_instant.atmospheric_pressure; vapour_pressure_equation).vapour_density)

    ctx = (; site, environment_instant, zenith_angle, ground_temperature, ground_emissivity, ground_relative_humidity,
        vapour_pressure_equation, leaf_water_potential, soil_hydraulic_model,
        canopy_top_air_temperature=wind_result.canopy_top_air_temperature,
        canopy_top_relative_humidity=wind_result.canopy_top_relative_humidity,
        friction_velocity=wind_result.friction_velocity, obukhov_length=wind_result.obukhov_length,
    )
    longwave_result, iter, canopy_top_air_temperature, canopy_top_relative_humidity, friction_velocity, obukhov_length =
        converge_canopy!(buffers, model, boundary_layer_model, ctx, model.convergence_model)

    @inbounds for layer in 1:n_layers
        # evaporation_mass_flow blends transpiration and wet-surface
        # evaporation; only the wet fraction should draw down intercepted
        # rain — transpiration draws from soil/plant water instead.
        wet_fraction = wet_canopy_fraction(model, buffers, layer)
        evaporated_mass = uconvert(u"kg", evaporation_mass_flow[layer] * wet_fraction * CANOPY_TIMESTEP) / 1.0u"m^2"
        deplete_canopy_water!(model, buffers, layer, evaporated_mass)
    end

    canopy_potential_transpiration = sum(potential_evaporation_mass_flow) / 1.0u"m^2"
    ground_heat_conductance = buffers.air_profile.ground_heat_conductance[]
    ground_vapor_conductance = buffers.air_profile.ground_vapor_conductance[]

    return (; shortwave_result..., longwave_result..., interception_result..., canopy_potential_transpiration,
        ground_heat_conductance, ground_vapor_conductance, iterations=iter,
        canopy_top_air_temperature, canopy_top_relative_humidity,
        friction_velocity, obukhov_length,
        displacement_height=buffers.wind.displacement_height)
end

"""
    _canopy_picard_pass!(buffers, model, boundary_layer_model, ctx, relaxation)

One relaxed Picard pass: re-anchor the above-canopy boundary on this pass's
own top-leaf temperature, longwave exchange, per-layer leaf temperature
(under-relaxed against `buffers.leaf.leaf_temperature_prev`, which the
caller must set beforehand), then the in-canopy air/vapor profile (also
under-relaxed). Called once per [`PicardCanopyConvergence`](@ref) iteration.
Returns `(longwave_result, any_clamped, canopy_top_air_temperature,
canopy_top_relative_humidity, friction_velocity, obukhov_length)` -- the last
two are this pass's own (`leaf_temperature_buffer[1]`-driven) values, the
ones `canopy_air_profile!` is actually driven with below, not `ctx`'s
one-time pre-loop seed.
"""
function _canopy_picard_pass!(buffers, model::MultilayerCanopy, boundary_layer_model, ctx, relaxation,
                               canopy_top_air_temperature_prev, canopy_top_relative_humidity_prev)
    (; site, environment_instant, zenith_angle, ground_temperature, ground_emissivity, ground_relative_humidity,
       vapour_pressure_equation, leaf_water_potential, soil_hydraulic_model) = ctx

    (; leaf_body, leaf_temperature_prev, sensible_heat_source, evaporation_mass_flow,
       potential_evaporation_mass_flow, net_balance, air_temperature_prev, relative_humidity_prev) = buffers.leaf
    leaf_temperature_buffer = buffers.leaf.leaf_temperature
    absorbed_radiation_buffer = buffers.leaf.absorbed_radiation
    (; absorbed_shortwave) = buffers.shortwave
    (; absorbed_longwave, layer_transmission) = buffers.longwave
    displacement_height = buffers.wind.displacement_height
    n_layers = length(leaf_temperature_buffer)
    (; atmospheric_pressure) = environment_instant
    leaf_emissivity = model.leaf_parameters.leaf_emissivity
    canopy_height = model.canopy_height

    # Wind/stability only here -- canopy_source_temperature (lagged top-leaf
    # temperature) sets stability, not the canopy-top boundary temperature
    # itself (see canopy_top_flux_boundary below).
    wind_result = canopy_wind_profile!(buffers, model, boundary_layer_model;
        site, environment_instant, canopy_source_temperature=leaf_temperature_buffer[1], vapour_pressure_equation)
    friction_velocity = wind_result.friction_velocity
    obukhov_length = wind_result.obukhov_length
    wind_speed = buffers.wind.wind_speed

    # A binding clamp can hold the value unchanged between iterations, falsely
    # satisfying the convergence check while net stays nonzero.
    any_clamped = false

    longwave_result = canopy_longwave!(buffers, model;
        leaf_temperature=leaf_temperature_buffer, ground_temperature, ground_emissivity,
        site, environment_instant, vapour_pressure_equation)

    # Raw leaf-temperature solve; Aitken-blended as a whole vector below, so
    # solve and flux/clamp are separate passes either side of the blend.
    @inbounds for layer in 1:n_layers
        # Exchange area is the layer's interception fraction (1-τ), not PAI
        # directly (self-shading means (1-τ)<PAI); canopy_longwave!'s own
        # emission scales the same way. 2.0 = both leaf surfaces.
        exchange_fraction = 1.0 - layer_transmission[layer]
        air_temperature_layer = buffers.air_profile.air_temperature[layer]

        if iszero(exchange_fraction)
            # No leaf area (zero PAI): avoid 0/0 below, no distinct leaf temperature.
            leaf_temperature_buffer[layer] = air_temperature_layer
            absorbed_radiation_buffer[layer] = zero(absorbed_radiation_buffer[layer])
            continue
        end

        relative_humidity_layer = buffers.air_profile.relative_humidity[layer]
        leaf_area = 2.0 * exchange_fraction * 1.0u"m^2"
        # absorbed_shortwave/absorbed_longwave are per unit ground area;
        # leaf_heat_balance/leaf_temperature want per unit leaf_area.
        absorbed_radiation = (absorbed_shortwave[layer] + absorbed_longwave[layer]) / (2.0 * exchange_fraction)
        dry_conductance = stomatal_conductance(model.stomatal_model, zenith_angle, leaf_water_potential, soil_hydraulic_model)
        conductance = blend_stomatal_conductance(dry_conductance, wet_canopy_fraction(model, buffers, layer))

        leaf_temperature_buffer[layer] = leaf_temperature(model.leaf_temperature_solver, absorbed_radiation,
            air_temperature_layer, relative_humidity_layer, wind_speed[layer], atmospheric_pressure,
            leaf_emissivity, conductance, leaf_water_potential, leaf_body;
            leaf_area, leaf_temperature_guess=leaf_temperature_buffer[layer],
            convection_model=model.leaf_convection_model)
        absorbed_radiation_buffer[layer] = absorbed_radiation
    end

    # Uniform across layers (weight_bottom==weight_top): layer heights aren't
    # available for every air_profile_model, so no bottom/top emphasis here.
    _aitken_relax!(leaf_temperature_buffer, leaf_temperature_prev, buffers.leaf.aitken_dummy_heights, 1.0u"m",
        buffers.leaf.aitken_omega, buffers.leaf.aitken_have_prev, buffers.leaf.aitken_residual_prev,
        0.02, 0.9, 1.0, 1.0, 0.0, u"K")

    @inbounds for layer in 1:n_layers
        exchange_fraction = 1.0 - layer_transmission[layer]
        if iszero(exchange_fraction)
            sensible_heat_source[layer] = zero(sensible_heat_source[layer])
            evaporation_mass_flow[layer] = zero(evaporation_mass_flow[layer])
            net_balance[layer] = zero(net_balance[layer])
            potential_evaporation_mass_flow[layer] = zero(potential_evaporation_mass_flow[layer])
            continue
        end

        air_temperature_layer = buffers.air_profile.air_temperature[layer]
        relative_humidity_layer = buffers.air_profile.relative_humidity[layer]
        leaf_area = 2.0 * exchange_fraction * 1.0u"m^2"
        absorbed_radiation = absorbed_radiation_buffer[layer]
        dry_conductance = stomatal_conductance(model.stomatal_model, zenith_angle, leaf_water_potential, soil_hydraulic_model)
        conductance = blend_stomatal_conductance(dry_conductance, wet_canopy_fraction(model, buffers, layer))

        # Clamp to a physically sane range (same bracket RootFindLeafTemperature
        # uses) as a hard backstop against a runaway Aitken step. Free/tunable,
        # uncited; asymmetric since radiative/convective heating under low
        # wind can push a leaf further above air temperature than
        # evaporative cooling pushes it below.
        clamp_lo = air_temperature_layer - 30.0u"K"
        clamp_hi = air_temperature_layer + 40.0u"K"
        any_clamped |= leaf_temperature_buffer[layer] < clamp_lo || leaf_temperature_buffer[layer] > clamp_hi
        leaf_temperature_buffer[layer] = clamp(leaf_temperature_buffer[layer], clamp_lo, clamp_hi)

        balance = leaf_heat_balance(leaf_temperature_buffer[layer], absorbed_radiation, air_temperature_layer,
            relative_humidity_layer, wind_speed[layer], atmospheric_pressure, leaf_emissivity,
            conductance, leaf_water_potential, leaf_body, leaf_area, model.leaf_convection_model)
        sensible_heat_source[layer] = balance.conv.convection_flow / 1.0u"m^2"
        evaporation_mass_flow[layer] = balance.evap.transpiration_mass_flow
        net_balance[layer] = balance.net

        # Unstressed (water_potential=0) reference transpiration, for a
        # soil-hydraulics demand term — not the leaf's own energy balance.
        # Conductance is recomputed at water_potential=0 too, so this stays a
        # true unstressed demand under MoistureResponsiveStomatalConductance.
        unstressed_dry_conductance = stomatal_conductance(model.stomatal_model, zenith_angle, 0.0u"J/kg", soil_hydraulic_model)
        atmos = AtmosphericConditions(relative_humidity_layer, wind_speed[layer], atmospheric_pressure)
        potential_evap = HeatExchange.evaporation(unstressed_dry_conductance, balance.conv.mass_transfer_coefficient,
            atmos, leaf_area, leaf_temperature_buffer[layer], air_temperature_layer; water_potential=0.0u"J/kg")
        potential_evaporation_mass_flow[layer] = potential_evap.transpiration_mass_flow
    end

    # Whole-canopy aggregated-flux boundary (canopy_top_flux_boundary,
    # wind/attenuation.jl) -- not a bare-ground MOST surface with a single
    # leaf's temperature standing in for the whole canopy's own boundary
    # condition. Relaxed across passes same as before.
    total_sensible_flux = sum(sensible_heat_source)
    total_latent_flux = sum(evaporation_mass_flow) / 1.0u"m^2"
    reference_height = buffers.wind.above_canopy_profile.heights[2]
    reference_temperature = environment_instant.reference_temperature
    reference_vapour_density = wet_air_properties(reference_temperature, environment_instant.reference_humidity,
        atmospheric_pressure; vapour_pressure_equation).vapour_density
    ground_vapour_density = wet_air_properties(ground_temperature, ground_relative_humidity,
        atmospheric_pressure; vapour_pressure_equation).vapour_density
    boundary = canopy_top_flux_boundary(boundary_layer_model, canopy_height, displacement_height,
        buffers.wind.roughness_length, reference_height, friction_velocity, obukhov_length,
        total_sensible_flux, total_latent_flux, ground_temperature, ground_vapour_density,
        reference_temperature, reference_vapour_density)
    saturation_vapour_density = wet_air_properties(boundary.top_temperature, 1.0, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    boundary_relative_humidity = clamp(
        ustrip(u"kg/m^3", boundary.top_vapour_density) / ustrip(u"kg/m^3", saturation_vapour_density), 0.0, 1.0)
    canopy_top_air_temperature = relaxation * boundary.top_temperature +
        (1.0 - relaxation) * canopy_top_air_temperature_prev
    # boundary.top_temperature is already bracketed by boundary.bound_lo/hi,
    # but canopy_top_air_temperature is a running blend that carries forward
    # pass to pass (and hour to hour) -- under heavy relaxation an already-
    # extreme carried-forward value barely moves per pass, so it needs the
    # same absolute bracket, not just a bounded input to blend towards.
    canopy_top_air_temperature = clamp(canopy_top_air_temperature, boundary.bound_lo, boundary.bound_hi)
    # Raupach :bulk further tightens this to reference_temperature/ground_temperature
    # (fixed ambient data) ± bulk_temperature_margin. Deliberately NOT widened by
    # leaf_temperature the way the interior profile's own bound is: leaf is clamped
    # relative to air_temperature, so under a sustained calm/cold-air runaway leaf
    # collapses in lockstep with it -- widening this bound via leaf lets the "floor"
    # collapse together with the thing it's meant to catch.
    air_profile_model = model.air_profile_model
    if air_profile_model isa RaupachLTheoryAirProfile && air_profile_model.far_field_mode isa Val{:bulk}
        margin = air_profile_model.bulk_temperature_margin
        tight_lo = min(reference_temperature, ground_temperature) - margin
        tight_hi = max(reference_temperature, ground_temperature) + margin
        canopy_top_air_temperature = clamp(canopy_top_air_temperature, tight_lo, tight_hi)
    end
    canopy_top_relative_humidity = relaxation * boundary_relative_humidity +
        (1.0 - relaxation) * canopy_top_relative_humidity_prev

    air_temperature_prev .= buffers.air_profile.air_temperature
    relative_humidity_prev .= buffers.air_profile.relative_humidity

    canopy_air_profile!(buffers, model, boundary_layer_model;
        canopy_height, displacement_height, friction_velocity, wind_attenuation=buffers.wind.wind_attenuation,
        canopy_top_air_temperature, canopy_top_relative_humidity, ground_temperature, ground_relative_humidity,
        sensible_heat_source, evaporation_mass_flow, obukhov_length, atmospheric_pressure, vapour_pressure_equation,
        leaf_temperature=leaf_temperature_buffer)

    _relax_air_profile!(buffers, model.air_profile_model, relaxation, air_temperature_prev, relative_humidity_prev)

    return longwave_result, any_clamped, canopy_top_air_temperature, canopy_top_relative_humidity, friction_velocity, obukhov_length
end

# Under-relax the air profile: it feeds back into next pass's leaf
# temperature via sensible_heat_source, so damping only the leaf side
# doesn't stop a joint leaf/air runaway. Skipped for air-profile models that
# already stabilize themselves internally -- RaupachLTheoryAirProfile runs
# its own adaptive Aitken acceleration inside canopy_air_profile!, and
# compounding a second fixed relaxation on top of that double-damps it.
function _relax_air_profile!(buffers, ::AbstractCanopyAirProfileModel, relaxation, air_temperature_prev, relative_humidity_prev)
    buffers.air_profile.air_temperature .= relaxation .* buffers.air_profile.air_temperature .+
        (1.0 - relaxation) .* air_temperature_prev
    buffers.air_profile.relative_humidity .= relaxation .* buffers.air_profile.relative_humidity .+
        (1.0 - relaxation) .* relative_humidity_prev
    return nothing
end
_relax_air_profile!(buffers, ::RaupachLTheoryAirProfile, relaxation, air_temperature_prev, relative_humidity_prev) = nothing

# Whether the canopy-top boundary (re-anchored each pass by
# _canopy_picard_pass!) has itself stopped moving -- checked alongside the
# per-layer leaf-temperature convergence so a still-drifting boundary can't
# hide behind an already-converged leaf array. FixedIterationConvergence
# always runs its full count regardless, so this is a no-op for it.
_canopy_top_converged(::FixedIterationConvergence, canopy_top_air_temperature, canopy_top_air_temperature_prev) = true
_canopy_top_converged(c::IterationToleranceConvergence, canopy_top_air_temperature, canopy_top_air_temperature_prev) =
    abs(canopy_top_air_temperature - canopy_top_air_temperature_prev) < c.tolerance

"""
    converge_canopy!(buffers, model::MultilayerCanopy, boundary_layer_model, ctx, convergence_model)

Drives the per-hour leaf/air-temperature fixed point to convergence, per
`convergence_model` ([`PicardCanopyConvergence`](@ref)). `ctx` bundles the
per-hour boundary conditions `canopy_energy_balance!` doesn't already expose
via `buffers`/`model` (site/environment/ground state). Re-anchors the
canopy-top boundary on the current pass's own leaf temperature each pass
(see `_canopy_picard_pass!`). Returns `(longwave_result, iterations,
canopy_top_air_temperature, canopy_top_relative_humidity, friction_velocity,
obukhov_length)` -- the last four are the final pass's own converged values.
"""
function converge_canopy!(buffers, model::MultilayerCanopy, boundary_layer_model, ctx, cm::PicardCanopyConvergence)
    leaf_temperature_buffer = buffers.leaf.leaf_temperature
    leaf_temperature_prev = buffers.leaf.leaf_temperature_prev
    niter = max_iterations(cm.convergence)
    local longwave_result, friction_velocity, obukhov_length
    canopy_top_air_temperature_prev = ctx.canopy_top_air_temperature
    canopy_top_relative_humidity_prev = ctx.canopy_top_relative_humidity
    iter = 0
    for i in 1:niter
        iter = i
        leaf_temperature_prev .= leaf_temperature_buffer
        longwave_result, any_clamped, canopy_top_air_temperature, canopy_top_relative_humidity, friction_velocity, obukhov_length =
            _canopy_picard_pass!(buffers, model, boundary_layer_model, ctx, cm.relaxation,
                                  canopy_top_air_temperature_prev, canopy_top_relative_humidity_prev)
        boundary_converged = _canopy_top_converged(cm.convergence, canopy_top_air_temperature, canopy_top_air_temperature_prev)
        canopy_top_air_temperature_prev = canopy_top_air_temperature
        canopy_top_relative_humidity_prev = canopy_top_relative_humidity
        !any_clamped && boundary_converged &&
            is_converged(cm.convergence, iter, niter, leaf_temperature_buffer, leaf_temperature_prev) && break
    end
    return longwave_result, iter, canopy_top_air_temperature_prev, canopy_top_relative_humidity_prev, friction_velocity, obukhov_length
end
