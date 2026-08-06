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

Returns `(; ground_absorbed_shortwave, canopy_absorbed_shortwave,
ground_absorbed_longwave, canopy_absorbed_longwave, ground_throughfall,
canopy_potential_transpiration, ground_heat_conductance, ground_vapor_conductance,
iterations)`. `canopy_potential_transpiration` is the canopy-summed unstressed
(water_potential=0) transpiration rate, for feeding a soil-hydraulics demand
term (e.g. `CampbellSoilHydraulics`). `ground_heat_conductance`/
`ground_vapor_conductance` are `model.air_profile_model`'s own ground-to-
lowest-layer conductances, for the soil model's surface exchange to reuse
instead of inverting the log law at a reference height too close to the
roughness length (see `ground_convection_conditions`).
"""
function canopy_energy_balance!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs)
    (; site, environment_instant, zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance,
       ground_reflectance, ground_temperature, ground_emissivity, ground_relative_humidity, canopy_source_temperature,
       rainfall, leaf_water_potential, vapour_pressure_equation) = inputs

    leaf_temperature_buffer = buffers.leaf.leaf_temperature  # avoids shadowing the leaf_temperature(solver, ...) function
    vapour_density_buffer = buffers.air_profile.vapour_density
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
    # vapour_density has no analogous "already warm" check -- it's read as
    # NonlinearSolveCanopyConvergence's starting guess every call, so leaving
    # it at its zero-allocation default on the first call is a degenerate
    # (zero-humidity) corner, not just a poor guess.
    iszero(vapour_density_buffer[1]) && fill!(vapour_density_buffer,
        wet_air_properties(wind_result.canopy_top_air_temperature, wind_result.canopy_top_relative_humidity,
            environment_instant.atmospheric_pressure; vapour_pressure_equation).vapour_density)

    ctx = (; site, environment_instant, zenith_angle, ground_temperature, ground_emissivity, ground_relative_humidity,
        vapour_pressure_equation, leaf_water_potential,
        canopy_top_air_temperature=wind_result.canopy_top_air_temperature,
        canopy_top_relative_humidity=wind_result.canopy_top_relative_humidity,
        friction_velocity=wind_result.friction_velocity, obukhov_length=wind_result.obukhov_length,
    )
    longwave_result, iter = converge_canopy!(buffers, model, boundary_layer_model, ctx, model.convergence_model)

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
        ground_heat_conductance, ground_vapor_conductance, iterations=iter)
end

"""
    _canopy_picard_pass!(buffers, model, boundary_layer_model, ctx, relaxation)

One relaxed Picard pass: longwave exchange, per-layer leaf temperature
(under-relaxed against `buffers.leaf.leaf_temperature_prev`, which the
caller must set beforehand), then the in-canopy air/vapor profile (also
under-relaxed). Shared by [`PicardCanopyConvergence`](@ref)'s own loop and
[`NonlinearSolveCanopyConvergence`](@ref)'s Picard warm-up fallback tier.
Returns `(longwave_result, any_clamped)`.
"""
function _canopy_picard_pass!(buffers, model::MultilayerCanopy, boundary_layer_model, ctx, relaxation)
    (; site, environment_instant, zenith_angle, ground_temperature, ground_emissivity, ground_relative_humidity,
       vapour_pressure_equation, leaf_water_potential, canopy_top_air_temperature, canopy_top_relative_humidity,
       friction_velocity, obukhov_length) = ctx

    (; leaf_body, leaf_temperature_prev, sensible_heat_source, evaporation_mass_flow,
       potential_evaporation_mass_flow, net_balance, air_temperature_prev, relative_humidity_prev) = buffers.leaf
    leaf_temperature_buffer = buffers.leaf.leaf_temperature
    absorbed_radiation_buffer = buffers.leaf.absorbed_radiation
    (; absorbed_shortwave) = buffers.shortwave
    (; absorbed_longwave, layer_transmission) = buffers.longwave
    (; wind_speed, displacement_height) = buffers.wind
    n_layers = length(leaf_temperature_buffer)
    (; atmospheric_pressure) = environment_instant
    leaf_emissivity = model.leaf_parameters.leaf_emissivity
    canopy_height = model.canopy_height

    # A binding clamp can hold the value unchanged between iterations, falsely
    # satisfying the convergence check while net stays nonzero.
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
            leaf_area, leaf_temperature_guess=leaf_temperature_buffer[layer],
            convection_model=model.leaf_convection_model)
        # Under-relax, then clamp to a physically sane range (same bracket
        # RootFindLeafTemperature uses) as a hard backstop. Free/tunable,
        # uncited; asymmetric since radiative/convective heating under low
        # wind can push a leaf further above air temperature than
        # evaporative cooling pushes it below.
        relaxed_leaf_temperature = relaxation * new_leaf_temperature +
            (1.0 - relaxation) * leaf_temperature_prev[layer]
        clamp_lo = air_temperature_layer - 30.0u"K"
        clamp_hi = air_temperature_layer + 40.0u"K"
        any_clamped |= relaxed_leaf_temperature < clamp_lo || relaxed_leaf_temperature > clamp_hi
        leaf_temperature_buffer[layer] = clamp(relaxed_leaf_temperature, clamp_lo, clamp_hi)

        balance = leaf_heat_balance(leaf_temperature_buffer[layer], absorbed_radiation, air_temperature_layer,
            relative_humidity_layer, wind_speed[layer], atmospheric_pressure, leaf_emissivity,
            conductance, leaf_water_potential, leaf_body, leaf_area, model.leaf_convection_model)
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
        canopy_height, displacement_height, friction_velocity, canopy_top_air_temperature,
        canopy_top_relative_humidity, ground_temperature, ground_relative_humidity, sensible_heat_source,
        evaporation_mass_flow, obukhov_length, atmospheric_pressure, vapour_pressure_equation)

    _relax_air_profile!(buffers, model.air_profile_model, relaxation, air_temperature_prev, relative_humidity_prev)

    return longwave_result, any_clamped
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

"""
    converge_canopy!(buffers, model::MultilayerCanopy, boundary_layer_model, ctx, convergence_model)

Drives the per-hour leaf/air-temperature fixed point to convergence, per
`convergence_model` ([`PicardCanopyConvergence`](@ref)/
[`NonlinearSolveCanopyConvergence`](@ref)). `ctx` bundles the per-hour
boundary conditions `canopy_energy_balance!` doesn't already expose via
`buffers`/`model` (site/environment/ground state, `canopy_wind_profile!`'s
canopy-top values). Returns `(longwave_result, iterations)`.
"""
function converge_canopy!(buffers, model::MultilayerCanopy, boundary_layer_model, ctx, cm::PicardCanopyConvergence)
    leaf_temperature_buffer = buffers.leaf.leaf_temperature
    leaf_temperature_prev = buffers.leaf.leaf_temperature_prev
    niter = max_iterations(cm.convergence)
    local longwave_result
    iter = 0
    for i in 1:niter
        iter = i
        leaf_temperature_prev .= leaf_temperature_buffer
        longwave_result, any_clamped = _canopy_picard_pass!(buffers, model, boundary_layer_model, ctx, cm.relaxation)
        !any_clamped && is_converged(cm.convergence, iter, niter, leaf_temperature_buffer, leaf_temperature_prev) && break
    end
    return longwave_result, iter
end

const _CANOPY_TEMPERATURE_SCALE = 10.0        # K
const _CANOPY_VAPOUR_DENSITY_SCALE = 0.03     # kg/m^3
const _CANOPY_LEAF_RESIDUAL_SCALE = 100.0     # W
const _CANOPY_AIR_RESIDUAL_SCALE = 100.0      # W/m^2
const _CANOPY_VAPOUR_RESIDUAL_SCALE = 1.0e-4  # kg/m^2/s
# Sentinel returned when a trial state makes the underlying physics throw
# (e.g. dry_air_properties outside its valid range) -- large enough that
# trust-region/line-search globalization rejects the step, without silently
# clamping the trial state and reporting a false residual for it.
const _CANOPY_RESIDUAL_DOMAIN_FAILURE = 1.0e6

"""
    _canopy_direct_conductances(buffers, model, boundary_layer_model, ctx)

Heat/vapor K-theory conductances for the reduced (layer-1-Dirichlet) system,
computed once per hour -- they depend only on wind/diffusivity, fixed across
every residual evaluation within an hour. Returns a NamedTuple consumed by
[`_canopy_direct_residual`](@ref).
"""
function _canopy_direct_conductances(buffers, model, boundary_layer_model, ctx)
    (; environment_instant, canopy_top_air_temperature, canopy_top_relative_humidity,
       ground_temperature, ground_relative_humidity, friction_velocity, vapour_pressure_equation) = ctx
    (; layer_spacing, relative_eddy_diffusivity, dl, du, dl_v, du_v) = buffers.air_profile
    (; atmospheric_pressure) = environment_instant
    ρ_cp = calc_ρ_cp(canopy_top_air_temperature)
    eddy_diffusivity_top = boundary_layer_model.karman_constant * friction_velocity *
        max(model.canopy_height - buffers.wind.displacement_height, 1.0e-3u"m")
    g_top, g_ground = _ktheory_conductances!(dl, du, relative_eddy_diffusivity, eddy_diffusivity_top,
        layer_spacing, ρ_cp, u"W/m^2/K")
    g_top_v, g_ground_v = _ktheory_conductances!(dl_v, du_v, relative_eddy_diffusivity, eddy_diffusivity_top,
        layer_spacing, 1.0, u"m/s")
    ρv_top = wet_air_properties(canopy_top_air_temperature, canopy_top_relative_humidity,
        atmospheric_pressure; vapour_pressure_equation).vapour_density
    ρv_ground = wet_air_properties(ground_temperature, ground_relative_humidity,
        atmospheric_pressure; vapour_pressure_equation).vapour_density
    return (; g_top, g_ground, dl=copy(dl), du=copy(du), g_top_v, g_ground_v, dl_v=copy(dl_v), du_v=copy(du_v),
        T_top=ustrip(u"K", canopy_top_air_temperature), T_ground=ustrip(u"K", ground_temperature),
        ρv_top=ustrip(u"kg/m^3", ρv_top), ρv_ground=ustrip(u"kg/m^3", ρv_ground))
end

"""
    _canopy_direct_residual(y, p)

Direct (non-fixed-point) canopy residual behind
[`NonlinearSolveCanopyConvergence`](@ref): `y` is the scaled state
`[leaf_temperature(n); air_temperature(2:n); vapour_density(2:n)]` (layer 1's
air/vapor are a Dirichlet boundary, see [`KTheoryAirProfile`](@ref), not
solved). Each block evaluates the actual physics equation directly -- the
leaf block is `leaf_heat_balance(...).net` (no nested root-find), the
air/vapor blocks are the K-theory finite-volume flux balance
`g_prev*(T_prev-T_i) - g_next*(T_i-T_next) + S_i` (no `LinearSolve` cache) --
so no hidden iterative or mutable-state dependency contaminates the Jacobian.
Scaled by [`_CANOPY_TEMPERATURE_SCALE`](@ref)/[`_CANOPY_VAPOUR_DENSITY_SCALE`](@ref)
(state) and [`_CANOPY_LEAF_RESIDUAL_SCALE`](@ref)/[`_CANOPY_AIR_RESIDUAL_SCALE`](@ref)/
[`_CANOPY_VAPOUR_RESIDUAL_SCALE`](@ref) (residual) -- fixed, physically
motivated scales rather than a state-dependent preconditioner, matching the
uniform residual magnitude confirmed empirically once the canopy-top
Dirichlet fix removed the boundary row's ~190x conductance outlier. No hard
clamp: an invalid trial point returns [`_CANOPY_RESIDUAL_DOMAIN_FAILURE`](@ref)
so the solver's own globalization rejects it, instead of silently
substituting a different (clamped) state and reporting its residual.
`p = (; buffers, model, ctx, conductances, n, nfree, call_count)`.
"""
function _canopy_direct_residual(y::AbstractVector, p)
    (; buffers, model, ctx, conductances, n, nfree, call_count) = p
    call_count[] += 1
    (; site, environment_instant, zenith_angle, ground_emissivity, ground_temperature,
       vapour_pressure_equation, leaf_water_potential) = ctx
    (; absorbed_shortwave) = buffers.shortwave
    (; layer_transmission) = buffers.longwave
    (; wind_speed) = buffers.wind
    (; atmospheric_pressure) = environment_instant
    leaf_emissivity = model.leaf_parameters.leaf_emissivity
    leaf_body = buffers.leaf.leaf_body
    leaf_temperature_buffer = buffers.leaf.leaf_temperature
    air_temperature = buffers.air_profile.air_temperature
    vapour_density = buffers.air_profile.vapour_density
    (; g_top, g_ground, dl, g_top_v, g_ground_v, dl_v, T_top, T_ground, ρv_top, ρv_ground) = conductances

    @inbounds for i in 1:n
        leaf_temperature_buffer[i] = clamp(y[i] * _CANOPY_TEMPERATURE_SCALE, 173.15, 373.15) * u"K"
    end
    air_temperature[1] = T_top * u"K"
    vapour_density[1] = ρv_top * u"kg/m^3"
    @inbounds for i in 1:nfree
        air_temperature[i + 1] = clamp(y[n + i] * _CANOPY_TEMPERATURE_SCALE, 173.15, 373.15) * u"K"
        vapour_density[i + 1] = clamp(y[n + nfree + i] * _CANOPY_VAPOUR_DENSITY_SCALE, 0.0, 0.1) * u"kg/m^3"
    end

    local longwave_result
    try
        longwave_result = canopy_longwave!(buffers, model; leaf_temperature=leaf_temperature_buffer,
            ground_temperature, ground_emissivity, site, environment_instant, vapour_pressure_equation)
    catch e
        e isa DomainError || rethrow()
        return fill!(similar(y), _CANOPY_RESIDUAL_DOMAIN_FAILURE), nothing
    end
    absorbed_longwave = buffers.longwave.absorbed_longwave

    R = similar(y)
    sensible = buffers.leaf.sensible_heat_source
    evap = buffers.leaf.evaporation_mass_flow
    try
        @inbounds for i in 1:n
            exchange_fraction = 1.0 - layer_transmission[i]
            air_temperature_i = air_temperature[i]
            rho_v_sat = wet_air_properties(air_temperature_i, 1.0, atmospheric_pressure; vapour_pressure_equation).vapour_density
            relative_humidity_i = clamp(vapour_density[i] / rho_v_sat, 0.0, 1.0)
            if iszero(exchange_fraction)
                R[i] = 0.0
                sensible[i] = zero(sensible[i])
                evap[i] = zero(evap[i])
                continue
            end
            leaf_area = 2.0 * exchange_fraction * 1.0u"m^2"
            absorbed_radiation = (absorbed_shortwave[i] + absorbed_longwave[i]) / (2.0 * exchange_fraction)
            dry_conductance = stomatal_conductance(model.stomatal_model, zenith_angle, leaf_water_potential)
            conductance = blend_stomatal_conductance(dry_conductance, wet_canopy_fraction(model, buffers, i))
            balance = leaf_heat_balance(leaf_temperature_buffer[i], absorbed_radiation, air_temperature_i,
                relative_humidity_i, wind_speed[i], atmospheric_pressure, leaf_emissivity,
                conductance, leaf_water_potential, leaf_body, leaf_area, model.leaf_convection_model)
            R[i] = ustrip(u"W", balance.net) / _CANOPY_LEAF_RESIDUAL_SCALE
            sensible[i] = balance.conv.convection_flow / 1.0u"m^2"
            evap[i] = balance.evap.transpiration_mass_flow
        end
    catch e
        e isa DomainError || rethrow()
        return fill!(similar(y), _CANOPY_RESIDUAL_DOMAIN_FAILURE), nothing
    end

    @inbounds for i in 1:nfree
        g_prev = i == 1 ? g_top : dl[i - 1]
        g_next = i == nfree ? g_ground : dl[i]
        T_prev = i == 1 ? T_top : ustrip(u"K", air_temperature[i])
        T_next = i == nfree ? T_ground : ustrip(u"K", air_temperature[i + 2])
        sensible_i = ustrip(u"W/m^2", sensible[i + 1])
        R[n + i] = (g_prev * (T_prev - ustrip(u"K", air_temperature[i + 1])) -
                    g_next * (ustrip(u"K", air_temperature[i + 1]) - T_next) + sensible_i) / _CANOPY_AIR_RESIDUAL_SCALE

        gv_prev = i == 1 ? g_top_v : dl_v[i - 1]
        gv_next = i == nfree ? g_ground_v : dl_v[i]
        ρv_prev = i == 1 ? ρv_top : ustrip(u"kg/m^3", vapour_density[i])
        ρv_next = i == nfree ? ρv_ground : ustrip(u"kg/m^3", vapour_density[i + 2])
        evap_i = ustrip(u"kg/m^2/s", evap[i + 1] / 1.0u"m^2")
        R[n + nfree + i] = (gv_prev * (ρv_prev - ustrip(u"kg/m^3", vapour_density[i + 1])) -
                            gv_next * (ustrip(u"kg/m^3", vapour_density[i + 1]) - ρv_next) + evap_i) / _CANOPY_VAPOUR_RESIDUAL_SCALE
    end

    return R, longwave_result
end

_canopy_direct_residual_only(y, p) = first(_canopy_direct_residual(y, p))

function converge_canopy!(buffers, model::MultilayerCanopy, boundary_layer_model, ctx, cm::NonlinearSolveCanopyConvergence)
    model.air_profile_model isa KTheoryAirProfile || throw(ArgumentError(
        "NonlinearSolveCanopyConvergence only supports KTheoryAirProfile: its direct residual is built " *
        "from _ktheory_conductances!'s tridiagonal network, not $(typeof(model.air_profile_model))'s " *
        "buffer layout. Use PicardCanopyConvergence with this air_profile_model instead."))
    leaf_temperature_buffer = buffers.leaf.leaf_temperature
    leaf_temperature_prev = buffers.leaf.leaf_temperature_prev
    air_temperature = buffers.air_profile.air_temperature
    vapour_density = buffers.air_profile.vapour_density
    n = length(leaf_temperature_buffer)
    nfree = n - 1

    conductances = _canopy_direct_conductances(buffers, model, boundary_layer_model, ctx)
    call_count = Ref(0)
    p = (; buffers, model, ctx, conductances, n, nfree, call_count)

    function current_y()
        y = Vector{Float64}(undef, n + 2nfree)
        @inbounds for i in 1:n
            y[i] = ustrip(u"K", leaf_temperature_buffer[i]) / _CANOPY_TEMPERATURE_SCALE
        end
        @inbounds for i in 1:nfree
            y[n + i] = ustrip(u"K", air_temperature[i + 1]) / _CANOPY_TEMPERATURE_SCALE
            y[n + nfree + i] = ustrip(u"kg/m^3", vapour_density[i + 1]) / _CANOPY_VAPOUR_DENSITY_SCALE
        end
        return y
    end

    function try_solve(y0)
        prob = SciMLBase.NonlinearProblem{false}(_canopy_direct_residual_only, y0, p)
        sol = try
            SciMLBase.solve(prob, cm.algorithm; abstol=cm.tolerance, maxiters=cm.max_iterations)
        catch e
            # SingularException: the trust-region dogleg step's internal
            # linear solve can hit an exactly-singular local Jacobian
            # (observed in production, not just synthetic tests).
            e isa DomainError || e isa SingularException || rethrow()
            nothing
        end
        (sol === nothing || !SciMLBase.successful_retcode(sol) || !all(isfinite, sol.u)) && return false, y0

        # Don't trust the solver's own retcode alone: verify the residual
        # independently (small and improved over y0), and verify the
        # accepted point stays near y0 (right branch, not a different valid
        # root). `y` is already scaled, so y0 itself is a valid distance
        # reference with no separate state-scale vector needed.
        r0 = _canopy_direct_residual_only(y0, p)
        r1 = _canopy_direct_residual_only(sol.u, p)
        residual_ok = all(isfinite, r1) && maximum(abs, r1) < cm.tolerance && maximum(abs, r1) <= maximum(abs, r0)
        branch_ok = maximum(abs, sol.u .- y0) <= cm.max_branch_distance
        accepted = residual_ok && branch_ok
        return accepted, accepted ? sol.u : y0
    end

    y0 = current_y()
    accepted, y_final = try_solve(y0)

    if !accepted
        # Restore buffers to y0 first: a failed try_solve leaves them at
        # whatever trial point the solver last evaluated internally, not y0
        # -- the warm-up must start from a known state, not that leftover.
        _canopy_direct_residual_only(y0, p)
        # Warm-up tier: a few cheap damped Picard passes to move the state
        # into a better basin (mutates buffers in place) before retrying.
        for _ in 1:cm.picard_warmup_iterations
            leaf_temperature_prev .= leaf_temperature_buffer
            _canopy_picard_pass!(buffers, model, boundary_layer_model, ctx, cm.picard_warmup_relaxation)
        end
        y0_warm = current_y()
        accepted, y_final = try_solve(y0_warm)
        # Even if the retry also fails, the Picard-warmed state is a better
        # accepted output than the original cold guess.
        accepted || (y_final = y0_warm)
    end

    # Re-evaluate at the accepted point: the solver's last internal call
    # isn't necessarily the accepted step (a rejected trust-region/line-
    # search trial), so buffers (sources, longwave_result) can't be trusted
    # without this, same convention as the old fixed-point-map path.
    _residual_at_final, longwave_result = _canopy_direct_residual(y_final, p)
    longwave_result === nothing && ((_residual_at_final, longwave_result) = _canopy_direct_residual(y0, p))

    # Clamp only as a final safety valve on the accepted output, never mid-
    # iteration.
    @inbounds for i in 1:n
        clamp_lo = air_temperature[i] - 30.0u"K"
        clamp_hi = air_temperature[i] + 40.0u"K"
        leaf_temperature_buffer[i] = clamp(y_final[i] * _CANOPY_TEMPERATURE_SCALE * u"K", clamp_lo, clamp_hi)
    end

    return longwave_result, call_count[]
end
