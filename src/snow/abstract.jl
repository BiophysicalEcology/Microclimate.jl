abstract type AbstractSnowModel end

"""
    n_snow_nodes(snow_model) -> Int

Number of snow nodes the model contributes to the ODE state.
"""
function n_snow_nodes end

"""
    init_snow_temperatures(snow_model, T_init)

Initial snow temperature container — `SVector{N}` for an `N`-node variant,
`nothing` for no-snow. Dispatch keeps the result concretely typed.
"""
function init_snow_temperatures end

"""
    initial_snow_state(snow_model, [depth, density])

Construct the initial `SnowState` (or `nothing` for no-snow).
"""
function initial_snow_state end

"""
    allocate_snow_scratch(snow_model, nsteps, num_ode_nodes, depths)

Allocate per-call scratch buffers for the snow hot path.
"""
function allocate_snow_scratch end

"""
    reset_snow_scratch!(snow_model, scratch)

Zero scratch buffers between independent days.
"""
function reset_snow_scratch! end

"""
    activate_snow_nodes!(snow_model, state, scratch, snow_temperature, step)

Mark which snow nodes are currently active given the current depth.
"""
function activate_snow_nodes! end

"""
    compute_effective_depths!(snow_model, state, scratch, soil_depths)

Combine snow node depths with soil depths into the ODE depth vector.
"""
function compute_effective_depths! end

"""
    snow_surface_overrides(snow_model, state, scratch, step)

Surface property overrides (albedo, emissivity, wetness, shade) when snow
covers the ground.
"""
function snow_surface_overrides end

"""
    snow_phase_transition(snow_model, state, scratch, T_snow, T_snow_before, atmospheric_pressure, T_soil_surface)

Phase-change accounting within the snowpack; returns updated state and
clamped temperatures.
"""
function snow_phase_transition end

"""
    update_snow!(snow_model, state, scratch, ...; rainfall_schedule, shade, day_of_year)

Mass-balance update: snowfall, sublimation, rain melt, thermal melt.
"""
function update_snow! end

"""
    adjust_snow_near_nodes!(snow_model, state, scratch, snow_temperature, step)

Tidy snow depth near node-threshold boundaries to avoid thrash.
"""
function adjust_snow_near_nodes! end

# Generic fallback applicable to every snow model variant: phase-change
# accounting in the soil column itself (not the snowpack), delegating to the
# configured freezing scheme. Variants that need bespoke behaviour can override.
function apply_phase_transition(::AbstractSnowModel, freezing_scheme::AbstractSoilFreezingScheme, soil_temperature, soil_temperature_past, buffers, accumulated_latent_heat, soil_moisture, depths)
    phase_transition!(freezing_scheme, buffers;
        temperatures=soil_temperature, temperatures_past=soil_temperature_past,
        accumulated_latent_heat, soil_moisture, depths,
    )
end
