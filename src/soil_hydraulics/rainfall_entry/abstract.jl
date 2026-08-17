abstract type AbstractRainfallEntryMode end

"""
    rainfall_flux_for_step(mode, pool, moisture_timestep, niter_moist)

Rainfall rate (kg/m²/s, always ≥ 0) to fold into `infiltration_step!`'s own
implicit solve as a flux boundary condition this sub-step. `0` for any mode
that instead delivers water by mutating `soil_moisture` directly in
[`apply_rainfall_entry!`](@ref).
"""
function rainfall_flux_for_step end

"""
    apply_rainfall_entry!(mode, soil_moisture, pool, sat, half_thickness,
                           rainfall_flux, moisture_timestep; depths, soil_profile)

Deliver `pool` into `soil_moisture` before `infiltration_step!` runs, by
whatever mechanism `mode` implements (may mutate `soil_moisture` in place).
Returns the remaining `pool`. A no-op on `soil_moisture` for a mode that
instead enters water via [`rainfall_flux_for_step`](@ref).
"""
function apply_rainfall_entry! end

"""
    post_infiltration_pool_update(mode, pool, rainfall_flux, moisture_timestep,
                                   water_flux, surf_evap, max_surface_pool)

Update `pool` after `infiltration_step!` has run this sub-step.
"""
function post_infiltration_pool_update end
