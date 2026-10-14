"""
    ImplicitFluxRainfall()

Rain enters as a flux boundary condition inside `infiltration_step!`'s own
implicit solve — the same role `evaporation_potential`/`vapor_flux[1]`
already plays — rather than pre-loading the top soil node from a `pool`.
Node 1 updates purely from the Newton solve, like every other layer, so a
single hour's rain can influence multiple layers within one implicit step
instead of needing several sub-steps to fill the surface node's own
capacity first.
"""
struct ImplicitFluxRainfall <: AbstractRainfallEntryMode end

rainfall_flux_for_step(::ImplicitFluxRainfall, pool, moisture_timestep, niter_moist) =
    pool / (niter_moist * moisture_timestep)

apply_rainfall_entry!(::ImplicitFluxRainfall, soil_moisture, pool, sat, half_thickness, rainfall_flux, moisture_timestep; kw...) = pool

# water_flux/surf_evap are already downstream consequences of the rain that
# entered via rainfall_flux this sub-step -- subtracting them here too would
# double-count it.
post_infiltration_pool_update(::ImplicitFluxRainfall, pool, rainfall_flux, moisture_timestep, water_flux, surf_evap, max_surface_pool) =
    max(0.0u"kg/m^2", pool - rainfall_flux * moisture_timestep)
