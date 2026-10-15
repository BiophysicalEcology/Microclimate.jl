"""
    ImplicitFluxRainfall()

Rain enters as a flux boundary condition in `infiltration_step!`'s implicit
solve, spread evenly over the hour's sub-steps, so it can wet several layers
in one step.
"""
struct ImplicitFluxRainfall <: AbstractRainfallEntryMode end

rainfall_flux_for_step(::ImplicitFluxRainfall, pool, moisture_timestep, niter_moist) =
    pool / (niter_moist * moisture_timestep)

rainfall_entry!(soil_moisture, ::ImplicitFluxRainfall, pool; kw...) = pool

# water_flux and surface_evaporation already result from the rain applied as
# flux this sub-step, so subtracting them as well would double-count it.
post_infiltration_pool_update(::ImplicitFluxRainfall, pool, water_flux, surface_evaporation;
    rainfall_flux, moisture_timestep,
) = max(0.0u"kg/m^2", pool - rainfall_flux * moisture_timestep)
