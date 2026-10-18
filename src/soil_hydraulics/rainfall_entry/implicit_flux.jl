"""
    ImplicitFluxRainfall()

Rain enters as a flux boundary condition in `infiltration_step!`'s implicit
solve, spread evenly over the hour's sub-steps, so it can wet several layers
in one step.
"""
struct ImplicitFluxRainfall <: AbstractRainfallEntryMode end

rainfall_flux_for_step(::ImplicitFluxRainfall, moisture_mode, pool) =
    pool / (substep_count(moisture_mode) * moisture_mode.moisture_timestep)

rainfall_entry!(soil_moisture, ::ImplicitFluxRainfall, soil_hydraulic_model, moisture_mode, pool; kw...) = pool

# water_flux and surface_evaporation already result from the rain applied as
# flux this sub-step, so subtracting them as well would double-count it.
post_infiltration_pool_update(::ImplicitFluxRainfall, moisture_mode, pool, water_flux, surface_evaporation; rainfall_flux) =
    max(0.0u"kg/m^2", pool - rainfall_flux * moisture_mode.moisture_timestep)
