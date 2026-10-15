"""
    PoolCapacityRainfall()

Default rainfall entry mode. Rain collects in the surface pool, which wets
the top soil node up to saturation each sub-step.
"""
struct PoolCapacityRainfall <: AbstractRainfallEntryMode end

rainfall_flux_for_step(::PoolCapacityRainfall, pool, moisture_timestep, niter_moist) = 0.0u"kg/m^2/s"

rainfall_entry!(soil_moisture, ::PoolCapacityRainfall, pool, step::RainfallEntryStep) =
    _wet_surface_node!(soil_moisture, pool, step.saturation, step.half_thickness, step.frozen_water_content[1])

post_infiltration_pool_update(::PoolCapacityRainfall, pool, step::RainfallEntryStep, water_flux, surface_evaporation) =
    clamp(pool - water_flux - surface_evaporation, 0.0u"kg/m^2", step.max_surface_pool)
