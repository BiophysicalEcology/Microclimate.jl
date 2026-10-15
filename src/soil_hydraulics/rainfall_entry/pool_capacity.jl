"""
    PoolCapacityRainfall()

Default rainfall entry mode. Rain collects in the surface pool, which wets
the top soil node up to saturation each sub-step.
"""
struct PoolCapacityRainfall <: AbstractRainfallEntryMode end

rainfall_flux_for_step(::PoolCapacityRainfall, pool, moisture_timestep, niter_moist) = 0.0u"kg/m^2/s"

rainfall_entry!(soil_moisture, ::PoolCapacityRainfall, pool; soil_profile, depths, moisture_timestep) =
    _wet_surface_node!(soil_moisture, pool, _saturation_water_content(soil_profile, 1), _layer_thickness(depths, 1))

post_infiltration_pool_update(::PoolCapacityRainfall, pool, water_flux, surface_evaporation;
    rainfall_flux, moisture_timestep, max_surface_pool,
) = clamp(pool - water_flux - surface_evaporation, 0.0u"kg/m^2", max_surface_pool)
