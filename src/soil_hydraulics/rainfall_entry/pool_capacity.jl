"""
    PoolCapacityRainfall()

Default rainfall-entry mode. Rain fills a surface `pool`, which
`_wet_surface_node!` (campbell.jl) then delivers directly into the top soil
node's water content, capped at that node's own storage capacity —
`infiltration_step!` never sees rain as a flux term.
"""
struct PoolCapacityRainfall <: AbstractRainfallEntryMode end

rainfall_flux_for_step(::PoolCapacityRainfall, pool, moisture_timestep, niter_moist) = 0.0u"kg/m^2/s"

apply_rainfall_entry!(::PoolCapacityRainfall, soil_moisture, pool, sat, half_thickness, rainfall_flux, moisture_timestep; kw...) =
    _wet_surface_node!(soil_moisture, pool, sat, half_thickness)

post_infiltration_pool_update(::PoolCapacityRainfall, pool, rainfall_flux, moisture_timestep, water_flux, surf_evap, max_surface_pool) =
    clamp(pool - water_flux - surf_evap, 0.0u"kg/m^2", max_surface_pool)
