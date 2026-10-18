"""
    PoolCapacityRainfall()

Default rainfall entry mode. Rain collects in the surface pool, which wets
the top soil node up to saturation each sub-step.
"""
struct PoolCapacityRainfall <: AbstractRainfallEntryMode end

rainfall_flux_for_step(::PoolCapacityRainfall, moisture_mode, pool) = 0.0u"kg/m^2/s"

function rainfall_entry!(soil_moisture, ::PoolCapacityRainfall, soil_hydraulic_model, moisture_mode, pool;
    soil_profile, depths, frozen_water_content,
)
    saturation = _saturation_water_content(soil_profile, 1)
    return _wet_surface_node!(soil_moisture, pool, saturation, _layer_thickness(depths, 1),
        conductivity_factor(soil_hydraulic_model.frozen_conductivity, frozen_water_content[1], saturation))
end

post_infiltration_pool_update(::PoolCapacityRainfall, moisture_mode, pool, water_flux, surface_evaporation; rainfall_flux) =
    max(0.0u"kg/m^2", pool - water_flux - surface_evaporation)
