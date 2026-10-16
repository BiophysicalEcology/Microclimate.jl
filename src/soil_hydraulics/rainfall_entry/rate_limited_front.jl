"""
    RateLimitedFrontRainfall(; fill_fraction=0.9)

Simplified wetting front after the SHAW model's `RAINSL` subroutine. Rain
fills layers from the surface down, each layer taking the smallest of the
remaining rain, its capacity up to `fill_fraction` of saturation, and
`gn * saturated_hydraulic_conductivity * timestep`.
"""
struct RateLimitedFrontRainfall{F} <: AbstractRainfallEntryMode
    fill_fraction::F
end
RateLimitedFrontRainfall(; fill_fraction=0.9) = RateLimitedFrontRainfall(fill_fraction)

rainfall_flux_for_step(::RateLimitedFrontRainfall, moisture_mode, pool) = 0.0u"kg/m^2/s"

function rainfall_entry!(soil_moisture, entry_mode::RateLimitedFrontRainfall, soil_hydraulic_model, moisture_mode, pool; soil_profile, depths)
    pool <= 0.0u"kg/m^2" && return pool
    fill_target = entry_mode.fill_fraction * _saturation_water_content(soil_profile, 1)
    saturated_conductivity = soil_profile.hydraulics.saturated_hydraulic_conductivity
    remaining = pool
    for i in eachindex(soil_moisture)
        remaining <= 0.0u"kg/m^2" && break
        layer_mass = _layer_thickness(depths, i) * water_density
        layer_capacity = max(0.0u"kg/m^2", (fill_target - soil_moisture[i]) * layer_mass)
        rate_limit = uconvert(u"kg/m^2", Unitful.gn * saturated_conductivity[i] * moisture_mode.moisture_timestep)
        infiltrated = min(remaining, layer_capacity, rate_limit)
        soil_moisture[i] += infiltrated / layer_mass
        remaining -= infiltrated
    end
    return remaining
end

post_infiltration_pool_update(::RateLimitedFrontRainfall, moisture_mode, pool, water_flux, surface_evaporation; rainfall_flux) =
    max(0.0u"kg/m^2", pool - water_flux - surface_evaporation)
