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

rainfall_flux_for_step(::RateLimitedFrontRainfall, pool, moisture_timestep, niter_moist) = 0.0u"kg/m^2/s"

function rainfall_entry!(soil_moisture, mode::RateLimitedFrontRainfall, pool, step::RainfallEntryStep)
    pool <= 0.0u"kg/m^2" && return pool
    (; depths, half_thickness, timestep, saturation, frozen_conductivity) = step
    fill_target = mode.fill_fraction * saturation
    saturated_conductivity = step.soil_profile.hydraulics.saturated_hydraulic_conductivity
    remaining = pool
    # soil_moisture is total (liquid+ice) content, so the ceiling stays
    # fill_target regardless of ice; ice instead reduces the rate limit.
    for i in eachindex(soil_moisture)
        remaining <= 0.0u"kg/m^2" && break
        ice_content = step.frozen_water_content[i]
        layer_mass = _layer_thickness(depths, half_thickness, i) * water_density
        layer_capacity = max(0.0u"kg/m^2", (fill_target - soil_moisture[i]) * layer_mass)
        conductivity = saturated_conductivity[i] * conductivity_factor(frozen_conductivity, ice_content, saturation)
        rate_limit = uconvert(u"kg/m^2", Unitful.gn * conductivity * timestep)
        infiltrated = min(remaining, layer_capacity, rate_limit)
        soil_moisture[i] += infiltrated / layer_mass
        remaining -= infiltrated
        conductivity_blocked(frozen_conductivity, ice_content, saturation) && break  # frozen layer halts the front
    end
    return remaining
end

post_infiltration_pool_update(::RateLimitedFrontRainfall, pool, step::RainfallEntryStep, water_flux, surface_evaporation) =
    clamp(pool - water_flux - surface_evaporation, 0.0u"kg/m^2", step.max_surface_pool)
