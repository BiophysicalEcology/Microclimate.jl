"""
    RateLimitedFrontRainfall(; fill_fraction=0.9)

Stateless, single-pass Green-Ampt-style wetting-front march (Shaw 1990's
`RAINSL`, simplified): rain advances layer by layer from the surface down,
each layer accepting the lesser of its remaining storage capacity (up to
`fill_fraction * sat`) and a conductivity-rate limit (`gn * saturated_hydraulic_conductivity`
— the same gravity-driven-flux relationship `infiltration_step!` already
uses for `drainage`), for as many layers as the timestep's rain and time
allow. Unlike [`PoolCapacityRainfall`](@ref), a strong pulse can wet several
layers in one call; unlike [`ImplicitFluxRainfall`](@ref), it never touches
`infiltration_step!`'s implicit solve. Front position isn't tracked across
calls — each call re-derives how far rain penetrates from the current
`soil_moisture` profile alone.
"""
struct RateLimitedFrontRainfall{F} <: AbstractRainfallEntryMode
    fill_fraction::F
end
RateLimitedFrontRainfall(; fill_fraction=0.9) = RateLimitedFrontRainfall(fill_fraction)

rainfall_flux_for_step(::RateLimitedFrontRainfall, pool, moisture_timestep, niter_moist) = 0.0u"kg/m^2/s"

function apply_rainfall_entry!(mode::RateLimitedFrontRainfall, soil_moisture, pool, sat, half_thickness, rainfall_flux, moisture_timestep; depths, soil_profile)
    pool <= 0.0u"kg/m^2" && return pool
    water_density = 1000.0u"kg/m^3"
    fill_target = mode.fill_fraction * sat
    saturated_conductivity = soil_profile.hydraulics.saturated_hydraulic_conductivity
    n = length(soil_moisture)
    remaining = pool
    @inbounds for i in 1:n
        remaining <= 0.0u"kg/m^2" && break
        layer_thickness = i == 1 ? half_thickness :
            i == n ? (depths[n] - depths[n-1]) / 2 : (depths[i+1] - depths[i-1]) / 2
        layer_capacity = max(0.0u"kg/m^2", (fill_target - soil_moisture[i]) * layer_thickness * water_density)
        rate_limit = uconvert(u"kg/m^2", Unitful.gn * saturated_conductivity[i] * moisture_timestep)
        infiltrated = min(remaining, layer_capacity, rate_limit)
        soil_moisture[i] += infiltrated / (layer_thickness * water_density)
        remaining -= infiltrated
    end
    return remaining
end

post_infiltration_pool_update(::RateLimitedFrontRainfall, pool, rainfall_flux, moisture_timestep, water_flux, surf_evap, max_surface_pool) =
    clamp(pool - water_flux - surf_evap, 0.0u"kg/m^2", max_surface_pool)
