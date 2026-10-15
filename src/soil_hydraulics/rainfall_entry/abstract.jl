abstract type AbstractRainfallEntryMode end

"""
    rainfall_flux_for_step(mode, pool, moisture_timestep, niter_moist)

Rainfall flux (kg/m²/s, ≥ 0) applied as the upper boundary condition of
`infiltration_step!`. Zero for modes that move water with [`rainfall_entry!`](@ref).
"""
function rainfall_flux_for_step end

"""
    rainfall_entry!(soil_moisture, mode, pool; soil_profile, depths, moisture_timestep)

Move water from `pool` into `soil_moisture` before `infiltration_step!`.
Returns the remaining `pool`.
"""
function rainfall_entry! end

"""
    post_infiltration_pool_update(mode, pool, water_flux, surface_evaporation; rainfall_flux, moisture_timestep, max_surface_pool)

Update `pool` after an `infiltration_step!`.
"""
function post_infiltration_pool_update end

# Saturated volumetric water content of soil node `i`
_saturation_water_content(soil_profile, i) =
    1 - soil_profile.bulk_density[i] / soil_profile.mineral_density[i]

# Control-volume thickness of soil node `i`, between the midpoints to its neighbours
@inline function _layer_thickness(depths, i)
    n = length(depths)
    i == 1 && return (depths[2] - depths[1]) / 2
    i == n && return (depths[n] - depths[n-1]) / 2
    return (depths[i+1] - depths[i-1]) / 2
end

"""
    soil_water_storage(soil_moisture, depths)

Total water stored in the soil profile, kg/m².
"""
soil_water_storage(soil_moisture, depths) =
    sum(i -> soil_moisture[i] * _layer_thickness(depths, i), eachindex(soil_moisture)) * water_density
