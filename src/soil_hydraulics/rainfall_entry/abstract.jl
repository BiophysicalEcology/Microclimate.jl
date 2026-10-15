abstract type AbstractRainfallEntryMode end

"""
    RainfallEntryStep(soil_hydraulic_model; soil_profile, depths, pool, moisture_timestep, niter_moist, max_surface_pool, frozen_water_content)

Values shared by every moisture sub-step of one hour of rainfall entry:
the top node's saturated water content and control-volume thickness, the
rainfall boundary flux, the sub-step length, the layers' frozen water
content and the model for its effect on conductivity.
"""
struct RainfallEntryStep{S,H,F,T,D,P,M,I,FC}
    saturation::S
    half_thickness::H
    rainfall_flux::F
    timestep::T
    depths::D
    soil_profile::P
    max_surface_pool::M
    # Fixed for the hour: re-deriving it mid-solve would freeze or thaw water
    # without the matching latent heat
    frozen_water_content::I
    frozen_conductivity::FC
end
function RainfallEntryStep(soil_hydraulic_model;
    soil_profile, depths, pool, moisture_timestep, niter_moist, max_surface_pool, frozen_water_content,
)
    (; rainfall_entry_mode, frozen_conductivity) = soil_hydraulic_model
    (; bulk_density, mineral_density) = soil_profile
    saturation = 1 - bulk_density[1] / mineral_density[1]
    # Node 1's midpoint control volume, matching infiltration_step!'s layer_water_mass[2]
    half_thickness = (depths[2] - depths[1]) / 2
    # Fixed from the hour's starting pool for all sub-steps, reduced by ice in the top node
    rainfall_flux = rainfall_flux_for_step(rainfall_entry_mode, pool, moisture_timestep, niter_moist) *
        conductivity_factor(frozen_conductivity, frozen_water_content[1], saturation)
    return RainfallEntryStep(saturation, half_thickness, rainfall_flux, moisture_timestep, depths,
        soil_profile, max_surface_pool, frozen_water_content, frozen_conductivity)
end

"""
    rainfall_flux_for_step(mode, pool, moisture_timestep, niter_moist)

Rainfall flux (kg/m²/s, ≥ 0) applied as the upper boundary condition of
`infiltration_step!`. Zero for modes that move water with [`rainfall_entry!`](@ref).
"""
function rainfall_flux_for_step end

"""
    rainfall_entry!(soil_moisture, mode, pool, step::RainfallEntryStep)

Move water from `pool` into `soil_moisture` before `infiltration_step!`.
Returns the remaining `pool`.
"""
function rainfall_entry! end

"""
    post_infiltration_pool_update(mode, pool, step::RainfallEntryStep, water_flux, surface_evaporation)

Update `pool` after an `infiltration_step!`.
"""
function post_infiltration_pool_update end

# Control-volume thickness of soil node `i`
@inline function _layer_thickness(depths, half_thickness, i)
    n = length(depths)
    i == 1 && return half_thickness
    i == n && return (depths[n] - depths[n-1]) / 2
    return (depths[i+1] - depths[i-1]) / 2
end

"""
    soil_water_storage(soil_moisture, step::RainfallEntryStep)

Total water stored in the soil profile, kg/m².
"""
soil_water_storage(soil_moisture, step::RainfallEntryStep) =
    sum(i -> soil_moisture[i] * _layer_thickness(step.depths, step.half_thickness, i), eachindex(soil_moisture)) * water_density
