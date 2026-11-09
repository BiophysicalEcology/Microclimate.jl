"""
    AbstractGroundSublayerModel

The thin layer of air on the ground surface, below the canopy's turbulent mixing.
Concrete models implement [`ground_sublayer_resistance`](@ref):
[`StantonGroundSublayer`](@ref), [`SakaguchiZengGroundSublayer`](@ref) and [`NoGroundSublayer`](@ref).
"""
abstract type AbstractGroundSublayerModel end

"""
    ground_sublayer_resistance(model, boundary_layer_model; lowest_wind_speed, lowest_layer_height,
        lowest_air_temperature, ground_temperature, atmospheric_pressure, canopy_height, plant_area_index, site)

The resistance of the air on the ground surface, as `(; heat, vapor)`. It adds to the air
profile's [`ground_turbulent_resistance`](@ref). `plant_area_index` is the canopy's total.
"""
function ground_sublayer_resistance end

"""
    StantonGroundSublayer()

A sublayer with the Stanton number [`sublayer_stanton`](@ref), and a friction velocity from
the log law between the ground's roughness height and the lowest layer's wind. Vapor crosses
it [`LEWIS_HEAT_TO_MASS_RATIO`](@ref) times faster than heat.
"""
struct StantonGroundSublayer <: AbstractGroundSublayerModel end

"""
    SakaguchiZengGroundSublayer(; dense_canopy_transfer_coefficient=0.004, bare_ground_factor=0.13,
        stability_factor=5.0, max_richardson_number=10.0)

The under-canopy ground resistance of Sakaguchi and Zeng (2009), as used in CLM. The transfer
coefficient blends a bare-ground value, `κ / (bare_ground_factor·Re^0.45)` with `Re` the
roughness Reynolds number, and `dense_canopy_transfer_coefficient`, weighted by
`exp(-plant_area_index)`. Air warmer than the ground is stable and divides the dense value by
`1 + stability_factor·Ri`, with `Ri` the Richardson number over the canopy height, at most
`max_richardson_number`. The resistance is `1 / (coefficient · lowest_wind_speed)`, the same
for heat and vapor.

# References
- Sakaguchi, K., & Zeng, X. (2009). Effects of soil wetness, plant litter, and under-canopy
  atmospheric stability on ground evaporation in the Community Land Model (CLM3.5).
  *Journal of Geophysical Research: Atmospheres*, 114(D1).
"""
@kwdef struct SakaguchiZengGroundSublayer{D,B,S,R} <: AbstractGroundSublayerModel
    dense_canopy_transfer_coefficient::D = 0.004
    bare_ground_factor::B = 0.13
    stability_factor::S = 5.0
    max_richardson_number::R = 10.0
end

"""
    NoGroundSublayer()

No sublayer: the ground exchanges with the lowest layer through turbulent mixing alone.
"""
struct NoGroundSublayer <: AbstractGroundSublayerModel end

# allocate_canopy checks the lowest layer is above the roughness height
function ground_sublayer_resistance(::StantonGroundSublayer, boundary_layer_model;
    lowest_wind_speed, lowest_layer_height, lowest_air_temperature, ground_temperature, atmospheric_pressure,
    canopy_height, plant_area_index, site,
)
    roughness_height = site.roughness_height
    friction_velocity = boundary_layer_model.karman_constant * lowest_wind_speed / log(lowest_layer_height / roughness_height)
    kinematic_viscosity = dry_air_properties((lowest_air_temperature + ground_temperature) / 2,
        atmospheric_pressure).kinematic_viscosity
    heat = 1 / (friction_velocity * sublayer_stanton(roughness_height, friction_velocity, kinematic_viscosity))
    return (; heat, vapor=heat / LEWIS_HEAT_TO_MASS_RATIO)
end

function ground_sublayer_resistance(model::SakaguchiZengGroundSublayer, boundary_layer_model;
    lowest_wind_speed, lowest_layer_height, lowest_air_temperature, ground_temperature, atmospheric_pressure,
    canopy_height, plant_area_index, site,
)
    (; dense_canopy_transfer_coefficient, bare_ground_factor, stability_factor, max_richardson_number) = model
    kinematic_viscosity = dry_air_properties((lowest_air_temperature + ground_temperature) / 2,
        atmospheric_pressure).kinematic_viscosity
    roughness_reynolds_number = uconvert(NoUnits, site.roughness_height * lowest_wind_speed / kinematic_viscosity)
    bare_ground_coefficient = boundary_layer_model.karman_constant / (bare_ground_factor * roughness_reynolds_number^0.45)
    dense_coefficient = if lowest_air_temperature > ground_temperature
        richardson_number = uconvert(NoUnits, g_n * canopy_height * (lowest_air_temperature - ground_temperature) /
            (lowest_air_temperature * lowest_wind_speed^2))
        dense_canopy_transfer_coefficient / (1 + stability_factor * min(richardson_number, max_richardson_number))
    else
        dense_canopy_transfer_coefficient
    end
    bare_weight = exp(-plant_area_index)
    transfer_coefficient = bare_weight * bare_ground_coefficient + (1 - bare_weight) * dense_coefficient
    resistance = 1 / (transfer_coefficient * lowest_wind_speed)
    return (; heat=resistance, vapor=resistance)
end

ground_sublayer_resistance(::NoGroundSublayer, boundary_layer_model;
    lowest_wind_speed, lowest_layer_height, lowest_air_temperature, ground_temperature, atmospheric_pressure,
    canopy_height, plant_area_index, site,
) = (; heat=0.0u"s/m", vapor=0.0u"s/m")
