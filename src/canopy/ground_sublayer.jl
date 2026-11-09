"""
    AbstractGroundSublayerModel

The thin layer of air on the ground surface, below the canopy's turbulent mixing.
Concrete models implement [`ground_sublayer_resistance`](@ref):
[`StantonGroundSublayer`](@ref) and [`NoGroundSublayer`](@ref).
"""
abstract type AbstractGroundSublayerModel end

"""
    ground_sublayer_resistance(model, boundary_layer_model; lowest_wind_speed, lowest_layer_height,
        lowest_air_temperature, ground_temperature, atmospheric_pressure, site)

The resistance of the air on the ground surface, as `(; heat, vapor)`. It adds to the air
profile's [`ground_turbulent_resistance`](@ref).
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
    NoGroundSublayer()

No sublayer: the ground exchanges with the lowest layer through turbulent mixing alone.
"""
struct NoGroundSublayer <: AbstractGroundSublayerModel end

# allocate_canopy checks the lowest layer is above the roughness height
function ground_sublayer_resistance(::StantonGroundSublayer, boundary_layer_model;
    lowest_wind_speed, lowest_layer_height, lowest_air_temperature, ground_temperature, atmospheric_pressure, site,
)
    roughness_height = site.roughness_height
    friction_velocity = boundary_layer_model.karman_constant * lowest_wind_speed / log(lowest_layer_height / roughness_height)
    kinematic_viscosity = dry_air_properties((lowest_air_temperature + ground_temperature) / 2,
        atmospheric_pressure).kinematic_viscosity
    heat = 1 / (friction_velocity * sublayer_stanton(roughness_height, friction_velocity, kinematic_viscosity))
    return (; heat, vapor=heat / LEWIS_HEAT_TO_MASS_RATIO)
end

ground_sublayer_resistance(::NoGroundSublayer, boundary_layer_model;
    lowest_wind_speed, lowest_layer_height, lowest_air_temperature, ground_temperature, atmospheric_pressure, site,
) = (; heat=0.0u"s/m", vapor=0.0u"s/m")
