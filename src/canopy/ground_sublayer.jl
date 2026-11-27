"""
    AbstractGroundSublayerModel

The thin layer of air on the ground surface, below the canopy's turbulent mixing.
Concrete models implement [`ground_sublayer_resistance`](@ref):
[`StantonGroundSublayer`](@ref) and [`NoGroundSublayer`](@ref).
"""
abstract type AbstractGroundSublayerModel end

"""
    StantonGroundSublayer()

A thin layer of air on the ground surface that heat and vapor must cross.
"""
struct StantonGroundSublayer <: AbstractGroundSublayerModel end

"""
    NoGroundSublayer()

No sublayer: the ground exchanges with the lowest layer through turbulent mixing alone.
"""
struct NoGroundSublayer <: AbstractGroundSublayerModel end


"""
    ground_sublayer_resistance(model::AbstractGroundSublayerModel; kw...)

The resistance of the air on the ground surface, as `(; heat, vapor)`.

## Keywords

- friction_velocity: at the ground, from the [`AbstractGroundTurbulenceModel`](@ref)
- roughness_height: of the ground
- lowest_air_temperature
- ground_temperature
- atmospheric_pressure
"""
function ground_sublayer_resistance end
function ground_sublayer_resistance(::StantonGroundSublayer;
    friction_velocity, roughness_height, lowest_air_temperature, ground_temperature, atmospheric_pressure,
)
    kinematic_viscosity = dry_air_properties((lowest_air_temperature + ground_temperature) / 2,
        atmospheric_pressure).kinematic_viscosity
    heat = 1 / (friction_velocity * sublayer_stanton(roughness_height, friction_velocity, kinematic_viscosity))
    return (; heat, vapor=heat / LEWIS_HEAT_TO_MASS_RATIO)
end
ground_sublayer_resistance(::NoGroundSublayer; kw...) = (; heat=0.0u"s/m", vapor=0.0u"s/m")
