"""
    AbstractGroundCoupling

Ground-to-lowest-layer resistance for any [`AbstractCanopyAirProfileModel`](@ref);
variants implement [`ground_resistance`](@ref).
"""
abstract type AbstractGroundCoupling end

"""
    ProfileGroundResistance()

The air-profile model's own ground resistance.
"""
struct ProfileGroundResistance <: AbstractGroundCoupling end

"""
    GroundExchangeForcing(; profile_heat_resistance, profile_vapor_resistance, layer_height, wind_speed,
        roughness_height, friction_velocity, plant_area_index, air_temperature, atmospheric_pressure, karman_constant)

Inputs to [`ground_resistance`](@ref). `friction_velocity` is canopy-top;
`layer_height` and `wind_speed` are the lowest layer's.
"""
@kwdef struct GroundExchangeForcing{RH,RV,Z,U,Z0,UF,PAI,T,P,K}
    profile_heat_resistance::RH
    profile_vapor_resistance::RV
    layer_height::Z
    wind_speed::U
    roughness_height::Z0
    friction_velocity::UF
    plant_area_index::PAI
    air_temperature::T
    atmospheric_pressure::P
    karman_constant::K
end

"""
    ground_resistance(coupling::AbstractGroundCoupling, forcing::GroundExchangeForcing)

Ground-to-lowest-layer resistance as `(; heat, vapor)`.
"""
function ground_resistance end

ground_resistance(::ProfileGroundResistance, forcing::GroundExchangeForcing) =
    (; heat=forcing.profile_heat_resistance, vapor=forcing.profile_vapor_resistance)

"""
    FloorTurbulence(; floor_velocity_ratio=0.2, free_convection_velocity=0.04u"m/s")

Log-law plus sublayer ([`sublayer_stanton`](@ref)) resistance with floor velocity
`u_g = √((floor_velocity_ratio·u★)² + free_convection_velocity²)`. Neutral.
"""
@kwdef struct FloorTurbulence{R,U} <: AbstractGroundCoupling
    floor_velocity_ratio::R = 0.2
    free_convection_velocity::U = 0.04u"m/s"
end

function ground_resistance(coupling::FloorTurbulence, forcing::GroundExchangeForcing)
    (; layer_height, roughness_height, friction_velocity, air_temperature, atmospheric_pressure, karman_constant) = forcing
    u_g = sqrt((coupling.floor_velocity_ratio * friction_velocity)^2 + coupling.free_convection_velocity^2)
    kinematic_viscosity = dry_air_properties(air_temperature, atmospheric_pressure).kinematic_viscosity
    turbulent = log(layer_height / roughness_height) / (karman_constant * u_g)
    sublayer = 1 / (u_g * sublayer_stanton(roughness_height, u_g, kinematic_viscosity))
    return (; heat=turbulent + sublayer, vapor=turbulent + sublayer / LEWIS_HEAT_TO_MASS_RATIO)
end

"""
    SakaguchiZeng2009(; dense_canopy_transfer_coefficient=0.004)

`1/(C_s u★)`, `C_s = C_bare·W + dense_canopy_transfer_coefficient·(1 - W)`, `W = exp(-PAI)`,
`C_bare = (κ/0.13)(z0 u★/ν)^-0.45` (Sakaguchi & Zeng 2009). Neutral.
"""
@kwdef struct SakaguchiZeng2009{C} <: AbstractGroundCoupling
    dense_canopy_transfer_coefficient::C = 0.004
end

function ground_resistance(coupling::SakaguchiZeng2009, forcing::GroundExchangeForcing)
    (; roughness_height, friction_velocity, plant_area_index, air_temperature, atmospheric_pressure, karman_constant) = forcing
    kinematic_viscosity = dry_air_properties(air_temperature, atmospheric_pressure).kinematic_viscosity
    bare = karman_constant / 0.13 * uconvert(NoUnits, roughness_height * friction_velocity / kinematic_viscosity)^-0.45
    W = exp(-plant_area_index)
    resistance = 1 / ((bare * W + coupling.dense_canopy_transfer_coefficient * (1 - W)) * friction_velocity)
    return (; heat=resistance, vapor=resistance)
end
