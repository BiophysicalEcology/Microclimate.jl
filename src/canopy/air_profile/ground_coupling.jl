"""
    AbstractGroundCoupling

Resistance to heat and vapor between the ground and the lowest canopy layer,
for any [`AbstractCanopyAirProfileModel`](@ref). Concrete variants implement
[`ground_resistance`](@ref).
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

Inputs to [`ground_resistance`](@ref): the air-profile model's own resistances,
the lowest canopy layer's height and wind speed, the ground's roughness height,
canopy-top friction velocity and total plant area index, and the air
temperature and pressure for air properties, and von Kármán's constant.
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

Heat and vapor resistance between the ground and the lowest canopy layer, as `(; heat, vapor)`.
"""
function ground_resistance end

ground_resistance(::ProfileGroundResistance, forcing::GroundExchangeForcing) =
    (; heat=forcing.profile_heat_resistance, vapor=forcing.profile_vapor_resistance)
