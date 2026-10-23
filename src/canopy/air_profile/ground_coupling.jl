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
        layer_air_temperature, ground_temperature, roughness_height, friction_velocity, plant_area_index,
        atmospheric_pressure)

Inputs to [`ground_resistance`](@ref). `friction_velocity` is canopy-top;
`layer_height`, `wind_speed` and `layer_air_temperature` are the lowest layer's.
"""
@kwdef struct GroundExchangeForcing{RH,RV,Z,U,TL,TG,Z0,UF,PAI,P}
    profile_heat_resistance::RH
    profile_vapor_resistance::RV
    layer_height::Z
    wind_speed::U
    layer_air_temperature::TL
    ground_temperature::TG
    roughness_height::Z0
    friction_velocity::UF
    plant_area_index::PAI
    atmospheric_pressure::P
end

_mean_air_temperature(forcing::GroundExchangeForcing) = (forcing.layer_air_temperature + forcing.ground_temperature) / 2

"""
    ground_resistance(coupling::AbstractGroundCoupling, boundary_layer_model, forcing::GroundExchangeForcing)

Ground-to-lowest-layer resistance as `(; heat, vapor)`.
"""
function ground_resistance end

ground_resistance(::ProfileGroundResistance, boundary_layer_model, forcing::GroundExchangeForcing) =
    (; heat=forcing.profile_heat_resistance, vapor=forcing.profile_vapor_resistance)

"""
    FloorTurbulence(; floor_velocity_ratio=0.2, free_convection_velocity=0.04u"m/s")

Log-law plus sublayer ([`sublayer_stanton`](@ref)) resistance with floor velocity
`u_g = √((floor_velocity_ratio·u★)² + free_convection_velocity²)`. The log-law takes
Monin-Obukhov stability from the floor's bulk Richardson number: `ζ = Ri` unstable,
`ζ = Ri/(1 - β·Ri)` stable (Dyer 1974).
"""
@kwdef struct FloorTurbulence{R,U} <: AbstractGroundCoupling
    floor_velocity_ratio::R = 0.2
    free_convection_velocity::U = 0.04u"m/s"
end

function ground_resistance(coupling::FloorTurbulence, boundary_layer_model, forcing::GroundExchangeForcing)
    (; layer_height, layer_air_temperature, ground_temperature, roughness_height, friction_velocity, atmospheric_pressure) = forcing
    (; karman_constant, dyer_constant, stable_beta, turbulent_prandtl_number) = boundary_layer_model
    air_temperature = _mean_air_temperature(forcing)
    u_g = sqrt((coupling.floor_velocity_ratio * friction_velocity)^2 + coupling.free_convection_velocity^2)
    kinematic_viscosity = dry_air_properties(air_temperature, atmospheric_pressure).kinematic_viscosity

    richardson_number = uconvert(NoUnits,
        g_n * layer_height * (layer_air_temperature - ground_temperature) / (air_temperature * u_g^2))
    # Ri capped below 1/β, where the stable ζ(Ri) diverges
    ζ = richardson_number < 0 ? richardson_number :
        richardson_number / (1 - stable_beta * min(richardson_number, 0.9 / stable_beta))
    ψ_top = calc_ψ_h(ζ, dyer_constant, stable_beta, turbulent_prandtl_number)
    ψ_bottom = calc_ψ_h(ζ * roughness_height / layer_height, dyer_constant, stable_beta, turbulent_prandtl_number)
    log_ratio = log(layer_height / roughness_height)
    # TODO drop the log-ratio floor with the smooth ψ (#137)
    turbulent = max(log_ratio - ψ_top + ψ_bottom, 0.1 * log_ratio) / (karman_constant * u_g)
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

function ground_resistance(coupling::SakaguchiZeng2009, boundary_layer_model, forcing::GroundExchangeForcing)
    (; roughness_height, friction_velocity, plant_area_index, atmospheric_pressure) = forcing
    kinematic_viscosity = dry_air_properties(_mean_air_temperature(forcing), atmospheric_pressure).kinematic_viscosity
    bare = boundary_layer_model.karman_constant / 0.13 * uconvert(NoUnits, roughness_height * friction_velocity / kinematic_viscosity)^-0.45
    W = exp(-plant_area_index)
    resistance = 1 / ((bare * W + coupling.dense_canopy_transfer_coefficient * (1 - W)) * friction_velocity)
    return (; heat=resistance, vapor=resistance)
end
