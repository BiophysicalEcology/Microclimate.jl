"""
    AbstractGroundCoupling

Ground-to-lowest-layer resistance for any [`AbstractCanopyAirProfileModel`](@ref);
variants implement [`ground_resistance`](@ref).
"""
abstract type AbstractGroundCoupling end

"""
    ProfileGroundResistance()

The air-profile model's own ground resistance, computed inside its
`canopy_heat_vapor_transport!` rather than by [`ground_resistance`](@ref).
"""
struct ProfileGroundResistance <: AbstractGroundCoupling end

"""
    GroundExchangeForcing(; layer_height, wind_speed, layer_air_temperature, ground_temperature,
        roughness_height, friction_velocity, plant_area_index, atmospheric_pressure)

Inputs to [`ground_resistance`](@ref). `friction_velocity` is canopy-top;
`layer_height`, `wind_speed` and `layer_air_temperature` are the lowest layer's.
"""
@kwdef struct GroundExchangeForcing{Z,U,TL,TG,Z0,UF,PAI,P}
    layer_height::Z
    wind_speed::U
    layer_air_temperature::TL
    ground_temperature::TG
    roughness_height::Z0
    friction_velocity::UF
    plant_area_index::PAI
    atmospheric_pressure::P
end

"""
    ground_resistance(coupling::AbstractGroundCoupling, boundary_layer_model, forcing::GroundExchangeForcing)

Ground-to-lowest-layer resistance as `(; heat, vapor)`, for every coupling
but [`ProfileGroundResistance`](@ref).
"""
function ground_resistance end

"""
    sublayer_resistance(roughness_height, friction_velocity, kinematic_viscosity)

Quasi-laminar ground sublayer resistance as `(; heat, vapor)`, from
[`sublayer_stanton`](@ref), with vapor through the Lewis factor.
"""
function sublayer_resistance(roughness_height, friction_velocity, kinematic_viscosity)
    heat = 1 / (friction_velocity * sublayer_stanton(roughness_height, friction_velocity, kinematic_viscosity))
    return (; heat, vapor=heat / LEWIS_HEAT_TO_MASS_RATIO)
end

"""
    FloorTurbulence(; gust_decay=4.18, leaf_drag_coefficient=0.25, shelter_factor=0.4,
        free_convection_velocity=0.04u"m/s")

Log-law plus sublayer ([`sublayer_stanton`](@ref)) resistance with floor velocity
`u_g = √((r·u★)² + free_convection_velocity²)`, sheltered by the foliage above:
`r = exp(-gust_decay·leaf_drag_coefficient·PAI/(1 + shelter_factor·PAI))` (Massman & Weil 1999). The log-law takes
Monin-Obukhov stability from the floor's bulk Richardson number: `ζ = Ri` unstable,
`ζ = Ri/(1 - β·Ri)` stable (Dyer 1974).
"""
@kwdef struct FloorTurbulence{G,C,S,U} <: AbstractGroundCoupling
    gust_decay::G = 4.18
    leaf_drag_coefficient::C = 0.25
    shelter_factor::S = 0.4
    free_convection_velocity::U = 0.04u"m/s"
end

function ground_resistance(coupling::FloorTurbulence, boundary_layer_model, forcing::GroundExchangeForcing)
    (; layer_height, layer_air_temperature, ground_temperature, roughness_height, friction_velocity, plant_area_index,
       atmospheric_pressure) = forcing
    (; karman_constant, dyer_constant, stable_beta, turbulent_prandtl_number) = boundary_layer_model
    (; gust_decay, leaf_drag_coefficient, shelter_factor, free_convection_velocity) = coupling

    # Mathematical aliases
    z = layer_height
    z₀ = roughness_height
    Tₐ = layer_air_temperature
    T_g = ground_temperature
    T = (Tₐ + T_g) / 2
    u★ = friction_velocity
    u_fc = free_convection_velocity
    PAI = plant_area_index
    κ = karman_constant
    γ = dyer_constant
    β = stable_beta
    Pr = turbulent_prandtl_number
    ν = dry_air_properties(T, atmospheric_pressure).kinematic_viscosity

    r = exp(-gust_decay * leaf_drag_coefficient * PAI / (1 + shelter_factor * PAI))
    u_g = sqrt((r * u★)^2 + u_fc^2)
    Ri = uconvert(NoUnits, g_n * z * (Tₐ - T_g) / (T * u_g^2))
    # Ri capped below 1/β, where the stable ζ(Ri) diverges
    ζ = Ri < 0 ? Ri : Ri / (1 - β * min(Ri, 0.9 / β))
    ψ_z = calc_ψ_h(ζ, γ, β, Pr)
    ψ_z₀ = calc_ψ_h(ζ * z₀ / z, γ, β, Pr)
    ln_z = log(z / z₀)
    # TODO drop the log-ratio floor with the smooth ψ (#137)
    turbulent_resistance = max(ln_z - ψ_z + ψ_z₀, 0.1 * ln_z) / (κ * u_g)
    sublayer = sublayer_resistance(z₀, u_g, ν)
    return (; heat=turbulent_resistance + sublayer.heat, vapor=turbulent_resistance + sublayer.vapor)
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
    (; layer_air_temperature, ground_temperature, roughness_height, friction_velocity, plant_area_index,
       atmospheric_pressure) = forcing

    # Mathematical aliases
    z₀ = roughness_height
    u★ = friction_velocity
    PAI = plant_area_index
    κ = boundary_layer_model.karman_constant
    C_dense = coupling.dense_canopy_transfer_coefficient
    ν = dry_air_properties((layer_air_temperature + ground_temperature) / 2, atmospheric_pressure).kinematic_viscosity

    C_bare = κ / 0.13 * uconvert(NoUnits, z₀ * u★ / ν)^-0.45
    W = exp(-PAI)
    resistance = 1 / ((C_bare * W + C_dense * (1 - W)) * u★)
    return (; heat=resistance, vapor=resistance)
end
