"""
    AbstractRoughnessSublayer

Mixing in the air just above a canopy; variants implement [`canopy_top_enhancement`](@ref).
"""
abstract type AbstractRoughnessSublayer end

"""
    NoRoughnessSublayer()

Monin-Obukhov theory down to canopy top.
"""
struct NoRoughnessSublayer <: AbstractRoughnessSublayer end

"""
    canopy_top_enhancement(model::AbstractRoughnessSublayer, canopy_height, plant_area_index, ground_roughness_height, karman_constant)

Ratio of the eddy diffusivity at canopy top to its Monin-Obukhov value in fully turbulent flow.
"""
function canopy_top_enhancement end

canopy_top_enhancement(::NoRoughnessSublayer, canopy_height, plant_area_index, ground_roughness_height, karman_constant) = 1.0

"""
    CanopyTopMixingEnhancement(; max_enhancement=3.0, drag_coefficient=0.3, max_friction_velocity_ratio=0.3,
        transition_reynolds_number=1000.0)

Canopy-top diffusivity `e = 1 + (max_enhancement - 1)·closure` times its Monin-Obukhov value,
fading with height via [`stretched_height`](@ref). `closure` is how far `u★/U_h` has risen from
bare ground toward `max_friction_velocity_ratio` (Raupach 1994). Each hour the enhancement
fades by `Re²/(Re² + transition_reynolds_number²)` with canopy Reynolds number `Re = U_h·h/ν`
(0 switches the fade off).

# References
- Raupach, M. R. (1994). Simplified expressions for vegetation roughness length and zero-plane
  displacement as functions of canopy height and area index. *Boundary-Layer Meteorology*, 71, 211-216.
- Raupach, M. R., Finnigan, J. J., & Brunet, Y. (1996). Coherent eddies and turbulence in vegetation
  canopies: the mixing-layer analogy. *Boundary-Layer Meteorology*, 78, 351-382.
"""
@kwdef struct CanopyTopMixingEnhancement{E,C,R,T} <: AbstractRoughnessSublayer
    max_enhancement::E = 3.0
    drag_coefficient::C = 0.3
    max_friction_velocity_ratio::R = 0.3
    transition_reynolds_number::T = 1000.0
end

function canopy_top_enhancement(model::CanopyTopMixingEnhancement, canopy_height, plant_area_index,
    ground_roughness_height, karman_constant,
)
    (; max_enhancement, drag_coefficient, max_friction_velocity_ratio) = model

    # Mathematical aliases
    h = canopy_height
    z₀ = ground_roughness_height
    PAI = plant_area_index
    κ = karman_constant
    e_max = max_enhancement
    C_R = drag_coefficient
    γ_max = max_friction_velocity_ratio

    # u★/U_h over bare ground, then with the canopy (Raupach 1994)
    γ_bare = min(κ / log(h / z₀), γ_max)
    # Bare ground already at the maximum ratio: the canopy closes no gap
    γ_bare >= γ_max && return 1.0
    γ = min(sqrt(γ_bare^2 + C_R * PAI / 2), γ_max)
    closure = (γ - γ_bare) / (γ_max - γ_bare)
    return 1.0 + (e_max - 1.0) * closure
end

"""
    turbulent_enhancement(model::AbstractRoughnessSublayer, enhancement, canopy_height, canopy_top_wind_speed, kinematic_viscosity)

The hour's canopy-top enhancement, from its fully turbulent value `enhancement`.
"""
function turbulent_enhancement end

turbulent_enhancement(::NoRoughnessSublayer, enhancement, canopy_height, canopy_top_wind_speed, kinematic_viscosity) = enhancement

function turbulent_enhancement(model::CanopyTopMixingEnhancement, enhancement, canopy_height, canopy_top_wind_speed,
    kinematic_viscosity,
)
    # Mathematical aliases
    e = enhancement
    h = canopy_height
    U_h = canopy_top_wind_speed
    ν = kinematic_viscosity
    Re_t = model.transition_reynolds_number

    iszero(Re_t) && return e
    Re = uconvert(NoUnits, U_h * h / ν)
    fade = Re^2 / (Re^2 + Re_t^2)
    return 1.0 + (e - 1.0) * fade
end

"""
    AboveCanopyColumn(; canopy_height, displacement_height, friction_velocity, inverse_obukhov_length, enhancement)

The hour's air column above the canopy.
"""
@kwdef struct AboveCanopyColumn{H,D,U,IL,E}
    canopy_height::H
    displacement_height::D
    friction_velocity::U
    inverse_obukhov_length::IL
    enhancement::E
end

# Canopy top above the displacement height, floored
_canopy_top_depth(column::AboveCanopyColumn) = max(column.canopy_height - column.displacement_height, 1.0e-3u"m")

"""
    stretched_height(column::AboveCanopyColumn, z)

Height above the displacement height at which Monin-Obukhov theory is applied:
`½(z′ + √(z′² + b²))`, `b = (h - d)√(e² - 1)`.
"""
function stretched_height(column::AboveCanopyColumn, z)
    d = column.displacement_height
    e = column.enhancement
    z′ = max(z - d, 1.0e-3u"m")
    b = _canopy_top_depth(column) * sqrt(max(e^2 - 1.0, 0.0))
    return (z′ + sqrt(z′^2 + b^2)) / 2
end

_column_ψ_h(boundary_layer_model, height, inverse_obukhov_length) =
    calc_ψ_h(stability_parameter(InverseObukhovLength(), height, inverse_obukhov_length), boundary_layer_model.dyer_constant,
        boundary_layer_model.stable_beta, boundary_layer_model.turbulent_prandtl_number)

"""
    column_heat_resistance(column::AboveCanopyColumn, boundary_layer_model, lower_height, upper_height)

Turbulent heat resistance between two heights at or above canopy top.
"""
function column_heat_resistance(column::AboveCanopyColumn, boundary_layer_model, lower_height, upper_height)
    (; friction_velocity, inverse_obukhov_length) = column

    # Mathematical aliases
    u★ = friction_velocity
    L⁻¹ = inverse_obukhov_length
    κ = boundary_layer_model.karman_constant

    s₁ = stretched_height(column, lower_height)
    s₂ = stretched_height(column, upper_height)
    ln_s = log(s₂ / s₁)
    ψ₁ = _column_ψ_h(boundary_layer_model, s₁, L⁻¹)
    ψ₂ = _column_ψ_h(boundary_layer_model, s₂, L⁻¹)
    # TODO drop the log-ratio floor with the smooth ψ (#137)
    return max(ln_s - ψ₂ + ψ₁, 0.1 * ln_s, 1.0e-6) / (κ * u★)
end
