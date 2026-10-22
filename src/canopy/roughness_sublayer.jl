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

Ratio of the eddy diffusivity at canopy top to its Monin-Obukhov value.
"""
function canopy_top_enhancement end

canopy_top_enhancement(::NoRoughnessSublayer, canopy_height, plant_area_index, ground_roughness_height, karman_constant) = 1.0

"""
    CanopyTopMixingEnhancement(; max_enhancement=3.0, drag_coefficient=0.3, max_friction_velocity_ratio=0.3)

Canopy-top diffusivity `e = 1 + (max_enhancement - 1)·closure` times its Monin-Obukhov value,
fading with height via [`stretched_height`](@ref). `closure` is how far `u★/U_h` has risen from
bare ground toward `max_friction_velocity_ratio` (Raupach 1994).

# References
- Raupach, M. R. (1994). Simplified expressions for vegetation roughness length and zero-plane
  displacement as functions of canopy height and area index. *Boundary-Layer Meteorology*, 71, 211-216.
- Raupach, M. R., Finnigan, J. J., & Brunet, Y. (1996). Coherent eddies and turbulence in vegetation
  canopies: the mixing-layer analogy. *Boundary-Layer Meteorology*, 78, 351-382.
"""
@kwdef struct CanopyTopMixingEnhancement{E,C,R} <: AbstractRoughnessSublayer
    max_enhancement::E = 3.0
    drag_coefficient::C = 0.3
    max_friction_velocity_ratio::R = 0.3
end

function canopy_top_enhancement(model::CanopyTopMixingEnhancement, canopy_height, plant_area_index,
    ground_roughness_height, karman_constant,
)
    (; max_enhancement, drag_coefficient, max_friction_velocity_ratio) = model
    bare = min(karman_constant / log(canopy_height / ground_roughness_height), max_friction_velocity_ratio)
    bare >= max_friction_velocity_ratio && return float(max_enhancement)
    ratio = min(sqrt(bare^2 + drag_coefficient * plant_area_index / 2), max_friction_velocity_ratio)
    closure = (ratio - bare) / (max_friction_velocity_ratio - bare)
    return 1.0 + (max_enhancement - 1.0) * closure
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
    canopy_top_stability_height(column::AboveCanopyColumn)

Height above the displacement height at which stability is evaluated for canopy-top mixing.
"""
canopy_top_stability_height(column::AboveCanopyColumn) = (1.0 + column.enhancement) / 2 * _canopy_top_depth(column)

"""
    stretched_height(column::AboveCanopyColumn, z)

Height above the displacement height at which Monin-Obukhov theory is applied:
`½(z′ + √(z′² + b²))`, `b = (h - d)√(e² - 1)`.
"""
function stretched_height(column::AboveCanopyColumn, z)
    z′ = max(z - column.displacement_height, 1.0e-3u"m")
    b = _canopy_top_depth(column) * sqrt(max(column.enhancement^2 - 1.0, 0.0))
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
    s₁ = stretched_height(column, lower_height)
    s₂ = stretched_height(column, upper_height)
    log_ratio = log(s₂ / s₁)
    ψ₁ = _column_ψ_h(boundary_layer_model, s₁, inverse_obukhov_length)
    ψ₂ = _column_ψ_h(boundary_layer_model, s₂, inverse_obukhov_length)
    # TODO drop the log-ratio floor with the smooth ψ (#137)
    return max(log_ratio - ψ₂ + ψ₁, 0.1 * log_ratio, 1.0e-6) / (boundary_layer_model.karman_constant * friction_velocity)
end
