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
