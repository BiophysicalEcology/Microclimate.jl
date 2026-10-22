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
    AboveCanopyColumn(; canopy_height, displacement_height, roughness_length,
        friction_velocity, inverse_obukhov_length, enhancement)

The hour's air column above the canopy.
"""
@kwdef struct AboveCanopyColumn{H,D,Z0,U,IL,E}
    canopy_height::H
    displacement_height::D
    roughness_length::Z0
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
