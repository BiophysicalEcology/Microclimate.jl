"""
    AbstractRoughnessSublayer

Mixing in the air just above a canopy; variants implement [`canopy_top_diffusivity_enhancement`](@ref).
"""
abstract type AbstractRoughnessSublayer end

"""
    NoRoughnessSublayer()

Monin-Obukhov theory down to canopy top.
"""
struct NoRoughnessSublayer <: AbstractRoughnessSublayer end

"""
    CanopyTopMixingEnhancement(; max_enhancement=3.0, drag_coefficient=0.3, max_friction_velocity_ratio=0.3,
        transition_reynolds_number=1000.0)

Canopy-top diffusivity `e = 1 + (max_enhancement - 1)·closure` times its Monin-Obukhov value,
fading with height via [`stretched_height`](@ref); `closure` is how far `u★/U_h` has risen from
bare ground to `max_friction_velocity_ratio` (Raupach 1994). Scaled by
`Re²/(Re² + transition_reynolds_number²)`, `Re = U_h·h/ν`.

The defaults are micropoint's: `max_enhancement` is 3 over a closed canopy (Raupach et al. 1996),
`drag_coefficient` and `max_friction_velocity_ratio` are Raupach's (1994), and
`transition_reynolds_number` is micropoint's order-of-magnitude value.

# References
- Raupach, M. R. (1994). Simplified expressions for vegetation roughness length and zero-plane
  displacement as functions of canopy height and area index. *Boundary-Layer Meteorology*, 71, 211-216.
- Raupach, M. R., Finnigan, J. J., & Brunet, Y. (1996). Coherent eddies and turbulence in vegetation
  canopies: the mixing-layer analogy. *Boundary-Layer Meteorology*, 78, 351-382.
- micropoint, `src/micropoint_new2.cpp`, `CANOPY_TOP_ENHANCEMENT` and `CANOPY_RE_TRANSITION`
  (https://github.com/ilyamaclean/micropoint).
"""
@kwdef struct CanopyTopMixingEnhancement{E,C,R,T} <: AbstractRoughnessSublayer
    max_enhancement::E = 3.0
    drag_coefficient::C = 0.3
    max_friction_velocity_ratio::R = 0.3
    transition_reynolds_number::T = 1000.0
end

"""
    canopy_top_diffusivity_enhancement(model::AbstractRoughnessSublayer; canopy_height, plant_area_index,
        ground_roughness_height, karman_constant, canopy_top_wind_speed, kinematic_viscosity)

Ratio of the eddy diffusivity at canopy top to its Monin-Obukhov value.
"""
function canopy_top_diffusivity_enhancement end

canopy_top_diffusivity_enhancement(::NoRoughnessSublayer; kw...) = 1.0

function canopy_top_diffusivity_enhancement(model::CanopyTopMixingEnhancement; canopy_height, plant_area_index,
    ground_roughness_height, karman_constant, canopy_top_wind_speed, kinematic_viscosity,
)
    (; max_enhancement, drag_coefficient, max_friction_velocity_ratio, transition_reynolds_number) = model
    # u★/U_h over bare ground, then with the canopy (Raupach 1994)
    bare_ratio = min(karman_constant / log(canopy_height / ground_roughness_height), max_friction_velocity_ratio)
    # Bare ground already at the maximum ratio leaves no gap for the canopy to close
    bare_ratio >= max_friction_velocity_ratio && return 1.0
    canopy_ratio = min(sqrt(bare_ratio^2 + drag_coefficient * plant_area_index / 2), max_friction_velocity_ratio)
    closure = (canopy_ratio - bare_ratio) / (max_friction_velocity_ratio - bare_ratio)
    turbulent_enhancement = 1.0 + (max_enhancement - 1.0) * closure
    reynolds_number = uconvert(NoUnits, canopy_top_wind_speed * canopy_height / kinematic_viscosity)
    fade = reynolds_number^2 / (reynolds_number^2 + transition_reynolds_number^2)
    return 1.0 + (turbulent_enhancement - 1.0) * fade
end

# Keeps log and power laws finite at the displacement height
const MIN_HEIGHT_ABOVE_DISPLACEMENT = 1.0e-3u"m"

height_above_displacement(height, displacement_height) = max(height - displacement_height, MIN_HEIGHT_ABOVE_DISPLACEMENT)

canopy_top_depth(canopy_height, displacement_height) = height_above_displacement(canopy_height, displacement_height)

"""
    stretched_height(height, canopy_height, displacement_height, enhancement)

Height above the displacement height at which Monin-Obukhov theory is applied:
`½(z′ + √(z′² + b²))`, with `z′` the height above the displacement height and
`b = (h - d)√(e² - 1)`.
"""
function stretched_height(height, canopy_height, displacement_height, enhancement)
    above_displacement = height_above_displacement(height, displacement_height)
    stretch = canopy_top_depth(canopy_height, displacement_height) * sqrt(max(enhancement^2 - 1.0, 0.0))
    return (above_displacement + sqrt(above_displacement^2 + stretch^2)) / 2
end
