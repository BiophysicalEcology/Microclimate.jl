"""
    AbstractGroundTurbulenceModel

Turbulent mixing between the lowest canopy layer and the ground. Concrete models implement
[`ground_turbulence`](@ref): [`ProfileGroundTurbulence`](@ref) and [`GroundTurbulence`](@ref).
"""
abstract type AbstractGroundTurbulenceModel end

"""
    ProfileGroundTurbulence()

The air profile's own mixing between its lowest layer and the ground
([`ground_turbulent_resistance`](@ref)). The ground friction velocity is the one the neutral
log law gives for the lowest layer's wind.
"""
struct ProfileGroundTurbulence <: AbstractGroundTurbulenceModel end

"""
    GroundTurbulence(; gust_decay=4.18, leaf_drag_coefficient=0.25, shelter_factor=0.4,
        free_convection_velocity=0.04u"m/s")

A neutral log law for heat from the ground to the lowest layer, at the ground friction velocity
`u_g = √((r·u★)² + free_convection_velocity²)`: the canopy's friction velocity `u★` sheltered by
the foliage above, `r = exp(-gust_decay·leaf_drag_coefficient·PAI/(1 + shelter_factor·PAI))`
(after Massman & Weil 1999). The default coefficients are estimates, not taken from a published table.

# References
- Massman, W. J. and Weil, J. C. (1999). An analytical one-dimensional second-order closure
  model of turbulence statistics and the Lagrangian time scale within and above plant
  canopies of arbitrary structure. *Boundary-Layer Meteorology*, 91, 81-107.
"""
@kwdef struct GroundTurbulence{G,C,S,U} <: AbstractGroundTurbulenceModel
    gust_decay::G = 4.18
    leaf_drag_coefficient::C = 0.25
    shelter_factor::S = 0.4
    free_convection_velocity::U = 0.04u"m/s"
end

"""
    ground_turbulence(model::AbstractGroundTurbulenceModel, air_profile_buffers, air_profile_model, boundary_layer_model; kw...)

The turbulent resistance between the lowest layer and the ground, the same for heat and vapor,
and the friction velocity at the ground, as `(; resistance, friction_velocity)`.

## Keywords

- friction_velocity: the canopy's
- inverse_obukhov_length: the canopy's
- canopy_top_enhancement: the roughness sublayer's [`canopy_top_diffusivity_enhancement`](@ref)
- plant_area_index: of the whole canopy
- lowest_wind_speed
- lowest_layer_height
- lowest_air_temperature
- ground_temperature
- roughness_height: of the ground
"""
function ground_turbulence end
# allocate_canopy checks the lowest layer is above the roughness height
function ground_turbulence(::ProfileGroundTurbulence, air_profile_buffers, air_profile_model, boundary_layer_model;
    friction_velocity, inverse_obukhov_length, canopy_top_enhancement, lowest_wind_speed, lowest_layer_height, roughness_height,
    kw...,
)
    resistance = ground_turbulent_resistance(air_profile_buffers, air_profile_model, boundary_layer_model;
        friction_velocity, inverse_obukhov_length, canopy_top_enhancement)
    ground_friction_velocity = boundary_layer_model.karman_constant * lowest_wind_speed / log(lowest_layer_height / roughness_height)
    return (; resistance, friction_velocity=ground_friction_velocity)
end
function ground_turbulence(model::GroundTurbulence, air_profile_buffers, air_profile_model, boundary_layer_model;
    friction_velocity, plant_area_index, lowest_layer_height, roughness_height, kw...,
)
    (; gust_decay, leaf_drag_coefficient, shelter_factor, free_convection_velocity) = model
    sheltering = exp(-gust_decay * leaf_drag_coefficient * plant_area_index / (1 + shelter_factor * plant_area_index))
    ground_friction_velocity = sqrt((sheltering * friction_velocity)^2 + free_convection_velocity^2)
    resistance = log(lowest_layer_height / roughness_height) / (boundary_layer_model.karman_constant * ground_friction_velocity)
    return (; resistance, friction_velocity=ground_friction_velocity)
end
