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
    AbstractGroundStability

The stability that corrects [`GroundTurbulence`](@ref)'s log law:
[`GroundRichardson`](@ref), [`LouisGroundRichardson`](@ref) or [`CanopyObukhovLength`](@ref).
"""
abstract type AbstractGroundStability end

"""
    GroundRichardson()

The stability parameter [`bulk_richardson_stability_parameter`](@ref) gives for the bulk
Richardson number of the air between the ground and the lowest layer, at the lowest layer's wind.
"""
struct GroundRichardson <: AbstractGroundStability end

"""
    CanopyObukhovLength()

The stability of the whole canopy's flux, at the canopy's Obukhov length, as micropoint does.
"""
struct CanopyObukhovLength <: AbstractGroundStability end

"""
    LouisGroundRichardson(; stability_coefficient=9.4, heat_coefficient=5.3)

Heat transfer scaled from neutral by Louis's (1979) function of the bulk Richardson number `Ri` of the air between the
ground and the lowest layer, with no critical `Ri`: `F = (1 + b·Ri/2)⁻²` when stable, and
`F = 1 - b·Ri/(1 + heat_coefficient·b·(κ/ln(z/z₀))²·√(-Ri·z/z₀))` when unstable. Both sides fall at the same rate
through neutral. `b` is `stability_coefficient`. Neutral is the log law, without Louis's `1/R`
factor for heat, so that every stability option shares it.

# References
- Louis, J.-F. (1979). A parametric model of vertical eddy fluxes in the atmosphere.
  *Boundary-Layer Meteorology*, 17(2), 187-202.
"""
@kwdef struct LouisGroundRichardson{B,C} <: AbstractGroundStability
    stability_coefficient::B = 9.4
    heat_coefficient::C = 5.3
end

"""
    GroundTurbulence(; gust_decay=4.18, leaf_drag_coefficient=0.25, shelter_factor=0.4,
        free_convection_velocity=0.04u"m/s", stability=GroundRichardson())

A log law for heat from the ground to the lowest layer, at the ground friction velocity
`u_g = √((r·u★)² + free_convection_velocity²)`: the canopy's friction velocity `u★` sheltered by
the foliage above, `r = exp(-gust_decay·leaf_drag_coefficient·PAI/(1 + shelter_factor·PAI))`
(after Massman & Weil 1999). The default coefficients are estimates, not taken from a published table.
`stability` ([`GroundRichardson`](@ref) or [`CanopyObukhovLength`](@ref))
corrects the log law.

# References
- Massman, W. J. and Weil, J. C. (1999). An analytical one-dimensional second-order closure
  model of turbulence statistics and the Lagrangian time scale within and above plant
  canopies of arbitrary structure. *Boundary-Layer Meteorology*, 91, 81-107.
"""
@kwdef struct GroundTurbulence{G,C,S,U,ST<:AbstractGroundStability} <: AbstractGroundTurbulenceModel
    gust_decay::G = 4.18
    leaf_drag_coefficient::C = 0.25
    shelter_factor::S = 0.4
    free_convection_velocity::U = 0.04u"m/s"
    stability::ST = GroundRichardson()
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
- ground_roughness_length: of the ground
"""
function ground_turbulence end
# allocate_canopy checks the lowest layer is above the roughness length
function ground_turbulence(::ProfileGroundTurbulence, air_profile_buffers, air_profile_model, boundary_layer_model;
    friction_velocity, inverse_obukhov_length, canopy_top_enhancement, lowest_wind_speed, lowest_layer_height, ground_roughness_length,
    kw...,
)
    resistance = ground_turbulent_resistance(air_profile_buffers, air_profile_model, boundary_layer_model;
        friction_velocity, inverse_obukhov_length, canopy_top_enhancement)
    ground_friction_velocity = boundary_layer_model.karman_constant * lowest_wind_speed / log(lowest_layer_height / ground_roughness_length)
    return (; resistance, friction_velocity=ground_friction_velocity)
end
function ground_turbulence(model::GroundTurbulence, air_profile_buffers, air_profile_model, boundary_layer_model;
    friction_velocity, inverse_obukhov_length, plant_area_index, lowest_wind_speed, lowest_layer_height, lowest_air_temperature,
    ground_temperature, ground_roughness_length, kw...,
)
    (; gust_decay, leaf_drag_coefficient, shelter_factor, free_convection_velocity, stability) = model

    sheltering = exp(-gust_decay * leaf_drag_coefficient * plant_area_index / (1 + shelter_factor * plant_area_index))
    ground_friction_velocity = sqrt((sheltering * friction_velocity)^2 + free_convection_velocity^2)
    # The Richardson number takes the mean wind, with the same floor in calm air
    ground_wind_speed = sqrt(lowest_wind_speed^2 + free_convection_velocity^2)
    resistance = ground_log_resistance(stability, boundary_layer_model; inverse_obukhov_length, lowest_layer_height,
        lowest_air_temperature, ground_temperature, ground_roughness_length, ground_friction_velocity, ground_wind_speed)
    return (; resistance, friction_velocity=ground_friction_velocity)
end

# The log-law resistance for heat from the ground to the lowest layer, corrected for stability
function ground_log_resistance(stability::AbstractGroundStability, boundary_layer_model;
    lowest_layer_height, ground_roughness_length, ground_friction_velocity, kw...,
)
    (; karman_constant) = boundary_layer_model
    stability_parameter = ground_stability_parameter(stability, boundary_layer_model; lowest_layer_height, ground_roughness_length, kw...)
    return (log(lowest_layer_height / ground_roughness_length) - heat_stability_correction(boundary_layer_model, stability_parameter) +
        heat_stability_correction(boundary_layer_model, stability_parameter * ground_roughness_length / lowest_layer_height)) /
        (karman_constant * ground_friction_velocity)
end
function ground_log_resistance(stability::LouisGroundRichardson, boundary_layer_model;
    lowest_layer_height, ground_roughness_length, ground_friction_velocity, ground_wind_speed, lowest_air_temperature, ground_temperature, kw...,
)
    (; stability_coefficient, heat_coefficient) = stability
    (; karman_constant) = boundary_layer_model
    richardson_number = ground_richardson_number(; lowest_layer_height, lowest_air_temperature, ground_temperature, ground_wind_speed)
    neutral_log = log(lowest_layer_height / ground_roughness_length)
    # Heat transfer as a fraction of neutral
    neutral_fraction = richardson_number > 0 ? 1 / (1 + stability_coefficient * richardson_number / 2)^2 :
        1 - stability_coefficient * richardson_number / (1 + heat_coefficient * stability_coefficient * (karman_constant / neutral_log)^2 *
            sqrt(-richardson_number * lowest_layer_height / ground_roughness_length))
    return neutral_log / (karman_constant * ground_friction_velocity * neutral_fraction)
end

# The bulk Richardson number of the air between the ground and the lowest layer, at the
# lowest layer's mean wind speed
function ground_richardson_number(; lowest_layer_height, lowest_air_temperature, ground_temperature, ground_wind_speed)
    mean_temperature = (lowest_air_temperature + ground_temperature) / 2
    return uconvert(NoUnits,
        g_n * lowest_layer_height * (lowest_air_temperature - ground_temperature) / (mean_temperature * ground_wind_speed^2))
end

# The stability parameter at the lowest layer's height
ground_stability_parameter(::CanopyObukhovLength, boundary_layer_model; inverse_obukhov_length, lowest_layer_height, kw...) =
    uconvert(NoUnits, lowest_layer_height * inverse_obukhov_length)
ground_stability_parameter(::GroundRichardson, boundary_layer_model; lowest_layer_height, ground_roughness_length,
    lowest_air_temperature, ground_temperature, ground_wind_speed, kw...) =
    bulk_richardson_stability_parameter(boundary_layer_model,
        ground_richardson_number(; lowest_layer_height, lowest_air_temperature, ground_temperature, ground_wind_speed),
        lowest_layer_height, ground_roughness_length)
