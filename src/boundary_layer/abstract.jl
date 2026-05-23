abstract type AbstractBoundaryLayerModel end

"""
    allocate_profile(boundary_layer_model, heights)

Allocate the per-variant scratch buffer needed by `atmospheric_surface_profile!`
(height arrays, output profile vectors, solver warm-starts).
"""
function allocate_profile end

"""
    atmospheric_surface_profile(boundary_layer_model; heights, site, environment_instant, surface_temperature, vapour_pressure_equation)

Compute vertical profiles of wind speed, air temperature, and relative
humidity above the surface. Variants supply the stability formulation.
Allocates the output buffer internally.
"""
function atmospheric_surface_profile end

"""
    atmospheric_surface_profile!(boundary_layer_model, buffers; site, environment_instant, surface_temperature, vapour_pressure_equation)

In-place version of `atmospheric_surface_profile` that writes into
pre-allocated buffers.
"""
function atmospheric_surface_profile! end

"""
    surface_fluxes(boundary_layer_model; surface_temperature, air_temperature, wind_speed, zenith_angle, roughness_height, reference_height, obukhov_length_prev)

Surface-layer convective heat flux and friction velocity. Variants apply
whichever atmospheric-stability parameterisation they use. Returns
`(; convective_heat_flux, friction_velocity, ΔT, ρ_cp)`.
"""
function surface_fluxes end

"""
    init_surface_fluxes!(boundary_layer_model, buffers, forcing, site, heights, T0, i)

Refresh any internal state the model needs for `surface_fluxes` on the
upcoming hour `i`. Default no-op for models that don't carry warm-start state.
"""
init_surface_fluxes!(::AbstractBoundaryLayerModel, args...) = nothing
