"""
    AbstractCanopyAirProfileModel

How a canopy resolves the in-canopy air-temperature profile from per-layer
leaf sensible-heat sources plus canopy-top and ground boundary values. A
canopy model (e.g. [`MultilayerCanopy`](@ref)) holds one of these in its
`air_profile_model` field.

Concrete variants implement [`allocate_air_profile`](@ref) and
[`canopy_air_profile!`](@ref). Currently: [`KTheoryAirProfile`](@ref)
(steady first-order-closure diffusion). The C++ reference's transient
Lagrangian near/far-field model (Raupach) is a more detailed alternative at
this same slot, should the steady-state simplification prove too coarse.
"""
abstract type AbstractCanopyAirProfileModel end

"""
    allocate_air_profile(model, canopy_height, plant_area_index, heights, n_layers, boundary_layer_model)

Pre-allocate and precompute the structural (once-per-run) buffers a
[`AbstractCanopyAirProfileModel`](@ref) needs.
"""
function allocate_air_profile end

"""
    canopy_air_profile!(buffers, model, boundary_layer_model; canopy_height, displacement_height,
                         friction_velocity, canopy_top_air_temperature, ground_temperature,
                         sensible_heat_source)

Solve for the in-canopy air-temperature profile this Picard pass, writing
into `buffers` and returning `(; air_temperature)`.
"""
function canopy_air_profile! end
