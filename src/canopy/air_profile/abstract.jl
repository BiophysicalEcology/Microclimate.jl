"""
    AbstractCanopyAirProfileModel

How a canopy resolves the in-canopy air-temperature profile from per-layer
leaf sensible-heat sources plus canopy-top and ground boundary values.
Concrete variants implement [`allocate_air_profile`](@ref) and
[`canopy_air_profile!`](@ref). Currently: [`KTheoryAirProfile`](@ref)
(steady first-order-closure diffusion).
"""
abstract type AbstractCanopyAirProfileModel end

"""
    allocate_air_profile(model, canopy_height, plant_area_index, heights, n_layers, boundary_layer_model)

Pre-allocate and precompute the structural buffers a model needs.
"""
function allocate_air_profile end

"""
    canopy_air_profile!(buffers, model, boundary_layer_model; canopy_height, displacement_height,
                         friction_velocity, canopy_top_air_temperature, ground_temperature,
                         sensible_heat_source)

Solve the in-canopy air-temperature profile for this Picard pass; writes
into `buffers`, returns `(; air_temperature)`.
"""
function canopy_air_profile! end

"""
    reset_air_profile_scratch!(model, buffers)

Clear any warm-start/relaxation state that shouldn't carry across a
day-boundary reset (see `reset_canopy_scratch!`) — e.g. Aitken Δ²
residual history, which would otherwise seed the new day's first Picard
pass from an unrelated previous day. No-op by default.
"""
reset_air_profile_scratch!(::AbstractCanopyAirProfileModel, buffers) = nothing
