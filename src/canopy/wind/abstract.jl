"""
    AbstractCanopyWindModel

How wind decays below the canopy top; above it is the boundary-layer model. See
[`ExponentialCanopyWindAttenuation`](@ref) (default) and [`MixingLengthCanopyWindAttenuation`](@ref).
"""
abstract type AbstractCanopyWindModel end

"""
    allocate_wind(wind_model, canopy_height, plant_area_index, boundary_layer_model, heights, n_layers)

Pre-allocate and precompute the structural buffers a model needs.
"""
function allocate_wind end
