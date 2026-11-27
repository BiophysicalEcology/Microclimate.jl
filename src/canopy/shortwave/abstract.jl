"""
    AbstractCanopyShortwaveModel

Layer-resolved shortwave radiative transfer in a canopy. Concrete variants implement [`allocate_shortwave`](@ref) and
[`canopy_shortwave!`](@ref). Currently: [`TwoStreamRadiation`](@ref)
(Dickinson/Sellers). See also its longwave sibling,
[`AbstractCanopyLongwaveModel`](@ref).
"""
abstract type AbstractCanopyShortwaveModel end

"""
    allocate_shortwave(shortwave_model, canopy_height, plant_area_index, n_layers, canopy_projection_ratio)

Pre-allocate and precompute the structural buffers a model needs.
"""
function allocate_shortwave end
