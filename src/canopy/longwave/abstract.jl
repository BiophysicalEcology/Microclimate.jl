"""
    AbstractCanopyLongwaveModel

Longwave exchange through canopy layers: [`LayeredLongwaveExchange`](@ref),
[`AllPairsLongwaveExchange`](@ref) or [`LayeredRadiosityExchange`](@ref).
"""
abstract type AbstractCanopyLongwaveModel end

"""
    allocate_longwave(longwave_model, plant_area_index, n_layers, canopy_projection_ratio)

Allocate a model's buffers. Models that depend on leaf angle use `canopy_projection_ratio`.
"""
function allocate_longwave end
