"""
    AbstractCanopyLongwaveModel

Layer-resolved canopy longwave exchange — the longwave
sibling of [`AbstractCanopyShortwaveModel`](@ref). Concrete variants
implement [`allocate_longwave`](@ref) and [`canopy_longwave!`](@ref).
Currently: [`LayeredLongwaveExchange`](@ref) (sequential, non-reflecting
gap-fraction cascade, O(n_layers)), [`AllPairsLongwaveExchange`](@ref)
(pairwise exchange-weight matrix, O(n_layers^2), ported from
micropoint/microclimlearn for direct comparison), and
[`LayeredRadiosityExchange`](@ref) (implicit reflecting two-stream slabs,
O(n_layers), following Flerchinger's SHAW model). Each is re-evaluated
every Picard pass rather than cached across iterations.
"""
abstract type AbstractCanopyLongwaveModel end

"""
    allocate_longwave(longwave_model, plant_area_index, n_layers, canopy_projection_ratio)

Pre-allocate and precompute the structural buffers a model needs.
`canopy_projection_ratio` (Campbell's ellipsoidal leaf-angle `x`, from
[`LeafParameters`](@ref)) is threaded through for models whose exchange
weights are leaf-angle dependent; models that don't need it just ignore it.
"""
function allocate_longwave end
