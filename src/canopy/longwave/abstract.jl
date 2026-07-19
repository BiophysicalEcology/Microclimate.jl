"""
    AbstractCanopyLongwaveModel

How a canopy resolves layer-resolved longwave exchange. A canopy model
(e.g. [`MultilayerCanopy`](@ref)) holds one of these in its `longwave_model`
field — the longwave sibling of [`AbstractCanopyShortwaveModel`](@ref).

Concrete variants implement [`allocate_longwave`](@ref) and
[`canopy_longwave!`](@ref). Currently: [`LayeredLongwaveExchange`](@ref)
(sequential gap-fraction exchange, re-evaluated with the latest leaf
temperatures each Picard pass rather than solved as one implicit multi-layer
system).
"""
abstract type AbstractCanopyLongwaveModel end

"""
    allocate_longwave(longwave_model, plant_area_index, n_layers)

Pre-allocate and precompute the structural (once-per-run) buffers a
[`AbstractCanopyLongwaveModel`](@ref) needs.
"""
function allocate_longwave end
