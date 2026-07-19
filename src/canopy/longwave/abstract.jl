"""
    AbstractCanopyLongwaveModel

How a canopy resolves layer-resolved longwave exchange — the longwave
sibling of [`AbstractCanopyShortwaveModel`](@ref). Concrete variants
implement [`allocate_longwave`](@ref) and [`canopy_longwave!`](@ref).
Currently: [`LayeredLongwaveExchange`](@ref) (sequential gap-fraction
exchange, re-evaluated each Picard pass rather than solved implicitly).
"""
abstract type AbstractCanopyLongwaveModel end

"""
    allocate_longwave(longwave_model, plant_area_index, n_layers)

Pre-allocate and precompute the structural buffers a model needs.
"""
function allocate_longwave end
