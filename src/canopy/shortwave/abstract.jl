"""
    AbstractCanopyShortwaveModel

How a canopy resolves layer-resolved shortwave radiative transfer. A canopy
model (e.g. [`MultilayerCanopy`](@ref)) holds one of these in its
`shortwave_model` field; canopy-height/plant-area-index are supplied
separately by the parent (they're solver geometry shared with every other
sub-model, not a shortwave-specific concept).

Concrete variants implement [`allocate_shortwave`](@ref) and
[`canopy_shortwave!`](@ref). Currently: [`TwoStreamRadiation`](@ref)
(Dickinson/Sellers). A simpler alternative (e.g. Beer-Lambert extinction,
in the style of Monteith's big-leaf canopy model) is a natural addition at
this same slot. See also [`AbstractCanopyLongwaveModel`](@ref), its
longwave sibling.
"""
abstract type AbstractCanopyShortwaveModel end

"""
    allocate_shortwave(shortwave_model, canopy_height, plant_area_index, n_layers)

Pre-allocate and precompute the structural (once-per-run) buffers a
[`AbstractCanopyShortwaveModel`](@ref) needs.
"""
function allocate_shortwave end
