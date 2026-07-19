"""
    AbstractCanopyWindModel

How a canopy resolves the in-canopy wind-speed profile below canopy top
(above canopy top is always the shared `AbstractBoundaryLayerModel`, e.g.
`MoninObukhov`, extended with the canopy's displacement height and
roughness length).

Concrete variants implement [`allocate_wind`](@ref) and
[`canopy_wind_profile!`](@ref). Currently: [`CanopyWindAttenuation`](@ref).
"""
abstract type AbstractCanopyWindModel end

"""
    allocate_wind(wind_model, canopy_height, plant_area_index, boundary_layer_model, heights, n_layers)

Pre-allocate and precompute the structural buffers a model needs.
"""
function allocate_wind end
