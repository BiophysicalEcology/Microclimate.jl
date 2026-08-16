"""
    NoInterception()

Default interception model: rain passes straight through, no leaf-surface
wetness ever accrues. Reproduces pre-canopy rainfall handling bit-for-bit.
"""
struct NoInterception <: AbstractCanopyInterceptionModel end

allocate_interception(::NoInterception, canopy_height, plant_area_index, heights, n_layers, boundary_layer_model) = nothing

canopy_interception!(buffers, ::NoInterception, canopy_projection_ratio; rainfall, wind_speed) =
    (; ground_throughfall = rainfall)

wet_canopy_fraction(::NoInterception, buffers, layer) = 0.0

available_canopy_water(::NoInterception, buffers, layer) = 0.0u"kg/m^2"

snapshot_interception!(::NoInterception, buffers) = nothing

restore_interception!(::NoInterception, buffers) = nothing

reset_interception!(::NoInterception, buffers) = nothing

deplete_canopy_water!(::NoInterception, buffers, layer, evaporated_mass) = nothing
