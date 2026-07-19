"""
    NoInterception()

Default interception model: rain passes straight through the canopy
untouched, and no leaf-surface wetness ever accrues. Reproduces the
package's previous (pre-canopy) rainfall handling bit-for-bit for
`MultilayerCanopy` users who don't opt into interception.
"""
struct NoInterception <: AbstractCanopyInterceptionModel end

allocate_interception(::NoInterception, canopy_height, plant_area_index, heights, n_layers, boundary_layer_model) = nothing

canopy_interception!(buffers, ::NoInterception, leaf_angle_distribution_parameter; rainfall, wind_speed) =
    (; ground_throughfall = rainfall)

wet_canopy_fraction(::NoInterception, buffers, layer) = 0.0

deplete_canopy_water!(::NoInterception, buffers, layer, evaporated_mass) = nothing
