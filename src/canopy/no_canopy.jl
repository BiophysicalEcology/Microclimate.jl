"""
    NoCanopy()

Default canopy model: no canopy-layer structure. Vegetation influence is a
scalar `environment_daily.shade` fraction consumed directly by the longwave
(`ViewFactorLongwave`) and soil-energy (`SoilHeatTransport1D`) budgets.
Reproduces existing behaviour bit-for-bit.
"""
struct NoCanopy <: AbstractCanopyModel end

allocate_canopy(::NoCanopy, heights, boundary_layer_model) = nothing

reset_canopy_scratch!(::NoCanopy, buffers) = nothing

allocate_canopy_inputs(::NoCanopy; kw...) = nothing

initial_ground_overrides(::NoCanopy) = (; ground_shortwave_transmission = nothing, ground_incoming_longwave = nothing)

n_canopy_layers(::NoCanopy, heights) = 0
