"""
    NoCanopy()

Default canopy model: no canopy-layer structure. Vegetation influence stays
exactly as it is today — a scalar `environment_daily.shade` fraction consumed
directly by the longwave (`ViewFactorLongwave`) and soil-energy
(`SoilHeatTransport1D`) budgets. Choosing `NoCanopy()` reproduces existing
behaviour bit-for-bit; nothing in the solve loop reads `allocate_canopy`
output for this model.
"""
struct NoCanopy <: AbstractCanopyModel end

allocate_canopy(::NoCanopy, heights, boundary_layer_model) = nothing
