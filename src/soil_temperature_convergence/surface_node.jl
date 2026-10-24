"""
    SurfaceNodePrediction(; tolerance=0.1u"K", max_iterations=20)

Canopy–soil coupling for `MicroConfig.canopy_soil_convergence` that integrates the soil once per
hour. Each pass solves the canopy, then predicts the end-of-hour surface temperature from a
backward-Euler step of the soil surface node, until it changes by less than `tolerance`.
Hours with snow use [`IterationToleranceConvergence`](@ref).
"""
@kwdef struct SurfaceNodePrediction{T} <: AbstractSoilTemperatureConvergence
    tolerance::T = 0.1u"K"
    max_iterations::Int = 20
end
