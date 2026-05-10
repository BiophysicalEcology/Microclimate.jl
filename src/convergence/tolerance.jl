"""
    SoilTemperatureConvergenceTolerance(; tolerance, max_iterations_per_day)

Iterate the soil temperature solver until the maximum nodal temperature change
between successive full-day passes falls below `tolerance`, or until
`max_iterations_per_day` passes have been completed.
"""
@kwdef struct SoilTemperatureConvergenceTolerance{T} <: AbstractSoilTemperatureConvergence
    tolerance::T
    max_iterations_per_day::Int
end

max_iterations(c::SoilTemperatureConvergenceTolerance) = c.max_iterations_per_day
may_iterate(::SoilTemperatureConvergenceTolerance) = true

function is_converged(c::SoilTemperatureConvergenceTolerance, iter, niter, T0, T0_prev)
    iter >= niter && return true
    iter <= 1 && return niter <= 1
    tol = ustrip(u"K", c.tolerance)
    max_change = maximum(abs.(ustrip.(u"K", T0) .- ustrip.(u"K", T0_prev)))
    return max_change < tol
end
