"""
    IterationToleranceConvergence(; tolerance, max_iterations_per_day)

Iterate until the max state change between passes is below `tolerance`, or
`max_iterations_per_day` passes have run. Generic -- used by both
`MicroConfig.convergence` and `PicardCanopyConvergence.convergence`.
"""
@kwdef struct IterationToleranceConvergence{T} <: AbstractSoilTemperatureConvergence
    tolerance::T
    max_iterations_per_day::Int
end

max_iterations(c::IterationToleranceConvergence) = c.max_iterations_per_day
may_iterate(::IterationToleranceConvergence) = true

function is_converged(c::IterationToleranceConvergence, iter, niter, T0, T0_prev)
    iter >= niter && return true
    iter <= 1 && return niter <= 1
    tol = ustrip(u"K", c.tolerance)
    max_change = 0.0
    @inbounds for i in eachindex(T0, T0_prev)
        change = abs(ustrip(u"K", T0[i]) - ustrip(u"K", T0_prev[i]))
        change > max_change && (max_change = change)
    end
    return max_change < tol
end
