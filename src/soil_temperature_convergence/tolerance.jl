"""
    IterationToleranceConvergence(; tolerance, max_iterations_per_day)

Iterate until the max state change between passes is below `tolerance`, or
`max_iterations_per_day` passes have run. Used by `MicroConfig.convergence`.
"""
@kwdef struct IterationToleranceConvergence{T} <: AbstractSoilTemperatureConvergence
    tolerance::T
    max_iterations_per_day::Int
end

max_iterations(c::IterationToleranceConvergence) = c.max_iterations_per_day
may_iterate(::IterationToleranceConvergence) = true

function is_converged(c::IterationToleranceConvergence, iter, niter, soil_temperatures, previous_soil_temperatures)
    iter >= niter && return true
    iter <= 1 && return niter <= 1
    tol = ustrip(u"K", c.tolerance)
    max_change = 0.0
    for i in eachindex(soil_temperatures, previous_soil_temperatures)
        change = abs(ustrip(u"K", soil_temperatures[i]) - ustrip(u"K", previous_soil_temperatures[i]))
        change > max_change && (max_change = change)
    end
    return max_change < tol
end
