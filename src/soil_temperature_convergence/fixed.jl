"""
    FixedIterationConvergence(iterations_per_day::Int)

Run a fixed number of Picard iterations per simulated day/hour, for
`MicroConfig.convergence` and `MicroConfig.canopy_soil_convergence`.
"""
struct FixedIterationConvergence <: AbstractSoilTemperatureConvergence
    iterations_per_day::Int
end

max_iterations(c::FixedIterationConvergence) = c.iterations_per_day
may_iterate(c::FixedIterationConvergence) = c.iterations_per_day > 1
is_converged(::FixedIterationConvergence, iter, niter, T0, T0_prev) = iter >= niter
