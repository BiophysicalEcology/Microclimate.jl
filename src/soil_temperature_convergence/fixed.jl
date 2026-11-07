"""
    FixedIterationConvergence(iterations_per_day::Int)

Run a fixed number of iterations per simulated day, for `MicroConfig.convergence`.
"""
struct FixedIterationConvergence <: AbstractSoilTemperatureConvergence
    iterations_per_day::Int
end

max_iterations(c::FixedIterationConvergence) = c.iterations_per_day
may_iterate(c::FixedIterationConvergence) = c.iterations_per_day > 1
is_converged(::FixedIterationConvergence, iter, niter, soil_temperatures, previous_soil_temperatures) = iter >= niter
