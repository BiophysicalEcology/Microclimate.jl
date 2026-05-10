"""
    FixedSoilTemperatureIterations(iterations_per_day::Int)

Run a fixed number of iterations of the soil temperature solver per simulated day.
"""
struct FixedSoilTemperatureIterations <: AbstractSoilTemperatureConvergence
    iterations_per_day::Int
end

max_iterations(c::FixedSoilTemperatureIterations) = c.iterations_per_day
may_iterate(c::FixedSoilTemperatureIterations) = c.iterations_per_day > 1
is_converged(::FixedSoilTemperatureIterations, iter, niter, T0, T0_prev) = iter >= niter
