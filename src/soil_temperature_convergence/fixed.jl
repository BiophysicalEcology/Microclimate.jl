"""
    FixedIterationConvergence(iterations_per_day::Int)

Run a fixed number of Picard iterations per simulated day/hour. Shared by both
the soil temperature solver's `MicroConfig.convergence` and
`MultilayerCanopy.convergence` (its own leaf-temperature Picard loop) -- the
name is generic on purpose, not soil-specific.
"""
struct FixedIterationConvergence <: AbstractSoilTemperatureConvergence
    iterations_per_day::Int
end

max_iterations(c::FixedIterationConvergence) = c.iterations_per_day
may_iterate(c::FixedIterationConvergence) = c.iterations_per_day > 1
is_converged(::FixedIterationConvergence, iter, niter, T0, T0_prev) = iter >= niter
