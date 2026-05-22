abstract type AbstractSoilTemperatureConvergence end

"""
    max_iterations(convergence) -> Int

Upper bound on soil-temperature iterations per day.
"""
function max_iterations end

"""
    may_iterate(convergence) -> Bool

Whether the strategy ever performs more than one iteration per day.
"""
function may_iterate end

"""
    is_converged(convergence, iter, niter, T0, T0_prev) -> Bool

Whether iteration should stop. `iter` is the current iteration, `niter` the
target (from `max_iterations`), `T0`/`T0_prev` the latest two temperature
states.
"""
function is_converged end
