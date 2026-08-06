"""
    NonlinearSolveCanopyConvergence(; algorithm=SimpleTrustRegion(; autodiff=AutoFiniteDiff()), tolerance=0.01,
                                     max_iterations=50, picard_warmup_iterations=10, picard_warmup_relaxation=0.01,
                                     max_branch_distance=1.0)

Solves the canopy leaf/air-temperature/vapour-density system directly (see
[`_canopy_direct_residual`](@ref)) with a proper globalized nonlinear solver
([`SimpleNonlinearSolve`](https://github.com/SciML/SimpleNonlinearSolve.jl))
instead of [`PicardCanopyConvergence`](@ref)'s hand-tuned relaxation. Default
`algorithm=SimpleTrustRegion(; autodiff=AutoFiniteDiff())` — trust-region
globalization with a finite-difference Jacobian of the actual physics
equations (no nested leaf root-find, no `LinearSolve` cache in the loop).
A future version may switch to an Enzyme-differentiated Jacobian once
`canopy_longwave!` is refactored into a purely-functional kernel -- its
current buffer-mutating form hides derivative information from Enzyme's
activity analysis.

`tolerance` applies to the *scaled* residual norm (state/residual scaled by
[`_CANOPY_TEMPERATURE_SCALE`](@ref) etc., not raw Kelvin) -- 0.01 is a tight
default since the scaled residual blocks are all O(1) by construction.
`max_iterations` plays the same role as `IterationToleranceConvergence`'s
field.

On failure (non-convergence, `DomainError`, or `SingularException` from the
trust-region's internal linear solve), retries once from a better starting
point: `picard_warmup_iterations` cheap, heavily-damped
([`_canopy_picard_pass!`](@ref) at `picard_warmup_relaxation`) passes move
the state into a better basin before handing it back to the solver.

A solver "success" is verified independently, not trusted at face value.
Two checks must both pass: the true residual (recomputed from scratch at
the accepted point) must be small and an improvement over the pre-solve
residual; and the accepted point must lie within `max_branch_distance` of
the pre-solve state (same scaled units as `y`) -- a stiff nonlinear system
can have multiple roots, and a small residual alone can't distinguish a
correct solve from a jump to a different valid root. If either check fails
(including after the Picard warm-up retry), the pre-solve state is kept,
mirroring [`RootFindLeafTemperature`](@ref)'s own fallback convention.
"""
@kwdef struct NonlinearSolveCanopyConvergence{ALG,TOL,MI,PWI,PWR,BD} <: AbstractCanopyConvergenceModel
    algorithm::ALG = SimpleTrustRegion(; autodiff=AutoFiniteDiff())
    tolerance::TOL = 0.01
    max_iterations::MI = 50
    picard_warmup_iterations::PWI = 10
    picard_warmup_relaxation::PWR = 0.01
    max_branch_distance::BD = 1.0
end
