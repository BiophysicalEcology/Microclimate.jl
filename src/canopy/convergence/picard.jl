"""
    PicardCanopyConvergence(; convergence=IterationToleranceConvergence(; tolerance=0.1u"K", max_iterations_per_day=20),
                               relaxation=0.05)

Hand-rolled, under-relaxed Picard fixed-point iteration for the canopy
leaf/air-temperature solve. An explicit, fully-supported alternative to
[`NonlinearSolveCanopyConvergence`](@ref) (the default) — useful as a
regression baseline or for reproducing this exact numerical behaviour.

- `convergence::AbstractSoilTemperatureConvergence` — reuses the soil
  spin-up loop's own dispatch (`IterationToleranceConvergence`/
  `FixedIterationConvergence`).
- `relaxation` — under-relaxation on each pass's leaf/air-temperature update
  (`x_new = relaxation*x_solved + (1-relaxation)*x_prev`). Free/tunable,
  uncited; prone to needing re-tuning per site/season since it has no
  principled step-size control of its own — see
  [`NonlinearSolveCanopyConvergence`](@ref).
"""
@kwdef struct PicardCanopyConvergence{CV,RL} <: AbstractCanopyConvergenceModel
    convergence::CV = IterationToleranceConvergence(; tolerance=0.1u"K", max_iterations_per_day=20)
    relaxation::RL = 0.05
end
