"""
    AbstractCanopyConvergenceModel

How `canopy_energy_balance!` drives its per-hour leaf-temperature solve to
convergence. See [`PicardCanopyConvergence`](@ref) and
[`NonlinearSolveCanopyConvergence`](@ref).
"""
abstract type AbstractCanopyConvergenceModel end
