"""
    AbstractCanopyConvergenceModel

How `canopy_energy_balance!` solves its per-hour leaf/air-temperature state.
See [`NewtonCanopyConvergence`](@ref).
"""
abstract type AbstractCanopyConvergenceModel end
