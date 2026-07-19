"""
    AbstractStomatalConductanceModel

How a canopy resolves the `HeatExchange.LeafEvaporationParameters` fed to
`leaf_heat_balance`/`leaf_temperature`. A canopy model holds one of these in
its `stomatal_model` field.

This is a genuine range of complexity, not just one formula — from fixed/
day-night gating, up to soil-moisture response and a full photosynthesis-
coupled demand function. Each level is its own concrete variant at this
same slot:

- [`PrescribedStomatalConductance`](@ref) (default): one fixed conductance
  in daylight (`zenith_angle < 90°`), a separate (by default cuticular-only)
  conductance at night — using the zenith angle the rest of the model
  already computes, not a new signal.
- [`MoistureResponsiveStomatalConductance`](@ref): the same closure curve
  `CampbellSoilHydraulics` (`soil_hydraulics/campbell.jl`) already uses,
  applied to a caller-supplied `leaf_water_potential`.
- Not yet implemented: a full photosynthesis-coupled demand function
  (Farquhar-style) — plugs in at this same slot without changing
  `leaf_heat_balance` or the leaf-temperature solvers.
"""
abstract type AbstractStomatalConductanceModel end

"""
    stomatal_conductance(model::AbstractStomatalConductanceModel, zenith_angle, leaf_water_potential)

The `HeatExchange.LeafEvaporationParameters` to use this hour, given the
solar `zenith_angle` and the current `leaf_water_potential`.
"""
function stomatal_conductance end
