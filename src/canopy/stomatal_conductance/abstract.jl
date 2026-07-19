"""
    AbstractStomatalConductanceModel

How a canopy resolves the `HeatExchange.LeafEvaporationParameters` fed to
`leaf_heat_balance`/`leaf_temperature`.

- [`PrescribedStomatalConductance`](@ref) (default): one fixed conductance
  in daylight (`zenith_angle < 90°`), a separate (by default cuticular-only)
  conductance at night.
- [`MoistureResponsiveStomatalConductance`](@ref): the same closure curve
  `CampbellSoilHydraulics` uses, applied to a caller-supplied
  `leaf_water_potential`.
- Not yet implemented: a full photosynthesis-coupled demand function
  (Farquhar-style) plugs in at this same slot.
"""
abstract type AbstractStomatalConductanceModel end

"""
    stomatal_conductance(model::AbstractStomatalConductanceModel, zenith_angle, leaf_water_potential)

The `HeatExchange.LeafEvaporationParameters` to use this hour, given the
solar `zenith_angle` and the current `leaf_water_potential`.
"""
function stomatal_conductance end
