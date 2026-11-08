"""
    AbstractStomatalConductanceModel

How a canopy sets the leaf conductances given to `leaf_heat_balance`.

- [`PrescribedStomatalConductance`](@ref) (default): one fixed conductance
  in daylight (`zenith_angle < 90°`), a separate (by default cuticular-only)
  conductance at night.
- [`MoistureResponsiveStomatalConductance`](@ref): the same closure curve
  `CampbellSoilHydraulics` uses, applied to a caller-supplied
  `leaf_water_potential`.
"""
abstract type AbstractStomatalConductanceModel end

"""
    stomatal_conductance(model::AbstractStomatalConductanceModel, zenith_angle,
                          leaf_water_potential, soil_hydraulic_model)

The `HeatExchange.LeafEvaporationParameters` to use this hour, given the
solar `zenith_angle` and the current `leaf_water_potential`.
`soil_hydraulic_model` is passed for models that read its parameters.
"""
function stomatal_conductance end
