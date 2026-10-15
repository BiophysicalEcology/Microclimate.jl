"""
    NoFrozenConductivity()

Ice has no effect on hydraulic conductivity.
"""
struct NoFrozenConductivity <: AbstractFrozenConductivity end

conductivity_factor(::NoFrozenConductivity, ice_content, porosity) = 1.0
conductivity_blocked(::NoFrozenConductivity, ice_content, porosity) = false
