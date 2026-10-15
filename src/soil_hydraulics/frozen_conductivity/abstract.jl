"""
    AbstractFrozenConductivity

How ice in the pores changes soil hydraulic conductivity.
"""
abstract type AbstractFrozenConductivity end

"""
    conductivity_factor(model, ice_content, porosity)

Fraction of hydraulic conductivity remaining with `ice_content` (m³/m³) in
soil of `porosity` (m³/m³).
"""
function conductivity_factor end

"""
    conductivity_blocked(model, ice_content, porosity)

Whether ice has cut off flow through the layer entirely.
"""
function conductivity_blocked end
