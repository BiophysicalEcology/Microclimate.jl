abstract type AbstractDiffuseFractionModel end

"""
    calc_diffuse_fraction(model, clearness_index)

Diffuse fraction of global radiation as a function of clearness index (the
ratio of global to extraterrestrial irradiance on a horizontal plane).
"""
function calc_diffuse_fraction end
