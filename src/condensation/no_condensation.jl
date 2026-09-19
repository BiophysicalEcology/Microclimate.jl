"""
    NoCondensation()

Ground dew/frost formation switched off: never condenses, regardless of
surface/air conditions. The model's own bulk-transfer evaporative loss
(driving `evaporation_potential`) is unaffected.
"""
struct NoCondensation <: AbstractCondensationModel end

condensation_energy_flux(::NoCondensation; kw...) = 0.0u"W/m^2"
