abstract type AbstractAtmosphericRadiationModel end

"""
    atmospheric_radiation(model, vapour_pressure, air_temperature)

Estimate downwelling longwave radiation from the clear-sky atmosphere using
the formulation `model`. Returns `(vapour_pressure, atmospheric_longwave)`.
"""
function atmospheric_radiation end
