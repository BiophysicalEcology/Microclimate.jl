abstract type AbstractAtmosphericRadiationModel end

"""
    atmospheric_radiation(model, vapour_pressure, air_temperature)

Estimate downwelling longwave radiation from the clear-sky atmosphere using
the formulation `model`. Returns `atmospheric_longwave` (W/m²).
"""
function atmospheric_radiation end
