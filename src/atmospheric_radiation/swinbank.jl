"""
    SwinbankAtmosphericRadiation()

Swinbank's clear-sky longwave estimate (Eq. 10.11 in Campbell & Norman 1998).
Depends on air temperature only.
"""
struct SwinbankAtmosphericRadiation <: AbstractAtmosphericRadiationModel end

function atmospheric_radiation(::SwinbankAtmosphericRadiation, vapour_pressure, air_temperature)
    return uconvert(u"W*m^-2", ((9.2e-6 * (u"K"(air_temperature))^2) * σ * (u"K"(air_temperature))^4) / 1u"K^2")
end
