"""
    CampbellNormanAtmosphericRadiation()

Clear-sky longwave estimate from Campbell & Norman (1998) eq. 10.10. Depends
on both air temperature and vapour pressure.
"""
struct CampbellNormanAtmosphericRadiation <: AbstractAtmosphericRadiationModel end

function atmospheric_radiation(::CampbellNormanAtmosphericRadiation, vapour_pressure, air_temperature)
    return u"W/m^2"((1.72 * (ustrip(u"kPa", vapour_pressure) / ustrip(u"K", air_temperature + 0.01u"K"))^(1//7)) * σ * (u"K"(air_temperature) + 0.01u"K")^4)
end
