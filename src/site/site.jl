"""
    Site(; latitude, longitude, elevation, slope, aspect, horizon_angles,
          sky_view_fraction, albedo, roughness_height, atmospheric_pressure)

Properties of the place being simulated.

`atmospheric_pressure` here is the *baseline* pressure used by the solar
radiation model. Time-varying hourly pressure lives on the hourly forcing
(`HourlyTimeseries.pressure`).
"""
@kwdef struct Site{LAT,LON,EL,SL,AS,HA,SVF,AB,RH,AP} <: AbstractSite
    latitude::LAT
    longitude::LON
    elevation::EL
    slope::SL
    aspect::AS
    horizon_angles::HA
    sky_view_fraction::SVF
    albedo::AB
    roughness_height::RH
    atmospheric_pressure::AP
end

function example_site(;
    latitude = 43.07305u"°",
    longitude = -89.40123u"°",
    elevation = 226.0u"m",
    slope = 0.0u"°",
    aspect = 0.0u"°",
    horizon_angles = fill(0.0u"°", 24),
    sky_view_fraction = 1.0,
    albedo = 0.15,
    roughness_height = 0.004u"m",
    atmos_pressure = atmospheric_pressure(elevation),
)
    Site(;
        latitude, longitude, elevation, slope, aspect, horizon_angles,
        sky_view_fraction, albedo, roughness_height,
        atmospheric_pressure=atmos_pressure,
    )
end

"""
    SolarTerrain(site::Site)

Build a `SolarRadiation.SolarTerrain` from a `Site`. Used internally when
calling into `SolarRadiation.solar_radiation!`.
"""
function SolarRadiation.SolarTerrain(site::Site)
    SolarRadiation.SolarTerrain(;
        elevation = site.elevation,
        horizon_angles = site.horizon_angles,
        slope = site.slope,
        aspect = site.aspect,
        albedo = site.albedo,
        atmospheric_pressure = site.atmospheric_pressure,
        latitude = site.latitude,
        longitude = site.longitude,
    )
end
