"""
    Site(; latitude, longitude, elevation, slope, aspect, horizon_angles, sky_view_fraction)

The geography of the place being simulated. Ground albedo and roughness length are daily
forcings ([`DailyTimeseries`](@ref)).
"""
@kwdef struct Site{LAT,LON,EL,SL,AS,HA,SVF} <: AbstractSite
    latitude::LAT
    longitude::LON
    elevation::EL
    slope::SL
    aspect::AS
    horizon_angles::HA
    sky_view_fraction::SVF
end

"""
    example_site

Example site, Madison Wisconsin, USA.
"""
function example_site(;
    latitude = 43.07305u"°",
    longitude = -89.40123u"°",
    elevation = 226.0u"m",
    slope = 0.0u"°",
    aspect = 0.0u"°",
    horizon_angles = fill(0.0u"°", 24),
    sky_view_fraction = 1.0,
)
    Site(; latitude, longitude, elevation, slope, aspect, horizon_angles, sky_view_fraction)
end

"""
    SolarTerrain(site::Site, ground_albedo)

A `SolarRadiation.SolarTerrain` for `site` with `ground_albedo`, at the standard atmosphere's
pressure for the site's elevation.
"""
function SolarRadiation.SolarTerrain(site::Site, ground_albedo)
    SolarRadiation.SolarTerrain(;
        elevation = site.elevation,
        horizon_angles = site.horizon_angles,
        slope = site.slope,
        aspect = site.aspect,
        albedo = ground_albedo,
        atmospheric_pressure = atmospheric_pressure(site.elevation),
        latitude = site.latitude,
        longitude = site.longitude,
    )
end
