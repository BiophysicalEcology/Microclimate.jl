"""
    KearneyLongwave()

NicheMapR's (Kearney et al.) surface-longwave budget. Combines the clear-sky
downwelling longwave from an `AbstractAtmosphericRadiationModel` with cloud
emissivity, sky view fraction, vegetation shade, and hillshade contributions
to give the net longwave at the surface.
"""
struct KearneyLongwave <: AbstractLongwaveScheme end

function precompute_longwave_sky(::KearneyLongwave, atmospheric_radiation_model;
    site,
    environment_instant,
    vapour_pressure_equation=GoffGratch(),
    # Optional explicit overrides so callers (e.g. snow surface) can pass
    # specific surface_emissivity / shade without rebuilding environment_instant.
    surface_emissivity=environment_instant.surface_emissivity,
    shade=environment_instant.shade,
)
    sky_view_fraction = site.sky_view_fraction
    (; atmospheric_pressure, reference_humidity, reference_temperature, cloud_emissivity, cloud_cover) = environment_instant

    wet_air_out = wet_air_properties(u"K"(reference_temperature), reference_humidity, atmospheric_pressure; vapour_pressure_equation)
    atmospheric_longwave = atmospheric_radiation(atmospheric_radiation_model, wet_air_out.vapour_pressure, reference_temperature)

    cloud_radiation = σ * cloud_emissivity * (u"K"(reference_temperature) - 2.0u"K")^4
    hillshade_radiation = σ * cloud_emissivity * (u"K"(reference_temperature))^4

    clear_sky_fraction = 1.0 - cloud_cover
    clear_component = atmospheric_longwave * clear_sky_fraction
    cloudy_component = cloud_radiation * cloud_cover
    longwave_radiation_sky = (clear_component + cloudy_component) * (1.0 - shade)
    longwave_radiation_vegetation = shade * hillshade_radiation
    longwave_radiation_hillshade = hillshade_radiation

    incoming_longwave = (longwave_radiation_sky + longwave_radiation_vegetation) * sky_view_fraction +
                        longwave_radiation_hillshade * (1.0 - sky_view_fraction)
    outgoing_coeff = (1.0 - shade) * σ * surface_emissivity
    ground_shade_term = shade * hillshade_radiation
    sky_temperature = (incoming_longwave / σ)^(1//4)

    return (;
        incoming_longwave,
        outgoing_coeff,
        ground_shade_term,
        sky_temperature=u"K"(sky_temperature),
        longwave_radiation_sky,
        longwave_radiation_vegetation,
        longwave_radiation_hillshade,
    )
end

function longwave_radiation(scheme::KearneyLongwave, atmospheric_radiation_model;
    site,
    environment_instant,
    surface_temperature,
    vapour_pressure_equation=GoffGratch(),
)
    sky = precompute_longwave_sky(scheme, atmospheric_radiation_model; site, environment_instant, vapour_pressure_equation)
    (; incoming_longwave, ground_shade_term) = sky

    surface_radiation = σ * environment_instant.surface_emissivity * (u"K"(surface_temperature))^4
    shade = environment_instant.shade
    longwave_radiation_ground = (1 - shade) * surface_radiation + ground_shade_term
    net_longwave_radiation = incoming_longwave - (1 - shade) * surface_radiation - ground_shade_term

    return (;
        sky.sky_temperature,
        sky.longwave_radiation_sky,
        sky.longwave_radiation_vegetation,
        sky.longwave_radiation_hillshade,
        longwave_radiation_ground,
        net_longwave_radiation,
    )
end
