"""
    precompute_longwave_sky([radiation_model]; site, environment_instant, vapour_pressure_equation)

Compute the per-hour sky/atmospheric longwave radiation terms that do not depend on surface
temperature. Returns a named tuple that can be cached for the duration of a timestep and passed
into the ODE solver via `SoilEnergyInputs`, avoiding repeated calls to `wet_air_properties`
on every internal ODE step.

The `net_longwave_radiation` at any surface temperature `T` can then be recovered cheaply as:
    net_Q = incoming_longwave - outgoing_coeff * T^4 - ground_shade_term
"""
function precompute_longwave_sky(radiation_model=CampbellNormanAtmosphericRadiation();
    site,
    environment_instant,
    vapour_pressure_equation=GoffGratch(),
)
    sky_view_fraction = site.sky_view_fraction
    (; atmospheric_pressure, reference_humidity, reference_temperature, surface_emissivity, cloud_emissivity, cloud_cover, shade) = environment_instant

    wet_air_out = wet_air_properties(u"K"(reference_temperature), reference_humidity, atmospheric_pressure; vapour_pressure_equation)
    _, atmospheric_longwave = atmospheric_radiation(radiation_model, wet_air_out.vapour_pressure, reference_temperature)

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
        radiation_model,
    )
end

function longwave_radiation(radiation_model=CampbellNormanAtmosphericRadiation();
    site,
    environment_instant,
    surface_temperature,
    vapour_pressure_equation=GoffGratch(),
)
    sky = precompute_longwave_sky(radiation_model; site, environment_instant, vapour_pressure_equation)
    (; incoming_longwave, ground_shade_term) = sky

    surface_radiation = σ * environment_instant.surface_emissivity * (u"K"(surface_temperature))^4
    shade = environment_instant.shade
    longwave_radiation_ground = (1 - shade) * surface_radiation + ground_shade_term
    net_longwave_radiation = incoming_longwave - (1 - shade) * surface_radiation - ground_shade_term

    return (;
        sky.sky_temperature,
        net_longwave_radiation,
        sky.longwave_radiation_sky,
        sky.longwave_radiation_vegetation,
        longwave_radiation_ground,
        sky.longwave_radiation_hillshade,
    )
end

"""
    cloud_adjust_radiation(output, cloud, diffuse_clear_sky, direct_clear_sky, zenith, doy; a=0.36, b=0.64, gamma=1.0)

Compute global, diffuse, and direct-beam solar radiation on a horizontal surface
given cloud cover fraction `cloud` (0–1), clear-sky diffuse `diffuse_clear_sky` and
direct `direct_clear_sky`, solar zenith angle `zenith` (radians), and day-of-year `doy`.

- Ångström scaling: global = (a + b*sunshine_fraction) * (diffuse_clear_sky + direct_clear_sky),
  with sunshine_fraction ≈ (1 - cloud)^gamma
- Diffuse fraction via Erbs (uses extraterrestrial horizontal irradiance) via
  a clearness index (Maxwell 1987) which is the ratio of global to extraterrestrial
  irradiance on a horizontal plane

Returns `(global_radiation, diffuse_fraction)`; works with arrays but needs to not use 'similar' if to work with scalars.

Reference
Maxwell, E. L., "A Quasi-Physical Model for Converting Hourly
           Global Horizontal to Direct Normal Insolation", Technical
           Report No. SERI/TR-215-3087, Golden, CO: Solar Energy Research
           Institute, 1987.
"""
function cloud_adjust_radiation(output, cloud::AbstractArray, diffuse_clear_sky, direct_clear_sky, zenith::AbstractArray, doy;
    diffuse_fraction_model::AbstractDiffuseFractionModel=ErbsDiffuseFraction(),
    cloud_adjust_model::AbstractCloudAdjustModel=Angstrom(),
)
    (; global_horizontal) = output.solar_radiation
    global_radiation = global_horizontal
    diffuse_fraction = output.diffuse_fraction
    solar_constant = 1367.0u"W/m^2"
    ϵ = 1e-9u"W/m^2"
    a, b = cloud_adjust_model.a, cloud_adjust_model.b
    @inbounds for i in eachindex(global_radiation)
        # 1) Extraterrestrial horizontal irradiance for this hour's day
        d = doy[i]
        θ1 = 360.0 * (d - 1) / 365.0
        θ2 = 2.0 * θ1
        ec = 1.00011 + 0.034221*cosd(θ1) + 0.00128*sind(θ1) +
                       0.000719*cosd(θ2) + 0.000077*sind(θ2)
        eth = solar_constant * ec
        # 2) Ångström–Prescott scaling
        sf = sunshine_fraction(cloud_adjust_model, cloud[i])
        t  = a + b * sf
        gcs = diffuse_clear_sky[i] + direct_clear_sky[i]
        gr = max(t * gcs, 0.0u"W/m^2")
        global_radiation[i] = gr
        # 3) Split global into diffuse/direct via clearness index
        ci = clamp(gr / max(eth, ϵ), 0.0, 1.2)
        diffuse_fraction[i] = clamp(calc_diffuse_fraction(diffuse_fraction_model, ci), 0.0, 1.0)
    end
    return (; global_radiation, diffuse_fraction)
end
