abstract type AbstractAtmosphericRadiationModel end
struct SwinbankAtmosphericRadiation <: AbstractAtmosphericRadiationModel end
struct CampbellNormanAtmosphericRadiation <: AbstractAtmosphericRadiationModel end

function atmospheric_radiation(::SwinbankAtmosphericRadiation, vapour_pressure, air_temperature)
    # Swinbank, Eq. 10.11 in Campbell and Norman 1998
    atmospheric_longwave = uconvert(u"W*m^-2", ((9.2e-6 * (u"K"(air_temperature))^2) * σ * (u"K"(air_temperature))^4) / 1u"K^2")
    return vapour_pressure, atmospheric_longwave
end
function atmospheric_radiation(::CampbellNormanAtmosphericRadiation, vapour_pressure, air_temperature)
    # Campbell and Norman 1998 eq. 10.10 to get emissivity of sky
    atmospheric_longwave = u"W/m^2"((1.72 * (ustrip(u"kPa", vapour_pressure) / ustrip(u"K", air_temperature + 0.01u"K"))^(1//7)) * σ * (u"K"(air_temperature) + 0.01u"K")^4)
    return vapour_pressure, atmospheric_longwave
end
"""
    precompute_longwave_sky([radiation_model]; micro_terrain, environment_instant, vapour_pressure_equation)

Compute the per-hour sky/atmospheric longwave radiation terms that do not depend on surface
temperature. Returns a named tuple that can be cached for the duration of a timestep and passed
into the ODE solver via `SoilEnergyInputs`, avoiding repeated calls to `wet_air_properties`
on every internal ODE step.

The `net_longwave_radiation` at any surface temperature `T` can then be recovered cheaply as:
    net_Q = incoming_longwave - outgoing_coeff * T^4 - ground_shade_term
"""
function precompute_longwave_sky(radiation_model=CampbellNormanAtmosphericRadiation();
    micro_terrain,
    environment_instant,
    vapour_pressure_equation=GoffGratch(),
)
    (; viewfactor) = micro_terrain
    (; atmospheric_pressure, reference_humidity, reference_temperature, surface_emissivity, cloud_emissivity, cloud_cover, shade) = environment_instant

    wet_air_out = wet_air_properties(u"K"(reference_temperature), reference_humidity, atmospheric_pressure; vapour_pressure_equation)
    _, atmospheric_longwave = atmospheric_radiation(radiation_model, wet_air_out.vapour_pressure, reference_temperature)

    cloud_radiation = σ * surface_emissivity * (u"K"(reference_temperature) - 2.0u"K")^4
    hillshade_radiation = σ * surface_emissivity * (u"K"(reference_temperature))^4

    clear_sky_fraction = 1.0 - cloud_cover
    clear_component = atmospheric_longwave * clear_sky_fraction
    cloudy_component = cloud_radiation * cloud_cover
    longwave_radiation_sky = (clear_component + cloudy_component) * (1.0 - shade)
    longwave_radiation_vegetation = shade * hillshade_radiation
    longwave_radiation_hillshade = hillshade_radiation

    incoming_longwave = (longwave_radiation_sky + longwave_radiation_vegetation) * viewfactor +
                        longwave_radiation_hillshade * (1.0 - viewfactor)
    outgoing_coeff = (1.0 - shade) * σ * cloud_emissivity
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

function longwave_radiation(radiation_model=CampbellNormanAtmosphericRadiation();
    micro_terrain,
    environment_instant,
    surface_temperature,
    vapour_pressure_equation=GoffGratch(),
)
    sky = precompute_longwave_sky(radiation_model; micro_terrain, environment_instant, vapour_pressure_equation)
    (; incoming_longwave, ground_shade_term) = sky

    surface_radiation = σ * environment_instant.cloud_emissivity * (u"K"(surface_temperature))^4
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
    a=0.36, b=0.64, gamma=1.0,
)
    (; global_horizontal, diffuse_horizontal, direct_horizontal) = output.solar_radiation
    global_radiation = global_horizontal
    # Solar geometry
    cos_zenith = cos.(zenith)
    cos_zenith_positive = max.(cos_zenith, 0.0)

    # 1) Extraterrestrial horizontal irradiance (W/m²)
    solar_constant = 1367.0u"W/m^2"
    eccentricity_correction = 1.00011 .+ 0.034221*cosd.(360.0 .* (doy .- 1) ./ 365.0) .+
                     0.00128*sind.(360.0 .* (doy .- 1) ./ 365.0) .+
                     0.000719*cosd.(2 .* 360.0 .* (doy .- 1) ./ 365.0) .+
                     0.000077*sind.(2 .* 360.0 .* (doy .- 1) ./ 365.0)
    extraterrestrial_horizontal = solar_constant .* eccentricity_correction

    # 2) Ångström–Prescott scaling of clear-sky global by cloud cover
    sunshine_fraction = (1 .- cloud).^gamma
    transmittance = a .+ b .* sunshine_fraction
    global_clear_sky = diffuse_clear_sky .+ direct_clear_sky
    global_radiation .= max.(transmittance .* global_clear_sky, 0.0u"W/m^2")

    # 3) Split global into diffuse/direct using clearness index
    ϵ = 1e-9u"W/m^2"
    clearness_index = global_radiation ./ max.(extraterrestrial_horizontal, ϵ)
    clearness_index = clamp.(clearness_index, 0.0, 1.2)
    diffuse_fraction = similar(clearness_index)
    for i in eachindex(clearness_index)
        diffuse_fraction[i] = calc_diffuse_fraction(diffuse_fraction_model, clearness_index[i])
    end
    diffuse_fraction .= clamp.(diffuse_fraction, zero(eltype(diffuse_fraction)), oneunit(eltype(diffuse_fraction)))

    return (;
        global_radiation,
        diffuse_fraction,
    )
end
