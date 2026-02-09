abstract type AbstractAtmosphericRadiationModel end
struct SwinbankAtmosphericRadiation <: AbstractAtmosphericRadiationModel end
struct CampbellNormanAtmosphericRadiation <: AbstractAtmosphericRadiationModel end

function atmospheric_radiation(::SwinbankAtmosphericRadiation, P_vap, tair)
    # Swinbank, Eq. 10.11 in Campbell and Norman 1998
    arad = uconvert(u"W*m^-2", ((9.2e-6 * (u"K"(tair))^2) * σ * (u"K"(tair))^4) / 1u"K^2")
    return P_vap, arad
end
function atmospheric_radiation(::CampbellNormanAtmosphericRadiation, P_vap, tair)
    # Campbell and Norman 1998 eq. 10.10 to get emissivity of sky
    arad = u"W/m^2"((1.72 * (ustrip(u"kPa", P_vap) / ustrip(u"K", tair + 0.01u"K"))^(1//7)) * σ * (u"K"(tair) + 0.01u"K")^4) 
    return P_vap, arad
end
function longwave_radiation(radiation_model=CampbellNormanAtmosphericRadiation();
    micro_terrain,
    environment_instant,
    surface_temperature,
)
    (; elevation, viewfactor) = micro_terrain
    (; P_atmos, reference_humidity, reference_temperature, surface_emissivity, cloud_emissivity, cloud_cover, shade) = environment_instant

    # Longwave radiation (handle both IR modes)
    wet_air_out = wet_air_properties(u"K"(reference_temperature), reference_humidity, P_atmos)

    # Atmospheric radiation
    P_vap, atmospheric_rad = atmospheric_radiation(radiation_model, wet_air_out.P_vap, reference_temperature)

    # Cloud radiation temperature (shade approximation, air temp - 2°C)
    cloud_radiation = σ * surface_emissivity * (u"K"(reference_temperature) - 2.0u"K")^4

    # Hillshade radiation temperature (approximated as air temperature)
    hillshade_radiation = σ * surface_emissivity * (u"K"(reference_temperature))^4

    # Ground surface radiation temperature
    surface_radiation = σ * cloud_emissivity * (u"K"(surface_temperature))^4

    # Clear sky fraction
    clear_sky_fraction = 1.0 - cloud_cover
    clear_component = atmospheric_rad * clear_sky_fraction
    cloudy_component = cloud_radiation * cloud_cover
    radiation_sky = (clear_component + cloudy_component) * (1.0 - shade)
    radiation_vegetation = shade * hillshade_radiation
    radiation_ground = (1 - shade) * surface_radiation + shade * hillshade_radiation
    radiation_hillshade = hillshade_radiation
    net_radiation = (radiation_sky + radiation_vegetation) * viewfactor + radiation_hillshade * (1.0 - viewfactor) - radiation_ground
    sky_temperature = (((radiation_sky + radiation_vegetation) * viewfactor + radiation_hillshade * (1.0 - viewfactor)) / σ)^(1//4)

    return (;
        Tsky=u"K"(sky_temperature),
        Qrad=net_radiation,
        Qrad_sky=radiation_sky,
        Qrad_veg=radiation_vegetation,
        Qrad_ground=radiation_ground,
        Qrad_hill=radiation_hillshade
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

Returns `(global_solar, diffuse_fraction)`; works with arrays but needs to not use 'similar' if to work with scalars.

Reference
Maxwell, E. L., "A Quasi-Physical Model for Converting Hourly
           Global Horizontal to Direct Normal Insolation", Technical
           Report No. SERI/TR-215-3087, Golden, CO: Solar Energy Research
           Institute, 1987.
"""
function cloud_adjust_radiation(output, cloud::AbstractArray, diffuse_clear_sky, direct_clear_sky, zenith::AbstractArray, doy;
    a=0.36, b=0.64, gamma=1.0,
)
    (; global_total, diffuse_total, direct_total) = output.solar_radiation
    global_radiation = global_total
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

    # 3) Split global into diffuse/direct using Erbs diffuse fraction vs clearness index
    ϵ = 1e-9u"W/m^2"
    clearness_index = global_radiation ./ max.(extraterrestrial_horizontal, ϵ)
    clearness_index = clamp.(clearness_index, 0.0, 1.2)
    diffuse_fraction = similar(clearness_index)
    for i in eachindex(clearness_index)
        if clearness_index[i] <= 0.22
            diffuse_fraction[i] = 1 - 0.09*clearness_index[i]
        elseif clearness_index[i] <= 0.80
            diffuse_fraction[i] = 0.9511 - 0.1604*clearness_index[i] + 4.388*clearness_index[i]^2 - 16.638*clearness_index[i]^3 + 12.336*clearness_index[i]^4
        else
            diffuse_fraction[i] = 0.165
        end
    end
    diffuse_fraction .= clamp.(diffuse_fraction, zero(eltype(diffuse_fraction)), oneunit(eltype(diffuse_fraction)))

    return (;
        global_solar=global_radiation,
        diffuse_fraction,
    )
end
