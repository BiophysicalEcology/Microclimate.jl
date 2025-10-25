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
    # TODO these are not the real names
    (; elevation, P_atmos, viewfactor) = micro_terrain
    (; reference_humidity, reference_temperature, surface_emissivity, cloud_emissivity, cloud_cover, shade) = environment_instant

    # Short names, hardly worth it
    tsurf = surface_temperature
    tair = reference_temperature
    rh = reference_humidity
    slep = surface_emissivity
    sle = cloud_emissivity
    cloud = cloud_cover

    # Longwave radiation (handle both IR modes)
    wet_air_out = wet_air_properties(u"K"(tair); rh, P_atmos)

    # Atmospheric radiation
    P_vap, arad = atmospheric_radiation(radiation_model, wet_air_out.P_vap, tair)

    # Cloud radiation temperature (shade approximation, TAIR - 2°C)
    crad = σ * slep * (u"K"(tair) - 2.0u"K")^4

    # Hillshade radiation temperature (approximated as air temperature)
    hrad = σ * slep * (u"K"(tair))^4

    # Ground surface radiation temperature
    srad = σ * sle * (u"K"(tsurf))^4

    # Clear sky fraction
    clr = 1.0 - cloud / 100.0
    clear = arad * clr
    clod = crad * (cloud / 100.0)
    qradsk = (clear + clod) * ((100 - shade) / 100.0)
    qradvg = (shade / 100.0) * hrad
    qradgr = ((100.0 - shade) / 100.0) * srad + (shade / 100.0) * hrad
    qradhl = hrad
    qrad = (qradsk + qradvg) * viewfactor + qradhl * (1.0 - viewfactor) - qradgr
    tsky = (((qradsk + qradvg) * viewfactor + qradhl * (1.0 - viewfactor)) / σ)^(1//4)

    return (;
        # TODO standardise these names with their target uses
        Tsky=u"K"(tsky),
        Qrad=qrad, # e.g. this is Q_infrared in `solar_radiation`
        Qrad_sky=qradsk,
        Qrad_veg=qradvg,
        Qrad_ground=qradgr,
        Qrad_hill=qradhl
    )
end

"""
    cloud_adjust_radiation(cloud, D_cs, B_cs, zenith, doy; a=0.25, b=0.5, gamma=1.0)

Compute global (G), diffuse (D), and direct-beam (B) solar on a horizontal surface
given cloud cover fraction `cloud` (0–1), clear-sky diffuse `D_cs` and direct `B_cs`,
solar zenith angle `zenith` (radians), and day-of-year `doy`.

- Ångström scaling: G = (a + b*S) * (D_cs + B_cs), with S ≈ (1 - cloud)^gamma
- Diffuse fraction via Erbs (uses extraterrestrial horizontal irradiance) via
    a clearness index (Maxwwell 1987) which is the ratio of global to extraterrestrial
    irradiance on a horizontal plane

Returns `(global_solar, diffuse_fraction)`; works with arrays but needs to not use 'similar' if to work with
    scalars.

Reference
Maxwell, E. L., "A Quasi-Physical Model for Converting Hourly
           Global Horizontal to Direct Normal Insolation", Technical
           Report No. SERI/TR-215-3087, Golden, CO: Solar Energy Research
           Institute, 1987.
"""
function cloud_adjust_radiation(output, cloud::AbstractArray, D_cs, B_cs, zenith::AbstractArray, doy; 
    a=0.36, b=0.64, gamma=1.0,
)
    (; global_total, diffuse_total, direct_total) = output.solar_radiation
    G, D, B = (global_total, diffuse_total, direct_total)
    # Solar geometry
    cosz     = cos.(zenith)
    cosz_pos = max.(cosz, 0.0)

    # 1) Extraterrestrial horizontal irradiance (W/m²)
    I_sc  = 1367.0u"W/m^2"
    E0    = 1.00011 .+ 0.034221*cosd.(360.0 .* (doy .- 1) ./ 365.0) .+
                     0.00128*sind.(360.0 .* (doy .- 1) ./ 365.0) .+
                     0.000719*cosd.(2 .* 360.0 .* (doy .- 1) ./ 365.0) .+
                     0.000077*sind.(2 .* 360.0 .* (doy .- 1) ./ 365.0)
    G0h   = I_sc .* E0 #.* cosz_pos  # on horizontal; 0 at night

    # 2) Ångström–Prescott scaling of clear-sky global by cloud cover
    S     = (1 .- cloud).^gamma                 # approx. sunshine fraction
    T     = a .+ b .* S                         # transmittance
    G_cs  = D_cs .+ B_cs                        # clear-sky global
    G     .= max.(T .* G_cs, 0.0u"W/m^2") #.* (cosz_pos .> 0)  # zero at night

    # 3) Split G into diffuse/direct using Erbs diffuse fraction vs clearness index K_t
    ϵ     = 1e-9u"W/m^2"
    Kt    = G ./ max.(G0h, ϵ)
    Kt    = clamp.(Kt, 0.0, 1.2)
    Fd = similar(Kt) # diffuse fraction
    for i in eachindex(Kt)
        if Kt[i] <= 0.22
            Fd[i] = 1 - 0.09*Kt[i]
        elseif Kt[i] <= 0.80
            Fd[i] = 0.9511 - 0.1604*Kt[i] + 4.388*Kt[i]^2 - 16.638*Kt[i]^3 + 12.336*Kt[i]^4
        else
            Fd[i] = 0.165
        end
    end
    # TODO probably this still allocates because of aliasing 
    Fd .= clamp.(Fd, zero(eltype(Fd)), oneunit(eltype(Fd)))

    #D .= Fd .* G
    #B .= G .- D

    # Zero everything at night
    # night = (cosz_pos .== 0)
    # D[night] .= 0.0u"W/m^2"
    # B[night] .= 0.0u"W/m^2"
    # G[night] .= 0.0u"W/m^2"

    return (;
        global_solar=G,
        diffuse_fraction=Fd,
    )
end