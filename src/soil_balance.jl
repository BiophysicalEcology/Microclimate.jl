function soil_energy_balance!(dT, T, i::MicroInput, t)
    
    # extract input
    p = i.params
    f = i.forcing

    # extract parameters
    ruf = p.ruf
    pctwet = p.pctwet
    sle = p.sle
    slep = p.slep
    refl = p.refl
    viewf = p.viewf
    elev = p.elev
    sabnew = 1.0 - refl
    slope = p.slope
    shade = p.shade
    dep = p.dep
    refhyt = p.refhyt
    d0 = p.d0
    zh = p.zh
    tdeep = p.tdeep
    nodes = p.nodes
    soilprops = p.soilprops

    θ_soil = p.θ_soil # parameter for now

    N = length(dep)
    #dT = fill(0.0u"K/minute", N)
    #dT .= (0.0u"K/minute")

    # get soil properties and convert to cal/cm/g/C
    λ_b, cp_b, ρ_b = soil_properties(T, θ_soil, nodes, soilprops, elev)
    λ_b = u"cal/cm/K/minute".(λ_b)
    cp_b = u"cal/g/K".(cp_b)
    ρ_b = u"g/cm^3".(ρ_b)

    # Get environmental data at time t
    tair = f.TAIRt(ustrip(t))
    vel = f.VELt(ustrip(t))
    zenr = min(90u"°", u"°"(round(f.ZENRt(ustrip(t)), digits=3)))
    solr = u"cal/cm^2/minute"(max(0.0u"W/m^2", f.SOLRt(ustrip(t))))
    cloud = f.CLDt(ustrip(t))
    rh = f.RHt(ustrip(t))
    zslr = f.ZSLt(ustrip(t))

    T[N] = tdeep # set boundary condition of deep soil temperature

    depp = fill(0.0u"cm", N + 1)
    depp[1:N] = dep
    # Compute soil layer properties
    wc = fill(1.0u"cal/K/cm^2", N)
    c = fill(1.0u"cal/K/cm^2/minute", N)
    for i in 1:N
        rcsp = ρ_b[i] * cp_b[i]
        if i == 1
            wc[i] = rcsp * depp[2] / 2.0
            sok = λ_b[1]
            c[i] = sok / depp[2]
        else
            wc[i] = rcsp * (depp[i+1] - depp[i-1]) / 2.0
            sok = λ_b[i]
            c[i] = sok / (depp[i+1] - depp[i])
        end
    end

    # Solar radiation
    if cloud > 0.0
        solr = solr * (0.36 + 0.64 * (1.0-(cloud / 100.0))) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
    end
    qsolar = sabnew * solr * ((100.0 - shade) / 100.0)
    if slope > 0 && zenr < 90u"°"
        cz = cosd(zenr)
        czsl = cosd(zslr)
        qsolar = (qsolar / cz) * czsl
    end

    # Longwave radiation (handle both IR modes)
    # Constants
    #σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ) # Stefan-Boltzmann constant, W/m^2/K^4
    σ = Unitful.uconvert(u"cal/cm^2/K^4/minute", Unitful.σ) # Stefan-Boltzmann constant
    P_atmos = get_pressure(elev)
    wet_air_out = wet_air(u"K"(tair); rh=rh, P_atmos=P_atmos)

    # Atmospheric radiation
    P_vap = wet_air_out.P_vap
    arad = u"cal/cm^2/minute"(ustrip(1.72 * (ustrip(u"kPa"(P_vap))/ustrip(u"K"(tair))) ^ (1.0/7.0)) * Unitful.uconvert(u"W/m^2/K^4", Unitful.σ) * (u"K"(tair)) ^ 4) # Campbell and Norman 1998 eq. 10.10 to get emissivity of sky
    #arad = ((0.0000092 * (u"K"(tair))^2) * σ * (u"K"(tair))^4) / 1u"K^2" # Swinbank, Eq. 10.11 in Campbell and Norman 1998

    # Cloud radiation temperature (shade approximation, TAIR - 2°C)
    crad = σ * slep * (u"K"(tair) - 2u"K")^4

    # Hillshade radiation temperature (approximated as air temperature)
    hrad = σ * slep * (u"K"(tair))^4

    # Ground surface radiation temperature
    srad = σ * sle * (u"K"(T[1]))^4

    # Clear sky fraction
    clr = 1.0 - cloud / 100.0
    clod = crad * (cloud / 100)
    qradsk = (arad * clr + clod) * ((100 - shade) / 100.0)
    qradvg = (shade / 100.0) * hrad
    qradgr = ((100.0 - shade) / 100.0) * srad + (shade / 100.0) * hrad
    qradhl = hrad
    qrad = (qradsk + qradvg) * viewf + qradhl * (1.0 - viewf) - qradgr
    # Conduction
    qcond = c[1] * (T[2] - T[1])

    # Convection
    profile_out = get_profile(
        refhyt = refhyt,
        ruf = ruf, 
        d0 = d0, 
        zh = zh, 
        D0cm=u"°C"(T[1]), 
        TAREF=u"°C"(tair), 
        VREF=vel, 
        ZEN=zenr, 
        heights=[0.01] .* u"m", 
        rh=rh, 
        elev=elev
        )
    qconv = profile_out.QCONV
    hc = max(abs(qconv / (T[1] - tair)), 0.5u"W/m^2/K")
    qconv = u"cal/cm^2/minute"(qconv) # now to cal/cm/g/C units

    # Evaporation
    wet_air_out = wet_air(u"K"(tair); rh=rh, P_atmos=P_atmos)
    cp_air = wet_air_out.cp
    ρ_air = wet_air_out.ρ_air
    hd = (hc / (cp_air * ρ_air)) * (0.71 / 0.60)^0.666
    qevap, gwsurf = evap(tsurf=u"K"(T[1]), tair=u"K"(tair), rh=rh, rhsurf=100.0, hd=hd, elev=elev, pctwet=pctwet, sat=false)
    qevap = u"cal/cm^2/minute"(qevap) # now to cal/cm/g/C units

    # Energy balance at surface
    dT[1] = (qsolar + qrad + qcond + qconv - qevap) / wc[1]

    # Soil conduction for internal nodes
    for i in 2:N-1
        dT[i] = (c[i-1] * (T[i-1] - T[i]) + c[i] * (T[i+1] - T[i])) / wc[i]
    end

    # Lower boundary condition
    dT[N] = 0.0u"K/minute"  # or set T[N] = T_surface from data
end

function evap(;tsurf, tair, rh, rhsurf, hd, elev, pctwet, sat)
    # Assumes all units are SI (Kelvin, Pascal, meters, seconds, kg, etc.)

    # Ground-level variables, shared via global or passed in as needed
    #global shayd, altt, maxshd, sabnew, ptwet, rainfall

    # Output variables
    qevap = 0.0
    gwsurf = 0.0

    tsurf = tsurf < u"K"(-81.0u"°C") ? u"K"(-81.0u"°C") : tsurf

    # Atmospheric pressure from elevation
    P_atmos = get_pressure(elev)

    # surface and air vapor densities
    ρ_vap_surf = wet_air(u"K"(tsurf); rh=rhsurf, P_atmos=P_atmos).ρ_vap
    ρ_vap_air = wet_air(u"K"(tair); rh=rh, P_atmos=P_atmos).ρ_vap

    # Effective wet surface fraction
    effsur = sat ? 1.0 : pctwet / 100.0

    # Water evaporated from surface (kg/s/m^2)
    water = effsur * hd * (ρ_vap_surf - ρ_vap_air)

    # Latent heat of vaporization (J/kg)
    htovpr = tsurf > u"K"(0.0u"°C") ?
        (2500.8 - 2.36 * ustrip(u"°C"(tsurf)) + 0.0016 * ustrip(u"°C"(tsurf))^2 - 0.00006 * ustrip(u"°C"(tsurf))^3) :
        (2834.1 - 0.29 * ustrip(u"°C"(tsurf)) - 0.004 * ustrip(u"°C"(tsurf))^2)
    htovpr = (htovpr * 1000.0 )u"J/kg" # convert to J/kg

    # Energy flux due to evaporation (W/m² or converted)
    qevap = water * htovpr  # SI units for water budget calcs

    # Mass flux (g/s)
    gwsurf = u"g/s/m^2"(water)

    # No water loss if TSURF ≤ 0 (e.g., melting snow only)
    if tsurf <= u"K"(0.0u"°C")
        gwsurf = 0.0u"g/s/m^2"
    end

    return qevap, gwsurf
end