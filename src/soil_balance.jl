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

    θ_soil = p.θ_soil

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
        # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
        solr = solr * (0.36 + 0.64 * (1.0-(cloud / 100.0))) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
    end
    qsolar = sabnew * solr * ((100.0 - shade) / 100.0)
    if slope > 0 && zenr < 90u"°"
        cz = cosd(zenr)
        czsl = cosd(zslr)
        qsolar = (qsolar / cz) * czsl
    end

    # Longwave radiation
    longwave_out = get_longwave(
        elev = elev, 
        rh = rh, 
        tair = tair, 
        tsurf = T[1], 
        slep = slep, 
        sle = sle, 
        cloud = cloud, 
        viewf = viewf
        )
    qrad = u"cal/cm^2/minute"(longwave_out.Qrad)

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
    λ_evap = get_λ_evap(tsurf) 

    # Energy flux due to evaporation (W/m² or converted)
    qevap = water * λ_evap  # SI units for water budget calcs

    # Mass flux (g/s)
    gwsurf = u"g/s/m^2"(water)

    # No water loss if TSURF ≤ 0 (e.g., melting snow only)
    if tsurf <= u"K"(0.0u"°C")
        gwsurf = 0.0u"g/s/m^2"
    end

    return qevap, gwsurf
end

function soil_water_balance(;
    PE = fill(1.1, 19)u"J/kg", # Air entry potential (J/kg) (19 values descending through soil for specified soil nodes in parameter DEP and points half way between)
    KS = fill(0.0037, 19)u"kg*s/m^3", # Saturated conductivity, (kg s/m3) (19 values descending through soil for specified soil nodes in parameter DEP and points half way between)
    BB = fill(4.5, 19), # Campbell's soil 'b' parameter (-) (19 values descending through soil for specified soil nodes in parameter DEP and points half way between)
    BD = fill(1.3, 19)u"Mg/m^3", # Soil bulk density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter DEP and points half way between)
    DD = fill(2.56, 19)u"Mg/m^3", # Soil density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter DEP and points half way between)
    rh_loc = 20.0,
    θ_soil = fill(0.2, 18),
    ET = 1.3e-5u"kg/m^2/s",
    T10 = fill(293.15u"K", 10),
    depth = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm",
    dt = 360u"s",
    elev = 0.0u"m",
    L = [0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0]*10000u"m/m^3", # root density, m m-3
    rw = 2.5E+10u"m^3/kg/s", # resistance per unit length of root, m3 kg-1 s-1
    pc = -1500.0u"J/kg", # critical leaf water potential for stomatal closure, J kg-1
    rl = 2000000.0u"m^4/kg/s", # leaf resistance, m4 kg-1 s-1
    sp = 10.0, # stability parameter, -
    r1 = 0.001u"m", # root radius, m
    lai = 0.1,
    im = 1e-6u"kg/m^2/s", # maximum overall mass balance error allowed, kg m-2 s-1
    maxcount=500
)

    A = zeros(Float64, 19)
    B = zeros(Float64, 19)
    C = zeros(Float64, 19)
    F = zeros(Float64, 19)
    P = zeros(Float64, 19)*u"J/kg"       # matric potential J/kg
    Z = zeros(Float64, 19)u"m"                # depth nodes
    V = zeros(Float64, 19)u"kg/m^2"
    DP = zeros(Float64, 19)
    W = zeros(Float64, 19)*u"m^3/m^3"    # water content m3/m3
    WN = zeros(Float64, 19)*u"m^3/m^3"   # water content m3/m3
    K = zeros(Float64, 19)*u"kg*s/m^3"   # hydraulic conductivity, kg s/m3
    CP = zeros(Float64, 19)
    H = zeros(Float64, 19)
    JV = zeros(Float64, 19)
    DJ = zeros(Float64, 19)
    T = zeros(Float64, 19)u"K"
    humid = zeros(Float64, 18)
    potent = zeros(Float64, 18)u"J/kg"
    rootpot = zeros(Float64, 18)u"J/kg"
    PR = zeros(Float64, 19)u"J/kg"
    PP = zeros(Float64, 19)u"J/kg"
    B1 = zeros(Float64, 19)
    N = zeros(Float64, 19)
    N1 = zeros(Float64, 19)
    WS = zeros(Float64, 19)
    RR = zeros(Float64, 19)
    E = zeros(Float64, 19)
    RS = zeros(Float64, 19)
    BZ = zeros(Float64, 19)

    M = 18 #number of elements
    P_atmos = get_pressure(elev)

    # Constants
    MW = 0.01801528u"kg/mol" # molar mass of water
    WD = 1000.0u"kg/m^3"     # kg/m³
    DV = 0.000024u"m^2/s"    # m²/s

    # Convert PE to negative absolute value
    PE .= -abs.(PE)

    # Initialize PP from PE
    PP .= PE
    PP .= -abs.(PP)

    # Saturation water content
    WS = 1.0 .- BD ./ DD  # WS = 1 - BD/DD

    # Depth to lower boundary (m)
    Z[M+1] = u"m"(depth[10])


    # Soil hydraulic properties
    B1 .= 1.0 ./ BB
    N .= 2.0 .+ 3.0 ./ BB
    N1 .= 1.0 .- N

    # Fill Z using provided depth vector
    j = 2
    for i in 3:18
        if isodd(i)
            Z[i] = depth[j]
            j += 1
        else
            Z[i] = Z[i-1] + (depth[j] - Z[i-1]) / 2
        end
    end

    # Interpolate T from temp
    j = 1
    for i in 1:19
        if isodd(i)
            T[i] = T10[j]
            j += 1
        else
            T[i] = T[i-1] + (T10[j] - T[i-1]) / 2
        end
    end

    # Set Z[1] and Z[2] to 0 m
    Z[1] = 0.0u"m"
    Z[2] = 0.0u"m"

    # Set initial water content and related variables
    for i in 2:M
        WN[i] = θ_soil[i-1]
        P[i] = PE[i] * (WS[i] / WN[i])^BB[i]
        H[i] = exp(MW * P[i] / (R * T[i-1]))
        K[i] = KS[i] * (PE[i] / P[i])^N[i]
        W[i] = WN[i]
    end

    # Bulk water mass per soil layer
    for i in 2:M
        V[i] = WD * (Z[i+1] - Z[i-1]) / 2
    end
    # Lower boundary condition
    P[M+1] = PE[M] * (WS[M+1] / WS[M+1])^BB[M]
    H[M+1] = 1.0
    W[M+1] = WS[M+1]
    WN[M+1] = WS[M+1]
    Z[1] = -1e10u"m"
    Z[M+1] = 1e20u"m"
    K[M+1] = KS[M] * (PE[M] / P[M+1])^N[M+1]

    # Initialize root water uptake variables
    RR = zeros(M + 1)u"m^4/kg/s"
    BZ = zeros(M + 1)u"m"
    for i in 2:M
        if L[i] > 0.0u"m/m^3"
            RR[i] = rw / (L[i] * (Z[i+1] - Z[i-1]) / 2.0)
            BZ[i] = (1.0 - M) * log(π * r1^2 * L[i]) / (4.0 * π * L[i] * (Z[i+1] - Z[i-1]) / 2.0)
        else
            RR[i] = 1e20u"m^4/kg/s"
            BZ[i] = 0.0u"m"
        end
    end

    P[1] = P[2]
    K[1] = 0.0u"kg*s/m^3"

    # Evapotranspiration
    EP = exp(-0.82 * lai) * ET
    TP = ET - EP

    # Plant water uptake
    PB1 = 0.0u"J*s/m^4" 
    RB1 = 0.0u"kg*s/m^4"
    PL = 0.0u"J/kg"
    RS = zeros(M + 1)u"m^4/kg/s"
    for i in 2:M
        RS[i] = BZ[i] / K[i]
        PB1 += P[i] / (RS[i] + RR[i])
        RB1 += 1.0 / (RS[i] + RR[i])
    end
    PB = PB1 / RB1
    RB = (1.0 / RB1)

    # Newton-Raphson to estimate PL
    count = 0
    while count < maxcount
        if PL > PB
            PL = PB - TP * (RB + rl) 
        end
        XP = (PL / pc)^sp
        SL = TP * (RB + rl) * sp * XP / (PL * (1.0 + XP)^2) - 1.0
        FF = PB - PL - TP * (RB + rl) / (1.0 + XP)
        PL -= FF / SL
        count += 1
        if abs(FF) <= 10.0u"J/kg"
            break
        end
    end
    PL = u"J/kg"(PL) # keep units in J/kg = m^2/s^2
    TR = TP / (1.0 + XP)
    E = zeros(M + 1)u"kg/m^2/s"
    for i in 2:M
        E[i] = (P[i] - PL - rl * TR) / (RR[i] + RS[i])
    end

    # Convergence loop
    SE = 0.0
    count = 0
    JV = zeros(M + 1)u"kg/m^2/s"
    DJ = zeros(M + 1)u"kg*s/m^4"
    CP = zeros(M + 1)u"kg*s/m^4"
    A = zeros(M + 1)u"kg*s/m^4"
    B = zeros(M + 1)u"kg*s/m^4"
    C = zeros(M + 1)u"kg*s/m^4"
    C2 = zeros(M + 1) # adding this to deal with ratio check
    F = zeros(M + 1)u"kg/m^2/s"
    F2 = zeros(M + 1)u"m^2/s^2"
    DP = zeros(M + 1)u"J/kg"
    while count < maxcount
        SE = 0.0u"kg/m^2/s"
        count += 1
        for i in 2:M
            K[i] = KS[i] * (PE[i] / P[i])^N[i]
        end

        JV[1] = EP * (H[2] - rh_loc) / (1.0 - rh_loc)
        DJ[1] = EP * MW * H[2] / (Unitful.R * T[1] * (1.0 - rh_loc))

        for i in 2:M
            VP = wet_air(u"K"(T[i]); rh=100.0, P_atmos=P_atmos).ρ_vap
            KV = 0.66 * DV * VP * (WS[i] - (WN[i] + WN[i+1]) / 2.0) / (Z[i+1] - Z[i])
            JV[i] = KV * (H[i+1] - H[i])
            DJ[i] = MW * H[i] * KV / (R * T[i])
            CP[i] = -1.0 * V[i] * WN[i] / (BB[i] * P[i] * dt)
            A[i] = -1.0 * K[i-1] / (Z[i] - Z[i-1]) + Unitful.gn * N[i] * K[i-1] / P[i-1]
            C[i] = -1.0 * K[i+1] / (Z[i+1] - Z[i])
            B[i] = K[i] / (Z[i] - Z[i-1]) + K[i] / (Z[i+1] - Z[i]) + CP[i] - Unitful.gn * N[i] * K[i] / P[i] + DJ[i-1] + DJ[i]
            F[i] = ((P[i] * K[i] - P[i-1] * K[i-1]) / (Z[i] - Z[i-1]) - (P[i+1] * K[i+1] - P[i] * K[i]) / (Z[i+1] - Z[i])) / N1[i] + V[i] * (WN[i] - W[i]) / dt - Unitful.gn * (K[i-1] - K[i]) + JV[i-1] - JV[i] + E[i]
            SE += abs(F[i])
        end

        for i in 2:M-1
            C2[i] = C[i] / B[i]
            C[i] = C2[i] < 1e-8 ? 0.0u"kg*s/m^4" : C[i]
            F2[i] = F[i] / B[i]
            B[i+1] -= A[i+1] * C2[i]
            F[i+1] -= A[i+1] * F2[i]
        end

        DP[M] = F[M] / B[M]
        P[M] -= DP[M]
        P[M] = min(P[M], PE[M])

        for i in (M-1):-1:2
            DP[i] = F2[i] - C2[i] * DP[i+1]
            P[i] -= DP[i]
            if P[i] > PE[i]
                P[i] = (P[i] + DP[i] + PE[i]) / 2.0
            end
        end

        for i in 2:M
            WN[i] = max(WS[i] * (PE[i] / P[i])^B1[i], 1e-7)
            P[i] = PE[i] * (WS[i] / WN[i])^BB[i]
            H[i] = exp(MW * P[i] / (Unitful.R * T[i]))
        end
        H[M+1] = H[M]

        if SE <= im
            break
        end
    end

    SW_out = ((P[2] * K[2] - P[3] * K[3]) / (N1[2] * (Z[3] - Z[2])) + Unitful.gn * K[2] + TR) * dt
    W .= WN
    for i in 2:M+1
        θ_soil[i-1] = WN[i]
    end

    FL_out = EP * (H[2] - rh_loc) / (1.0 - rh_loc) * dt
    humid .= H[2:19]
    potent .= P[2:19]

    for i in 2:M
        PR[i] = -1.0 * (TR * RS[i] - P[i])
    end
    rootpot = PR[2:19]
    leafpot = PL
    trans = TR

    return(
        flux_soil = FL_out,
        flux_leaf = trans,
        θ_soil = θ_soil,
        ψ_soil = potent,
        ψ_root = rootpot,
        ψ_leaf = leafpot,
        rh_soil = humid
    ) 
end