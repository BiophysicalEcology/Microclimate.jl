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
    longwave_out = get_longwave(elev, rh, tair, tsurf, slep, sle, cloud, viewf)
    qrad = u"cal/cm^2/minute"(longwave_out.qrad)

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

function soil_water_balance(
    rh_loc,
    θ_soil,
    ET,
    TEMP,
    depth,
    dt,
    elev,
    rw=2.5E+10u"m^3/kg/s", # resistance per unit length of root, m3 kg-1 s-1
    pc=-1500.0u"J/kg", # critical leaf water potential for stomatal closure, J kg-1
    rl=2000000.0u"m^3/kg/s", # resistance per unit length of leaf, m3 kg-1 s-1
    sp=10.0, # stability parameter, -
    r1=0.001u"m", # root radius, m
    im=1e-6, # maximum overall mass balance error allowed, kg
    maxcount=500
)

    A = zeros(Float64, 19)
    B = zeros(Float64, 19)
    C = zeros(Float64, 19)
    F = zeros(Float64, 19)
    P = zeros(Float64, 19) .* u"J/kg"       # matric potential J/kg
    Z = zeros(Float64, 19)                # depth nodes
    V = zeros(Float64, 19)
    DP = zeros(Float64, 19)
    W = zeros(Float64, 19) .* u"m^3/m^3"    # water content m3/m3
    WN = zeros(Float64, 19) .* u"m^3/m^3"    # water content m3/m3
    K = zeros(Float64, 19) .* u"kg*s/m^3"  # hydraulic conductivity, kg s/m3
    CP = zeros(Float64, 19)
    H = zeros(Float64, 19)
    JV = zeros(Float64, 19)
    DJ = zeros(Float64, 19)
    temp = zeros(Float64, 10)
    θ_soil = zeros(Float64, 18)
    T = zeros(Float64, 19)
    depth = zeros(Float64, 10)
    humid = zeros(Float64, 18)
    potent = zeros(Float64, 18)
    PE = zeros(Float64, 19) .* u"J/kg" # air entry potential J/kg
    KS = zeros(Float64, 19) .* u"kg*s/m^3" # saturated conductivity, kg s/m3
    BB = zeros(Float64, 19) # soil 'b' parameter
    PP = zeros(Float64, 19)
    B1 = zeros(Float64, 19)
    N = zeros(Float64, 19)
    N1 = zeros(Float64, 19)
    WS = zeros(Float64, 19)
    rootpot = zeros(Float64, 18)
    RR = zeros(Float64, 19)
    L = zeros(Float64, 19)
    E = zeros(Float64, 19)
    RS = zeros(Float64, 19)
    PR = zeros(Float64, 19)
    BZ = zeros(Float64, 19)
    BD = zeros(Float64, 19) .* u"Mg/m^3" # soil bulk density, Mg/m3
    DD = zeros(Float64, 19) .* u"Mg/m^3" # soil mineral density, Mg/m3
    M = 18 #number of elements
    P_atmos = get_pressure(elev)

    # Constants
    MW = 0.01801528u"kg/mol" # molar mass of water

    # Convert PE to negative absolute value
    PE .= -abs.(PE)

    # Initialize PP from PE
    PP .= PE
    PP .= -abs.(PP)

    # Saturation water content
    WS .= 1 .- BD ./ DD  # WS = 1 - BD/DD

    # Depth to lower boundary (m)
    Z[M+1] = depth[10] / 100

    # Constants
    WD = 1000.0                  # kg/m³
    DV = 0.000024                # m²/s

    # Soil hydraulic properties
    B1 .= 1.0 ./ BB
    N .= 2.0 .+ 3.0 ./ BB
    N1 .= 1.0 .- N

    # Preparation for wetair call
    WB = 0.0
    DPP = 999.0
    PSTD = 101325.0
    BP = PSTD * (1.0 - (0.0065 * ALTT / 288.0))^(1.0 / 0.190284)

    # Fill Z using provided depth vector
    j = 2
    for I in 3:18
        if isodd(I)
            Z[I] = depth[j] / 100
            j += 1
        else
            Z[I] = Z[I-1] + (depth[j] / 100 - Z[I-1]) / 2
        end
    end

    # Interpolate T from temp
    j = 1
    for I in 1:19
        if isodd(I)
            T[I] = temp[j]
            j += 1
        else
            T[I] = T[I-1] + (temp[j] - T[I-1]) / 2
        end
    end

    # Convert T to Kelvin
    T .+= 273.0

    # Set Z[1] and Z[2] to 0
    Z[1] = 0.0
    Z[2] = 0.0

    # Set initial water content and related variables
    for I in 2:M
        WN[I] = θ_soil[I-1]
        P[I] = PE[I] * (WS[I] / WN[I])^BB[I]
        H[I] = exp(MW * P[I] / (R * T[I-1]))
        K[I] = KS[I] * (PE[I] / P[I])^N[I]
        W[I] = WN[I]
    end

    # Bulk water mass per soil layer
    for I in 2:M
        V[I] = WD * (Z[I+1] - Z[I-1]) / 2
    end
    # Lower boundary condition
    P[M+1] = PE[M] * (WS[M+1] / WS[M+1])^BB[M]
    H[M+1] = 1.0
    W[M+1] = WS[M+1]
    WN[M+1] = WS[M+1]
    Z[1] = -1e10
    Z[M+1] = 1e20
    K[M+1] = KS[M] * (PE[M] / P[M+1])^N[M+1]

    # Initialize root water uptake variables
    RR = zeros(M + 1)
    BZ = zeros(M + 1)
    for I in 2:M
        if L[I] > 0.0
            RR[I] = rw / (L[I] * (Z[I+1] - Z[I-1]) / 2.0)
            BZ[I] = (1.0 - M) * log(π * r1^2 * L[I]) / (4.0 * π * L[I] * (Z[I+1] - Z[I-1]) / 2.0)
        else
            RR[I] = 1e20
            BZ[I] = 0.0
        end
    end

    P[1] = P[2]
    K[1] = 0.0

    # Evapotranspiration
    EP = exp(-0.82 * LAI) * ET
    TP = ET - EP

    # Plant water uptake
    PB = 0.0
    RB = 0.0
    PL = 0.0
    RS = zeros(M + 1)
    for i in 2:M
        RS[i] = BZ[i] / K[i]
        PB += P[i] / (RS[i] + RR[i])
        RB += 1.0 / (RS[i] + RR[i])
    end
    PB /= RB
    RB = 1.0 / RB

    # Newton-Raphson to estimate PL
    count = 0
    while count < maxcount
        if PL > PB
            PL = PB - TP * (RB + rl)
        end
        XP = (PL / pc)^pc
        SL = TP * (RB + rl) * pc * XP / (PL * (1.0 + XP)^2) - 1.0
        FF = PB - PL - TP * (RB + rl) / (1.0 + XP)
        PL -= FF / SL
        count += 1
        if abs(FF) <= 10.0
            break
        end
    end

    TR = TP / (1.0 + XP)
    E = zeros(M + 1)
    for I in 2:M
        E[I] = (P[I] - PL - rl * TR) / (RR[I] + RS[I])
    end

    # Convergence loop
    SE = 0.0
    count = 0
    JV = zeros(M + 1)
    DJ = zeros(M + 1)
    CP = zeros(M + 1)
    A = zeros(M + 1)
    B = zeros(M + 1)
    C = zeros(M + 1)
    F = zeros(M + 1)
    DP = zeros(M + 1)
    while count < maxcount
        SE = 0.0
        count += 1
        for I in 2:M
            K[I] = KS[I] * (PE[I] / P[I])^N[I]
        end

        JV[1] = EP * (H[2] - rh_loc) / (1.0 - rh_loc)
        DJ[1] = EP * MW * H[2] / (Unitful.R * T[1] * (1.0 - rh_loc))

        for I in 2:M
            VP = wet_air(u"K"(T[I]); rh=100.0, P_atmos=P_atmos).ρ_vap
            KV = 0.66 * DV * VP * (WS[I] - (WN[I] + WN[I+1]) / 2.0) / (Z[I+1] - Z[I])
            JV[I] = KV * (H[I+1] - H[I])
            DJ[I] = MW * H[I] * KV / (R * T[I])
            CP[I] = -1.0 * V[I] * WN[I] / (BB[I] * P[I] * dt)
            A[I] = -1.0 * K[I-1] / (Z[I] - Z[I-1]) + Unitful.gn * N[I] * K[I-1] / P[I-1]
            C[I] = -1.0 * K[I+1] / (Z[I+1] - Z[I])
            B[I] = K[I] / (Z[I] - Z[I-1]) + K[I] / (Z[I+1] - Z[I]) + CP[I] - Unitful.gn * N[I] * K[I] / P[I] + DJ[I-1] + DJ[I]
            F[I] = ((P[I] * K[I] - P[I-1] * K[I-1]) / (Z[I] - Z[I-1]) - (P[I+1] * K[I+1] - P[I] * K[I]) / (Z[I+1] - Z[I])) / N1[I] + V[I] * (WN[I] - W[I]) / dt - Unitful.gn * (K[I-1] - K[I]) + JV[I-1] - JV[I] + E[I]
            SE += abs(F[I])
        end

        for I in 2:M-1
            C[I] /= B[I]
            C[I] = C[I] < 1e-8 ? 0.0 : C[I]
            F[I] /= B[I]
            B[I+1] -= A[I+1] * C[I]
            F[I+1] -= A[I+1] * F[I]
        end

        DP[M] = F[M] / B[M]
        P[M] -= DP[M]
        P[M] = min(P[M], PE[M])

        for I in (M-1):-1:2
            DP[I] = F[I] - C[I] * DP[I+1]
            P[I] -= DP[I]
            if P[I] > PE[I]
                P[I] = (P[I] + DP[I] + PE[I]) / 2.0
            end
        end

        for I in 2:M
            WN[I] = max(WS[I] * (PE[I] / P[I])^B1[I], 1e-7)
            P[I] = PE[I] * (WS[I] / WN[I])^BB[I]
            H[I] = exp(MW * P[I] / (Unitful.R * T[I]))
        end
        H[M+1] = H[M]

        if SE <= im
            break
        end
    end

    SW_out[] = ((P[2] * K[2] - P[3] * K[3]) / (N1[2] * (Z[3] - Z[2])) + Unitful.gn * K[2] + TR) * dt
    W .= WN
    for I in 2:M+1
        θ_soil[I-1] = WN[I]
    end

    FL_out[] = EP * (H[2] - rh_loc) / (1.0 - rh_loc) * dt
    humid .= H[2:19]
    potent .= P[2:19]

    for I in 2:M
        rootpot[I-1] = -1.0 * (TR * RS[I] - P[I])
    end
    leafpot[] = PL
    trans[] = TR

    return FL_out, θ_soil, potent, humid, rootpot, leafpot, trans
end