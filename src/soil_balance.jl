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

    # check for unstable conditions of ground surface temperature
    T .= clamp.(T, (-81.0+273.15)u"K", (95.0+273.15)u"K")

    N = length(dep)

    # get soil properties and convert to cal/cm/g/C
    λ_b, cp_b, ρ_b = soil_properties(T, θ_soil, nodes, soilprops, elev)

    # Get environmental data at time t
    tair = f.TAIRt(ustrip(t))
    vel = max(0.1u"m/s", f.VELt(ustrip(t)))
    zenr = min(90u"°", u"°"(round(f.ZENRt(ustrip(t)), digits=3)))
    #solr = u"cal/cm^2/minute"(max(0.0u"W/m^2", f.SOLRt(ustrip(t))))
    solr = max(0.0u"W/m^2", f.SOLRt(ustrip(t)))
    cloud = max(0.0, f.CLDt(ustrip(t)))
    rh = max(0.0, f.RHt(ustrip(t)))
    zslr = min(90u"°", f.ZSLt(ustrip(t)))

    T[N] = tdeep # set boundary condition of deep soil temperature

    depp = fill(0.0u"cm", N + 1)
    depp[1:N] = dep
    # Compute soil layer properties
    wc = fill(1.0u"J/K/m^2", N)
    c = fill(1.0u"W/K/m^2", N)
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
        solr = solr * (0.36 + 0.64 * (1.0 - (cloud / 100.0))) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
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
    qrad = longwave_out.Qrad

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
    P_atmos = get_pressure(elev)
    
    # Evaporation
    wet_air_out = wet_air(u"K"(tair); rh=rh, P_atmos=P_atmos)
    cp_air = wet_air_out.cp
    ρ_air = wet_air_out.ρ_air
    hd = (hc / (cp_air * ρ_air)) * (0.71 / 0.60)^0.666
    qevap, gwsurf = evap(tsurf=u"K"(T[1]), tair=u"K"(tair), rh=rh, rhsurf=100.0, hd=hd, elev=elev, pctwet=pctwet, sat=false)

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
    M = 18,
    PE = fill(1.1, M+1)u"J/kg", # Air entry potential (J/kg) (M+1 values descending through soil for specified soil nodes in parameter DEP and points half way between)
    KS = fill(0.0037, M+1)u"kg*s/m^3", # Saturated conductivity, (kg s/m3) (M+1 values descending through soil for specified soil nodes in parameter DEP and points half way between)
    BB = fill(4.5, M+1), # Campbell's soil 'b' parameter (-) (M+1 values descending through soil for specified soil nodes in parameter DEP and points half way between)
    BD = fill(1.3, M+1)u"Mg/m^3", # Soil bulk density (Mg/m3)  (M+1 values descending through soil for specified soil nodes in parameter DEP and points half way between)
    DD = fill(2.56, M+1)u"Mg/m^3", # Soil density (Mg/m3)  (M+1 values descending through soil for specified soil nodes in parameter DEP and points half way between)
    rh_loc = 0.2,
    θ_soil = fill(0.2, M),
    ET = 1.3e-5u"kg/m^2/s",
    depth = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm",
    T10 = fill(293.15u"K", div(M+2,2)),
    dt = 360u"s",
    elev = 0.0u"m",
    L = [0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0]*1e4u"m/m^3", # root density, m m-3
    rw = 2.5e+10u"m^3/kg/s", # resistance per unit length of root, m3 kg-1 s-1
    pc = -1500.0u"J/kg", # critical leaf water potential for stomatal closure, J kg-1
    rl = 2.0e6u"m^4/kg/s", # leaf resistance, m4 kg-1 s-1
    sp = 10.0, # stability parameter, -
    r1 = 0.001u"m", # root radius, m
    lai = 0.1,
    im = 1e-6u"kg/m^2/s", # maximum overall mass balance error allowed, kg m-2 s-1
    maxcount=500
)

    P = zeros(Float64, M+1)*u"J/kg"       # matric potential J/kg
    Z = zeros(Float64, M+1)u"m"                # depth nodes
    V = zeros(Float64, M+1)u"kg/m^2"
    W = zeros(Float64, M+1)*u"m^3/m^3"    # water content m3/m3
    WN = zeros(Float64, M+1)*u"m^3/m^3"   # water content m3/m3
    K = zeros(Float64, M+1)*u"kg*s/m^3"   # hydraulic conductivity, kg s/m3
    H = zeros(Float64, M+1)
    T = zeros(Float64, M+1)u"K"
    rh_soil = zeros(Float64, M)
    ψ_soil = zeros(Float64, M)u"J/kg"
    ψ_root = zeros(Float64, M)u"J/kg"
    PR = zeros(Float64, M+1)u"J/kg"
    PP = zeros(Float64, M+1)u"J/kg"
    B1 = zeros(Float64, M+1)
    N = zeros(Float64, M+1)
    N1 = zeros(Float64, M+1)
    WS = zeros(Float64, M+1)

    P_atmos = get_pressure(elev)

    # Constants
    MW = 0.01801528u"kg/mol" # molar mass of water, kg/mol
    WD = 1000.0u"kg/m^3"     # density of water, kg/m³
    DV = 2.4e-5u"m^2/s"      # diffusivity of water vapour, m²/s

    # Convert PE to negative absolute value
    PE .= -abs.(PE) # air entry potential J/kg

    # Initialize PP from PE
    PP .= PE
    PP .= -abs.(PP)

    # Saturation water content,  m3/m3
    WS = 1.0 .- BD ./ DD

    # Depth to lower boundary (m)
    Z[M+1] = u"m"(depth[div(M+2,2)])


    # Soil hydraulic properties
    B1 .= 1.0 ./ BB
    N .= 2.0 .+ 3.0 ./ BB
    N1 .= 1.0 .- N

    # Fill Z using provided depth vector
    j = 2
    for i in 3:M
        if isodd(i)
            Z[i] = depth[j]
            j += 1
        else
            Z[i] = Z[i-1] + (depth[j] - Z[i-1]) / 2
        end
    end

    # Interpolate T from temp
    j = 1
    for i in 1:M+1
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
        P[i] = PE[i] * (WS[i] / WN[i])^BB[i] # matric water potential, EQ5.9 (note thetas=W are inverted so not raised to -BB)
        H[i] = exp(MW * P[i] / (Unitful.R * T[i-1])) # fractional humidity, EQ5.14
        K[i] = KS[i] * (PE[i] / P[i])^N[i] # hydraulic conductivity, EQ6.14
        W[i] = θ_soil[i-1] #  water content
    end

    # Bulk water mass per soil layer
    for i in 2:M
        V[i] = WD * (Z[i+1] - Z[i-1]) / 2 # bulk density x volume per unit area, kg/m²
    end
    # Lower boundary condition set to saturated (stays constant)
    P[M+1] = PE[M] * (WS[M+1] / WS[M+1])^BB[M] # water potential
    H[M+1] = 1.0 # fractional humidity
    W[M+1] = WS[M+1] # water content
    WN[M+1] = WS[M+1] # water content
    Z[1] = -1e10u"m" # depth at node 1, m
    Z[M+1] = 1e20u"m" # depth at deepest node, m
    K[M+1] = KS[M] * (PE[M] / P[M+1])^N[M+1] # lower boundary conductivity

    # Initialize root water uptake variables
    RR = zeros(M + 1)u"m^4/kg/s"
    BZ = zeros(M + 1)u"m"
    for i in 2:M
        if L[i] > 0.0u"m/m^3"
            RR[i] = rw / (L[i] * (Z[i+1] - Z[i-1]) / 2.0) # root resistance
            BZ[i] = N1[i] * log(π * r1^2 * L[i]) / (4.0 * π * L[i] * (Z[i+1] - Z[i-1]) / 2.0)
        else
            RR[i] = 1e20u"m^4/kg/s" # root resistance
            BZ[i] = 0.0u"m"
        end
    end

    P[1] = P[2]
    K[1] = 0.0u"kg*s/m^3"

    # Evapotranspiration
    EP = exp(-0.82 * lai) * ET # partition potential evaporation from potential evapotranspiration, EQ12.30
    TP = ET - EP # now get potential transpiration

    # Plant water uptake
    PB1 = 0.0u"J*s/m^4"  # numerator of first term on left of EQ11.18, J * s / m⁴
    RB1 = 0.0u"kg*s/m^4" # weighted mean root-soil resistance, R_bar, m4 /(s kg)
    PL = 0.0u"J/kg"      # leaf water potential, J/kg
    RS = zeros(M + 1)u"m^4/kg/s" # soil resistance, m4 /(s kg)
    for i in 2:M
        RS[i] = BZ[i] / K[i] # soil resistance, simplification of EQ11.14, assuming conductivity constant in the rhizosphere
        PB1 += P[i] / (RS[i] + RR[i]) # summing over layers
        RB1 += 1.0 / (RS[i] + RR[i]) # summing over layers
    end
    PB = PB1 / RB1 # final step in evaluating psi_bar, weighted mean soil water potential, first term on right in EQ11.18
    RB = (1.0 / RB1) # denominator of first and second terms on right in EQ11.18

    # Newton-Raphson to estimate PL
    counter = 0
    while counter < maxcount
        if PL > PB
            PL = PB - TP * (RB + rl) # variation on EQ11.18
        end
        XP = (PL / pc)^sp # part of EQ12.28 determining stomatal closure
        SL = TP * (RB + rl) * sp * XP / (PL * (1.0 + XP)^2) - 1.0 # derivative of stomatal function 
        FF = PB - PL - TP * (RB + rl) / (1.0 + XP) # transpiration mass balance (variation on EQ11.18)
        PL -= FF / SL
        counter += 1
        if abs(FF) <= 10.0u"J/kg"
            break
        end
    end
    PL = u"J/kg"(PL) # keep units in J/kg = m^2/s^2
    XP = (PL / pc)^sp
    TR = TP / (1.0 + XP)
    E = zeros(M + 1)u"kg/m^2/s"
    for i in 2:M
        E[i] = (P[i] - PL - rl * TR) / (RR[i] + RS[i]) # root water uptake, EQ11.15
    end

    # Convergence loop
    JV = zeros(M + 1)u"kg/m^2/s"
    DJ = zeros(M + 1)u"kg*s/m^4"
    CP = zeros(M + 1)u"kg*s/m^4"
    A = zeros(M + 1)u"kg*s/m^4"
    B = zeros(M + 1)u"kg*s/m^4"
    C = zeros(M + 1)u"kg*s/m^4"
    C2 = zeros(M + 1) # adding this to deal with ratio check
    F = zeros(M + 1)u"kg/m^2/s"
    F2 = zeros(M + 1)u"J/kg"
    DP = zeros(M + 1)u"J/kg"
    counter = 0
    while counter < maxcount
        SE = 0.0u"kg/m^2/s"
        counter += 1
        for i in 2:M
            K[i] = KS[i] * (PE[i] / P[i])^N[i]
        end

        JV[1] = EP * (H[2] - rh_loc) / (1.0 - rh_loc) # vapour flux at soil surface, EQ9.14
        DJ[1] = EP * MW * H[2] / (Unitful.R * T[1] * (1.0 - rh_loc)) # derivative of vapour flux at soil surface, combination of EQ9.14 and EQ5.14

        for i in 2:M
            VP = wet_air(u"K"(T[i]); rh=100.0, P_atmos=P_atmos).ρ_vap # VP is vapour density = c'_v in EQ9.7
            KV = 0.66 * DV * VP * (WS[i] - (WN[i] + WN[i+1]) / 2.0) / (Z[i+1] - Z[i]) # vapour conductivity, EQ9.7, assuming epsilon(psi_g) = b*psi_g^m (eq. 3.10) where b = 0.66 and m = 1 (p.99)
            JV[i] = KV * (H[i+1] - H[i]) # fluxes of vapour within soil, EQ9.14
            DJ[i] = MW * H[i] * KV / (Unitful.R * T[i-1]) # derivatives of vapour fluxes within soil, combination of EQ9.14 and EQ5.14
            CP[i] = -1.0 * V[i] * WN[i] / (BB[i] * P[i] * dt) # hydraulic capacity = capacitance, d_theta/d_psi
            # Jacobian components
            A[i] = -1.0 * K[i-1] / (Z[i] - Z[i-1]) + Unitful.gn * N[i] * K[i-1] / P[i-1] # sub-diagonal element in tridagonal matrix
            C[i] = -1.0 * K[i+1] / (Z[i+1] - Z[i]) # super-diagonal element in tridagonal matrix
            B[i] = K[i] / (Z[i] - Z[i-1]) + K[i] / (Z[i+1] - Z[i]) + CP[i] - Unitful.gn * N[i] * K[i] / P[i] + DJ[i-1] + DJ[i] # diagonal element in tridagonal matrix
            # mass balance including vapour fluxes and root water uptake
            # version of equation 8.28 that additionally conatins vapour fluxes and root water uptake
            F[i] = ((P[i] * K[i] - P[i-1] * K[i-1]) / (Z[i] - Z[i-1]) - (P[i+1] * K[i+1] - P[i] * K[i]) / (Z[i+1] - Z[i])) / N1[i] + V[i] * (WN[i] - W[i]) / dt - Unitful.gn * (K[i-1] - K[i]) + JV[i-1] - JV[i] + E[i]
            SE += abs(F[i]) # total mass balance error
        end
        
        # Thomas algorithm (Gauss elimination)
        for i in 2:M-1
            C2[i] = C[i] / B[i]
            #C[i] = C2[i] < 1e-8 ? 1e-8u"kg*s/m^4" : C[i]
            #C2[i] = C2[i] < 1e-8 ? 1e-8 : C2[i]
            F2[i] = F[i] / B[i]
            B[i+1] -= A[i+1] * C2[i]
            F[i+1] -= A[i+1] * F2[i]
        end

        DP[M] = F[M] / B[M]
        P[M] -= DP[M]
        P[M] = min(P[M], PE[M])

        for i in (M-1):-1:2
            DP[i] = F2[i] - C2[i] * DP[i+1] # change in matric potential in an interation step, J/kg
            P[i] -= DP[i]                   # matric potential, J/kg
            if P[i] > PE[i]
                P[i] = (P[i] + DP[i] + PE[i]) / 2.0
            end
        end

        for i in 2:M
            WN[i] = max(WS[i] * (PE[i] / P[i])^B1[i], 1e-7)
            P[i] = PE[i] * (WS[i] / WN[i])^BB[i]
            H[i] = exp(MW * P[i] / (Unitful.R * T[i-1]))
        end
        H[M+1] = H[M]
        if SE <= im
            break
        end
    end

    SW = ((P[2] * K[2] - P[3] * K[3]) / (N1[2] * (Z[3] - Z[2])) + Unitful.gn * K[2] + TR) * dt
    W .= WN
    for i in 2:M+1
        θ_soil[i-1] = WN[i]
    end

    for i in 2:M
        PR[i] = -1.0 * (TR * RS[i] - P[i])
    end
    
    return(
        evap = EP * (H[2] - rh_loc) / (1.0 - rh_loc) * dt,
        trans = TR,
        θ_soil = θ_soil,
        ψ_soil = P[2:(M+1)],
        ψ_root = PR[2:(M+1)],
        ψ_leaf = PL,
        rh_soil = H[2:(M+1)],
        Δ_H2O = SW,
        drain = Unitful.gn * K[M]
    ) 
end