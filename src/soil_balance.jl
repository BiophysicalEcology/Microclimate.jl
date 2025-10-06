# This is a 3-parameters OrdinaryDiffEq function
function soil_energy_balance(
    T::U,                # state
    i::SoilEnergyInputs, # "tarameters"
    t::Quantity,          # timestep
) where U <: SVector{N} where N
    #t_min = t / 60 * u"minute"  # convert Float64 time back to unitful
    #T_K = T .* u"K"  # convert Float64 time back to unitful
    #dT_K = dT .* 60 .* u"K/minute"  # convert Float64 time back to unitful
    # extract prameters
    (; forcings, buffers, heights, depths, nodes, soilprops, environment_instant, terrain, runmoist) = params
    (; depp, wc, c) = buffers.soil_energy
    # Get environmental data at time t
    (; tair, vel, zenr, solr, cloud, rh, zslr) = interpolate_forcings(forcings, t)

    reference_height = last(heights)
    sabnew = 1.0 - environment_instant.albedo

    # check for unstable conditions of ground surface temperature
    T1 = map(t -> clamp(t, (-81.0+273.15)u"K", (85.0+273.15)u"K"), T)::U

    # get soil properties
    (; λ_b, cp_b, ρ_b) = soil_props_vector!(buffers.soil_properties; T_soil=T1, θ_soil=θ_soil, soilprops, terrain)
    # TODO Why do we reset the last value
    T1m = MVector(T1)
    T1m[N] = tdeep # boundary condition
    T2 = SVector(T1m)

    # set boundary condition of deep soil temperature
    depp[1:N] = depths
    # Compute soil layer properties
    @inbounds for i in 1:N
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
    Q_solar = sabnew * solr * ((100.0 - shade) / 100.0)
    if slope > 0 && zenr < 90u"°"
        cz = cosd(zenr)
        czsl = cosd(zslr)
        Q_solar = (Q_solar / cz) * czsl
    end

    # Longwave radiation
    longwave_out = get_longwave(; terrain, tsurf=T2[1], environment_instant)
    Q_infrared = longwave_out.Qrad

    # Conduction
    Q_conduction = c[1] * (T2[2] - T2[1])

    # Convection
    log_z_ratio = log(reference_height / roughness_height + 1)
    T_ref_height = tair
    T_surface = T2[1]
    ΔT = T_ref_height - T_surface
    T_mean = (T_surface + T_ref_height) / 2
    # TODO call calc_ρ_cp method specific to elevation and RH in final version but do it this way for NicheMapR comparison
    ρ_cp = calc_ρ_cp(T_mean)#, elevation, relative_humidity)
    if T_ref_height ≥ T_surface || zenr ≥ 90°
        u_star = calc_u_star(; reference_wind_speed=vel, log_z_ratio, κ)
        Q_convection = calc_convection(; u_star, log_z_ratio, ΔT, ρ_cp, z0=roughness_height)
    else
        # compute ρcpTκg (was a constant in original Fortran version)
        # dry_air_out = dry_air_properties(u"K"(reference_temperature), elevation=elevation)
        # wet_air_out = wet_air_properties(u"K"(reference_temperature), rh = relative_humidity)
        #ρ = dry_air_out.ρ_air
        #c_p = wet_air_out.c_p
        # TODO make this work with SI units
        #ρcpTκg = u"cal*minute^2/cm^4"(ρ * c_p * T_ref_height / (κ * g_n))
        ρcpTκg = 6.003e-8u"cal*minute^2/cm^4"
        L_Obukhov = -30.0u"cm" # initialise Obukhov length
        Obukhov_out = calc_Obukhov_length(T_ref_height, T_surface, vel, roughness_height, reference_height, ρcpTκg, κ, log_z_ratio, ΔT, ρ_cp)
        L_Obukhov = Obukhov_out.L_Obukhov
        Q_convection = Obukhov_out.Q_convection
    end
    hc = max(abs(Q_convection / (T2[1] - tair)), 0.5u"W/m^2/K")

    # Evaporation
    wet_air_out = wet_air_properties(u"K"(tair); rh, terrain.P_atmos)
    c_p_air = wet_air_out.c_p
    ρ_air = wet_air_out.ρ_air
    hd = (hc / (c_p_air * ρ_air)) * (0.71 / 0.60)^0.666
    Q_evaporation, gwsurf = evaporation(; tsurf=u"K"(T[1]), tair=u"K"(tair), rh, rhsurf=100.0, hd, terrain, pctwet, sat=false)

    # Construct static vector of change in soil temperature, to return
    # Energy balance at surface
    surface = u"K/minute"((Q_solar + Q_infrared + Q_conduction + Q_convection - Q_evaporation) / wc[1])
    # Soil conduction for internal nodes
    internal = ntuple(Val{N-2}()) do i
        u"K/minute".((c[i] * (T2[i] - T2[i+1]) + c[i+1] * (T2[i+2] - T2[i+1])) / wc[i+1])
    end
    # Lower boundary condition
    lower_boundary = 0.0u"K/minute"  # or set T[N] = T_surface from data
    dT = SVector{N}((surface, internal..., lower_boundary))

    return dT
end

function interpolate_forcings(f, t)
    t_m = ustrip(u"minute", t)
    return (; 
        tair = f.TAIRt(t_m)
        vel = max(0.1u"m/s", f.VELt(t_m))
        zenr = min(90.0u"°", u"°"(round(f.ZENRt(t_m), digits=3)))
        solr = max(0.0u"W/m^2", f.SOLRt(t_m))
        cloud = clamp(f.CLDt(t_m), 0.0, 100.0)
        rh = clamp(f.RHt(t_m), 0.0, 100.0)
        zslr = min(90.0u"°", f.ZSLt(t_m))
    )
end

function evaporation(; terain, tsurf, tair, rh, rhsurf, hd, pctwet, sat)
    (; elevation, P_atmos) = terrain
    # Assumes all units are SI (Kelvin, Pascal, meters, seconds, kg, etc.)

    # Ground-level variables, shared via global or passed in as needed
    #global shayd, altt, maxshd, sabnew, ptwet, rainfall

    # Output variables
    Q_evaporation = 0.0
    gwsurf = 0.0

    tsurf = tsurf < u"K"(-81.0u"°C") ? u"K"(-81.0u"°C") : tsurf

    # surface and air vapor densities
    ρ_vap_surf = wet_air_properties(u"K"(tsurf); rh=rhsurf, P_atmos).ρ_vap
    ρ_vap_air = wet_air_properties(u"K"(tair); rh, P_atmos).ρ_vap

    # Effective wet surface fraction
    effsur = sat ? 1.0 : pctwet / 100.0

    # Water evaporated from surface (kg/s/m^2)
    water = effsur * hd * (ρ_vap_surf - ρ_vap_air)

    # Latent heat of vaporization (J/kg)
    λ_evap = enthalpy_of_vaporisation(tsurf) 

    # Energy flux due to evaporation (W/m² or converted)
    Q_evaporation = water * λ_evap  # SI units for water budget calcs

    # Mass flux (g/s)
    gwsurf = u"g/s/m^2"(water)

    # No water loss if TSURF ≤ 0 (e.g., melting snow only)
    if tsurf <= u"K"(0.0u"°C")
        gwsurf = 0.0u"g/s/m^2"
    end

    return Q_evaporation, gwsurf
end

function allocate_soil_water_balance(M)
    (;
        P = zeros(M+1) .* u"J/kg",
        Z = zeros(M+1) .* u"m",
        V = zeros(M+1) .* u"kg/m^2",
        W = zeros(M+1) .* u"m^3/m^3",
        WN = zeros(M+1) .* u"m^3/m^3",
        K = zeros(M+1) .* u"kg*s/m^3",
        H = zeros(M+1),
        T = zeros(M+1) .* u"K",
        rh_soil = zeros(M),
        ψ_soil = zeros(M) .* u"J/kg",
        ψ_root = zeros(M) .* u"J/kg",
        PR = zeros(M+1) .* u"J/kg",
        PP = zeros(M+1) .* u"J/kg",
        B1 = zeros(M+1),
        N = zeros(M+1),
        N1 = zeros(M+1),
        WS = zeros(M+1),
        RR = zeros(M+1) .* u"m^4/kg/s",
        BZ = zeros(M+1) .* u"m",
        JV = zeros(M+1) .* u"kg/m^2/s",
        DJ = zeros(M+1) .* u"kg*s/m^4",
        CP = zeros(M+1) .* u"kg*s/m^4",
        A = zeros(M+1) .* u"kg*s/m^4",
        B = zeros(M+1) .* u"kg*s/m^4",
        C = zeros(M+1) .* u"kg*s/m^4",
        C2 = zeros(M+1),
        F = zeros(M+1) .* u"kg/m^2/s",
        F2 = zeros(M+1) .* u"J/kg",
        DP = zeros(M+1) .* u"J/kg",
        RS = zeros(typeof(0.0u"m^4/kg/s"), M + 1), # soil resistance, m4 /(s kg)
        E = zeros(typeof(0.0u"kg/m^2/s"), M + 1),
        # The output is resized for these
        P_out = Vector{typeof(1.0u"J/kg")}(undef, M),
        H_out = Vector{Float64}(undef, M),
        PR_out = Vector{typeof(1.0u"J/kg")}(undef, M),
    )  
end

function soil_water_balance!(buffers, smm::SoilMoistureModel;
    depth,
    terrain,
    soil_moisture,
    rh_loc,
    ET,
    T10
    lai,
)
    dt = smm.moist_step
    PE = smm.air_entry_water_potential
    KS = smm.saturated_hydraulic_conductivity
    BB = smm.Campbells_b_parameter
    BD = smm.soil_bulk_density2
    DD = smm.soil_mineral_density2
    L = smm.root_density
    rw = smm.root_resistance
    pc = smm.stomatal_closure_potential
    rl = smm.leaf_resistance
    sp = smm.stomatal_stability_parameter
    r1 = smm.root_radius
    im = smm.moist_error
    moist_count = smm.moist_count

    # TODO: some of these are actually buffers, and some user data??
    (; P, Z, V, W, WN, K, H, T, rh_soil, ψ_soil, ψ_root, PR, PP, B1, N, N1, WS,
       RR, BZ, JV, DJ, CP, A, B, C, C2, F, F2, DP, RS, E) = buffers
    M = length(P)

    # Constants
    MW = 0.01801528u"kg/mol" # molar mass of water, kg/mol
    WD = 1000.0u"kg/m^3"     # density of water, kg/m³
    DV = 2.4e-5u"m^2/s"      # diffusivity of water vapour, m²/s

    # Convert PE to negative absolute value
    map!(x -> -abs(x), PE, PE) # air entry potential J/kg

    # Initialize PP from PE
    PP .= PE

    # Saturation water content,  m3/m3
    WS .= 1.0 .- BD ./ DD

    # Depth to lower boundary (m)
    Z[M+1] = depth[div(M+2,2)]


    # Soil hydraulic properties
    B1 .= 1.0 ./ BB
    N .= 2.0 .+ 3.0 ./ BB
    N1 .= 1.0 .- N

    # Fill Z using provided depth vector
    j = 2
    @inbounds for i in 3:M
        if isodd(i)
            Z[i] = depth[j]
            j += 1
        else
            Z[i] = Z[i-1] + (depth[j] - Z[i-1]) / 2
        end
    end

    # Interpolate T from temp
    j = 1
    @inbounds for i in 1:M+1
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
    @inbounds for i in 2:M
        WN[i] = θ_soil[i-1]
        P[i] = PE[i] * (WS[i] / WN[i])^BB[i] # matric water potential, EQ5.9 (note thetas=W are inverted so not raised to -BB)
        H[i] = exp(MW * P[i] / (R * T[i-1])) # fractional humidity, EQ5.14
        K[i] = KS[i] * (PE[i] / P[i])^N[i] # hydraulic conductivity, EQ6.14
        W[i] = θ_soil[i-1] #  water content
    end

    # Bulk water mass per soil layer
    @inbounds for i in 2:M
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
    @inbounds for i in 2:M
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
    EP = exp(-0.82 * ustrip(lai)) * ET # partition potential evaporation from potential evapotranspiration, EQ12.30
    TP = ET - EP # now get potential transpiration

    # Plant water uptake
    PB1 = 0.0u"J*s/m^4"  # numerator of first term on left of EQ11.18, J * s / m⁴
    RB1 = 0.0u"kg*s/m^4" # weighted mean root-soil resistance, R_bar, m4 /(s kg)
    PL = 0.0u"J/kg"      # leaf water potential, J/kg
    @inbounds for i in 2:M
        RS[i] = BZ[i] / K[i] # soil resistance, simplification of EQ11.14, assuming conductivity constant in the rhizosphere
        PB1 += P[i] / (RS[i] + RR[i]) # summing over layers
        RB1 += 1.0 / (RS[i] + RR[i]) # summing over layers
    end
    PB = PB1 / RB1 # final step in evaluating psi_bar, weighted mean soil water potential, first term on right in EQ11.18
    RB = (1.0 / RB1) # denominator of first and second terms on right in EQ11.18

    # Newton-Raphson to estimate PL
    counter = 0
    while counter < moist_count
        if PL > PB
            # Seems we need to force the units here or PL is type unstable in the loop
            PL = uconvert(u"J/kg", PB - TP * (RB + rl)) # variation on EQ11.18
        end
        XP = (PL / pc)^sp # part of EQ12.28 determining stomatal closure
        SL = TP * (RB + rl) * sp * XP / (PL * (1.0 + XP)^2) - 1.0 # derivative of stomatal function 
        FF = PB - PL - TP * (RB + rl) / (1.0 + XP) # transpiration mass balance (variation on EQ11.18)
        PL = uconvert(u"J/kg", PL - (FF / SL))
        counter += 1
        if abs(FF) <= 10.0u"J/kg"
            break
        end
    end
    XP = (PL / pc)^sp
    TR = TP / (1.0 + XP)
    @inbounds for i in 2:M
        E[i] = (P[i] - PL - rl * TR) / (RR[i] + RS[i]) # root water uptake, EQ11.15
    end

    # Convergence loop
    counter = 0
    while counter < moist_count
        SE = 0.0u"kg/m^2/s"
        counter += 1
        @inbounds @simd for i in 2:M
            K[i] = KS[i] * (PE[i] / P[i])^N[i]
        end

        JV[1] = EP * (H[2] - rh_loc) / (1.0 - rh_loc) # vapour flux at soil surface, EQ9.14
        DJ[1] = EP * MW * H[2] / (R * T[1] * (1.0 - rh_loc)) # derivative of vapour flux at soil surface, combination of EQ9.14 and EQ5.14

        @inbounds @simd for i in 2:M
            VP = wet_air_properties(u"K"(T[i]); rh=100.0, P_atmos).ρ_vap # VP is vapour density = c'_v in EQ9.7
            KV = 0.66 * DV * VP * (WS[i] - (WN[i] + WN[i+1]) / 2.0) / (Z[i+1] - Z[i]) # vapour conductivity, EQ9.7, assuming epsilon(psi_g) = b*psi_g^m (eq. 3.10) where b = 0.66 and m = 1 (p.99)
            JV[i] = KV * (H[i+1] - H[i]) # fluxes of vapour within soil, EQ9.14
            DJ[i] = MW * H[i] * KV / (R * T[i-1]) # derivatives of vapour fluxes within soil, combination of EQ9.14 and EQ5.14
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
        @inbounds @simd for i in 2:M-1
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

        @inbounds for i in (M-1):-1:2
            DP[i] = F2[i] - C2[i] * DP[i+1] # change in matric potential in an interation step, J/kg
            P[i] -= DP[i]                   # matric potential, J/kg
            if P[i] > PE[i]
                P[i] = (P[i] + DP[i] + PE[i]) / 2.0
            end
        end

        @inbounds @simd for i in 2:M
            WN[i] = max(WS[i] * (PE[i] / P[i])^B1[i], 1e-7)
            P[i] = PE[i] * (WS[i] / WN[i])^BB[i]
            H[i] = exp(MW * P[i] / (R * T[i-1]))
        end
        H[M+1] = H[M]
        if SE <= im
            break
        end
    end

    SW = ((P[2] * K[2] - P[3] * K[3]) / (N1[2] * (Z[3] - Z[2])) + Unitful.gn * K[2] + TR) * dt
    W .= WN
    @inbounds for i in 2:M+1
        θ_soil[i-1] = WN[i]
    end

    @inbounds for i in 2:M
        PR[i] = -1.0 * (TR * RS[i] - P[i])
    end
    
    evap = EP * (H[2] - rh_loc) / (1.0 - rh_loc) * dt
    return (;
        evap,
        trans = TR,
        θ_soil,
        ψ_leaf = PL,
        # These need the first value removed. Why?
        ψ_soil = (ml.P_out .= view(P, 2:(M+1))),
        ψ_root = (ml.PR_out .= view(PR, 2:(M+1))),
        rh_soil = (ml.H_out .= view(H, 2:(M+1))),
        Δ_H2O = SW,
        drain = Unitful.gn * K[M]
    ) 
end


get_soil_water_balance(; M=18, kw...) = 
    get_soil_water_balance!(allocate_soil_water_balance(M); kw...)

function get_soil_water_balance!(buffers, soil_moisture_model::SoilMoistureModel;
    heights,
    depths,
    terrain,
    environment_instant,
    T0,
    pool,
    niter_moist,
    pctwet,
)
    P_atmos = terrain.P_atmos
    ei = environment_instant
    tair = ei.air_temperature
    tsurf = ei.surface_temerature
    rh = ei.relative_humidity

    enviroment_instant = (;
        reference_temperature = TAIRs[step],
        reference_wind_speed = VELs[step],
        relative_humidity = RHs[step],
        surface_temperature = u"°C"(T0[1]),  # top layer temp
        zenith_angle = ZENRs[step],
    )
    # compute scalar profiles
    profile_out = get_profile!(buffers.profile;
        terrain,
        environment_instant,
    )

    # convection
    Q_convection = profile_out.Q_convection

    # evaporation
    rh_loc = min(0.99, profile_out.humidities[2] / 100)
    hc = max(abs(Q_convection / (tsurf - tair), 0.5u"W/m^2/K")
    wet_air_out = wet_air_properties(tair; rh, P_atmos)
    c_p_air = wet_air_out.c_p
    ρ_air = wet_air_out.ρ_air
    hd = (hc / (c_p_air * ρ_air)) * (0.71 / 0.60)^0.666
    Q_evaporation, gwsurf = evaporation(; tsurf, tair, rh, rhsurf=100.0, hd, terrain, pctwet, sat=true)
    λ_evap = enthalpy_of_vaporisation(tsurf)
    EP = max(1e-7u"kg/m^2/s", Q_evaporation / λ_evap) # evaporation potential, mm/s (kg/m2/s)

    if pool > 0.0u"kg/m^2" # surface is wet - saturate it for infiltration
        θ_soil0_b[1] = 1 - BD[1] / DD[1]
    end
    # run infiltration algorithm
    infil_out = soil_water_balance!(buffers.soil_water_balance;
        terrain,
        rh_loc, 
        θ_soil=θ_soil0_b,
        ET=EP,
        T10=T0,
        depth=depths,
        dt=moist_step,
    )
    θ_soil0_b = infil_out.θ_soil
    surf_evap = max(0.0u"kg/m^2", infil_out.evap)
    Δ_H2O = max(0.0u"kg/m^2", infil_out.Δ_H2O)
    pool = clamp(pool - Δ_H2O - surf_evap, 0.0u"kg/m^2", maxpool) # pooling surface water
    if pool > 0.0u"kg/m^2" # surface is wet - saturate it for infiltration
        θ_soil0_b[1] = 1 - BD[1] / DD[1]
    end
    for _ in 1:(niter_moist-1)
        infil_out = soil_water_balance!(buffers.soil_water_balance;
            terrain, 
            rh_loc,
            θ_soil=θ_soil0_b,
            ET=EP,
            T10=T0,
            depth,
            dt=moist_step,
        )
        θ_soil0_b = infil_out.θ_soil
        surf_evap = max(0.0u"kg/m^2", infil_out.evap)
        Δ_H2O = max(0.0u"kg/m^2", infil_out.Δ_H2O)
        pool = clamp(pool - Δ_H2O - surf_evap, 0.0u"kg/m^2", maxpool)
        if pool > 0.0u"kg/m^2"
            θ_soil0_b[1] = 1 - BD[1] / DD[1]
        end
    end
    pctwet = clamp(abs(surf_evap / (EP * moist_step) * 100), 0, 100)

    return (; infil_out, pctwet, pool, θ_soil0_b)
end

function allocate_phase_transition(nodes)
    layermass = zeros(Float64, nodes)u"kg"
    qphase = zeros(Float64, nodes)u"J"
    return (; layermass, qphase)
end

phase_transition(; depths, kw...) = 
    phase_transition!(allocate_phase_transition(length(depths)); depths, kw...)

function phase_transition!(buffers::NamedTuple; 
    Ts::AbstractVector,       # current temps at nodes
    T_past::AbstractVector,  # temps at previous step
    ∑phase::AbstractVector,  # accumulated latent heat
    θ::AbstractVector,       # soil moisture by layer
    depths::AbstractVector   # soil depth boundaries (cm)
)
    (; layermass, qphase) = buffers
    HTOFN = 333500.0u"J/kg" # latent heat of fusion of waterper unit mass
    c_p = 4186.0u"J/kg/K" # specific heat of water
    nodes = length(depths)
    meanT = similar(Ts)
    meanTpast = similar(Ts)
    T = MVector(Ts)

    for j in 1:nodes
        if θ[j] > 0.0
            if j < nodes
                meanT[j] = (T[j] + T[j+1]) / 2.0
                meanTpast[j] = (T_past[j] + T_past[j+1]) / 2.0
            else
                meanT[j] = T[j]
                meanTpast[j] = T_past[j]
            end

            if meanTpast[j] > 273.15u"K" && meanT[j] <= 273.15u"K"
                if j < nodes
                    layermass[j] = u"m"(depths[j+1] - depths[j]) * 1000.0u"kg/m" * θ[j]
                else
                    layermass[j] = u"m"(depths[j] + 100.0u"cm" - depths[j]) * 1000.0u"kg/m" * θ[j]
                end

                qphase[j] = (meanTpast[j] - meanT[j]) * layermass[j] * c_p
                ∑phase[j] += qphase[j]

                if ∑phase[j] > HTOFN * layermass[j]
                    # Fully frozen
                    T[j] = 273.14u"K"
                    if j < nodes
                        T[j+1] = 273.14u"K"
                    end
                    ∑phase[j] = 0.0u"J"
                else
                    # In the process of freezing
                    T[j] = 273.16u"K"
                    if j < nodes
                        T[j+1] = 273.16u"K"
                    end
                end
            end
        end
    end

    return (; ∑phase, qphase, T=SVector(T))
end
