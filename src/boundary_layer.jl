function get_profile(;
    refhyt=1.2u"m",
    ruf=0.004u"m",
    zh=0.0u"m",
    d0=0.0u"m",
    κ=0.4, # Kármán constant
    TAREF=27.77818u"°C",
    VREF=2.749575u"m/s",
    rh=49.0415,
    D0cm=48.58942u"°C",
    maxsurf=95.0u"°C",
    ZEN=21.50564u"°",
    a=0.15,
    heights=[0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0] * u"m",
    elev=0.0u"m",
    warn=false)

    if minimum(heights) < ruf
        error("ERROR: the minimum height is not greater than the roughness height (ruf).")
    end

    addheight = false
    heights_extra = nothing
    if minimum(heights) > refhyt
        addheight = true
        heights = vcat([0.01], heights)
    end

    if maximum(heights) >= refhyt && warn
        println("Warning: some heights are ≥ reference height. Using constant air temperature and adjusting wind speed.")
    end

    heights_orig = copy(heights)
    heights = filter(h -> h < refhyt, heights_orig)

    function RHOCP(TAVE)
        return u"(cal*g)/(g*cm^3*K)"*(0.08472 / ustrip(TAVE))
    end
    function RHOCP(TAVE, elev, rh)
        dry_air_out = dry_air(u"K"(TAVE), elev=elev)
        wet_air_out = wet_air(u"K"(TAVE), rh=rh)
        ρ = dry_air_out.ρ_air
        cp = wet_air_out.cp
        return u"(cal*g)/(g*cm^3*K)"(ρ * cp)
    end
    function PHI(z, GAM, AMOL)
        return (1 - min(1, GAM * ustrip(z / AMOL)))^0.25
    end

    function PSI1(X)
        return 2 * log((1 + X) / 2) + log((1 + X^2) / 2) - 2 * atan(X) + π / 2
    end

    function PSI2(X)
        return 2 * log((1 + X^2) / 2)
    end

    function get_Obukhov(TA, TS, V, z, z0, rcptkg, κ)

        AMOL = -30.0u"cm" # initial Monin-Obukhov length cm
        GAM = 16.0 # -
        #RCPTKG = 6.003e-8u"cal/minute/cm/K" #CAL-MIN-CM-C
        z = u"cm"(z)
        z0 = u"cm"(z0)
        ZRATIO = z / z0 + 1
        DUM = log(ZRATIO)
        TA = u"K"(TA)
        TS = u"K"(TS)
        DIFFT = TA - TS
        TAVE = (TA + TS) / 2.0
        RCP = RHOCP(TAVE)
        DEL = 1.0
        count = 0
        USTAR = 0.0u"cm/minute"
        QC = 0.0u"cal/minute/cm^2"
        STO = 0.0
        STB = 0.0
        STS = 0.0

        while DEL > 1e-2 && count < 500
            count += 1
            X = PHI(z, GAM, AMOL)
            Y = PSI1(X)
            USTAR = κ * V / (log(z / z0) - Y)

            if AMOL > 0.0u"cm"
                STS = 0.62 / (ustrip(z0) * ustrip(USTAR) / 12)^0.45
                STB = 0.64 / DUM
                QC = RCP * DIFFT * USTAR * STB / (1.0 + STB / STS)
            else
                STS = 0.62 / (ustrip(z0) * ustrip(USTAR) / 12)^0.45
                STB = (0.64 / DUM) * (1 - 0.1 * z / AMOL)
                STO = STB / (1 + STB / STS)
                QC = RCP * DIFFT * USTAR * STO
            end

            AMOLN = rcptkg * USTAR^3 / QC
            DEL = abs((AMOLN - AMOL) / AMOL)
            AMOL = AMOLN
        end

        return (; AMOL=u"m"(AMOL), STS, STO, STB, USTAR, QC)
    end

    T1 = u"K"(TAREF)
    T3 = u"K"(D0cm)

    # Units: m to cm
    z = u"cm"(refhyt)
    z0 = u"cm"(ruf)
    zh_cm = u"cm"(zh)
    d0_cm = u"cm"(d0)
    V = u"cm/minute"(VREF)
    # define air heights
    AIRDP = vcat(z, reverse(u"cm".(heights)))
    ZZ = AIRDP
    NAIR = length(AIRDP)
    VV = (zeros(Float64, NAIR)) .* 1u"cm/minute" # output wind speeds
    T = Vector{typeof(0.0u"K")}(undef, NAIR) # output temperatures, need to do this otherwise get InexactError
    RHs = zeros(Float64, NAIR) # output relative humidities
    VV[1] = V
    T[1] = T1

    # compute rcptkg (was a constant in original Fortran version)
    dry_air_out = dry_air(u"K"(TAREF), elev=elev)
    wet_air_out = wet_air(u"K"(TAREF), rh=rh)
    ρ = dry_air_out.ρ_air
    cp = wet_air_out.cp
    g = 9.80665u"m/s^2"
    TREF = u"K"(TAREF)
    rcptkg = u"cal*minute^2/cm^4"(ρ * cp * TREF / (κ * g))
    #rcptkg = 6.003e-8u"cal*minute^2/cm^4"
    GAM = 16
    ZRATIO = z / z0 + 1.0
    DUM = log(ZRATIO)
    USTAR = κ * V / DUM
    DIFFT = T1 - T3
    TAVE = (T3 + T1) / 2
    RCP = RHOCP(TAVE, elev, rh)
    AMOL = -30.0u"cm"
    if zh > 0.0u"m"
        STS = 0.62 / (ustrip(z0) * ustrip(USTAR) / 12)^0.45
        STB = 0.64 / DUM
        QC = RCP * DIFFT * USTAR * STB / (1.0 + STB / STS)

        for i in 2:NAIR
            if T1 ≥ T3 || T3 ≤ u"K"(maxsurf) || ZEN ≥ 90°
                VV[i] = (USTAR / κ) * log(ZZ[i] / z0 + 1)
            else
                X1 = PHI(ZZ[i], GAM, AMOL)
                Y1 = PSI1(X1)
                ADUM = ZZ[i] / z0 - Y1
                VV[i] = (USTAR / κ) * log(ADUM)
            end

            A = (T1 - T3) / (1 - log((z - d0_cm) / zh_cm))
            T0 = T1 + A * log((z - d0_cm) / zh_cm)
            T[i] = T0 - A * log((ZZ[i] - d0_cm) / zh_cm)
        end
    else
        if T1 ≥ T3 || T3 ≤ u"K"(maxsurf) || ZEN ≥ 90°
            STS = 0.62 / (ustrip(z0) * ustrip(USTAR) / 12.)^0.45
            STB = 0.64 / DUM
            QC = RCP * DIFFT * USTAR * STB / (1.0 + STB / STS)

            for i in 2:NAIR
                VV[i] = (USTAR / κ) * log(ZZ[i] / z0 + 1.0)
                TZO = (T1 * STB + T3 * STS) / (STB + STS)
                T[i] = TZO + (T1 - TZO) * log(ZZ[i] / z0 + 1.0) / DUM
            end
        else
            for i in 2:NAIR
                X1 = PHI(ZZ[i], GAM, AMOL)
                Y1 = PSI1(X1)
                YY2 = PSI2(X1)
                X = PHI(z, GAM, AMOL)
                #Y = PSI1(X)
                YY = PSI2(X)
                ADUM = ZZ[i] / z0 - Y1
                VV[i] = (USTAR / κ) * log(ADUM)

                Obukhov_out = get_Obukhov(T1, T3, V, ZZ[i], z0, rcptkg, κ)
                TZO = (T1 * Obukhov_out.STB + T3 * Obukhov_out.STS) / (Obukhov_out.STB + Obukhov_out.STS)
                T[i] = TZO + (T1 - TZO) * log(ZZ[i] / z0 - YY2) / log(z / z0 - YY)
            end
        end
    end

    heights = [0.0u"cm"; reverse(ZZ); u"cm"(refhyt)]
    VV = [0.0u"cm/minute"; reverse(VV)]
    T = [T3; reverse(T)]

    e = wet_air(T1, rh=rh).P_vap
    wet_air_out = wet_air.(T; rh=rh)
    es = getproperty.(wet_air_out, :P_vap_sat)
    RHs = clamp.(e ./ es .* 100, 0, 100)

    if heights_extra !== nothing
        VV_extra = V .* (heights_extra ./ refhyt) .^ a
        T_extra = fill(TAREF, length(heights_extra))
        RH_extra = fill(rh, length(heights_extra))

        heights = vcat(heights, heights_extra)
        VV = vcat(VV, VV_extra)
        T = vcat(T, T_extra)
        RHs = vcat(RHs, RH_extra)
    end

    if addheight
        VV = deleteat!(VV, 2)
        T = deleteat!(T, 2)
        RHs = deleteat!(RHs, 2)
        heights = deleteat!(heights, 2)
    end

    return (
        heights=unique(heights),
        VELs=u"m/s".(VV),      # m/s
        TAs=u"°C".(T),         # deg C
        RHs=RHs,               # %
        QCONV=u"W/m^2"(QC),    # W
        USTAR=u"m/s"(USTAR)    # m/s
    )
end