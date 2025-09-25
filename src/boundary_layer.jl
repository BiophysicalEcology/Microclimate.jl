const DEFAULT_HEIGHTS = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0] .* u"m"

function allocate_profile(heights, reference_height)
    if minimum(heights) > reference_height
        throw(ArgumentError("`reference_height` is lower than minimum height in `heights`"))
        # addheight = true
        # newheights = Vector{eltype(heights)}(length(heights) + 1)
        # newheights[1] = 0.01u"m"
        # newheights[2:end] .= heights
        # heights = newheights
    elseif maximum(heights) >= reference_height
        throw(ArgumentError("Some values in `heights` are larger than `reference_height`"))
    end
    heights = filter(h -> h < reference_height, heights)
    NAIR = length(heights) + 1
    AIRDP = Vector{typeof(reference_height)}(undef, NAIR)
    AIRDP[1] = reference_height
    AIRDP[end:-1:2] .= u"m".(heights)
    VV = zeros(typeof(1.0u"m/minute"), NAIR) # output wind speeds
    T = zeros(typeof(0.0u"K"), NAIR) # output temperatures, need to do this otherwise get InexactError
    RHs = zeros(Float64, NAIR) # output relative humidities
    # heights_orig = copy(heights)
    # heights = fill(0.0u"m", NAIR + 1)
    # heights[end:-1:2] .= AIRDP

    return (; heights, reference_height, AIRDP, VV, T, RHs)
end

get_profile(; heights=DEFAULT_HEIGHTS, reference_height=1.2u"m", kw...) =
    get_profile!(allocate_profile(heights, reference_height); kw...)
function get_profile!(buffers;
    z0=0.004u"m",
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
    elevation=0.0u"m",
)
    (; heights, reference_height, AIRDP, VV, T, RHs) = buffers

    minimum(heights) < z0 && _minimum_heigth_error(heights, z0)

    NAIR = length(T)
    T1 = u"K"(TAREF)
    T3 = u"K"(D0cm)

    # Units: m to cm
    z = u"cm"(reference_height)
    z0 = u"cm"(z0)
    zh_cm = u"cm"(zh)
    d0_cm = u"cm"(d0)
    V = u"cm/minute"(VREF)
    # define air heights
    ZZ = AIRDP # TODO why rename
    VV[1] = V
    T[1] = T1

    # compute rcptkg (was a constant in original Fortran version)
    dry_air_out = dry_air_properties(u"K"(TAREF); elevation)
    wet_air_out = wet_air_properties(u"K"(TAREF); rh)
    ρ = dry_air_out.ρ_air
    c_p = wet_air_out.c_p
    TREF = u"K"(TAREF)
    rcptkg = u"cal*minute^2/cm^4"(ρ * c_p * TREF / (κ * g_n))
    #rcptkg = 6.003e-8u"cal*minute^2/cm^4"
    GAM = 16.0
    ZRATIO = z / z0 + 1.0
    DUM = log(ZRATIO)
    USTAR = κ * V / DUM
    DIFFT = T1 - T3
    TAVE = (T3 + T1) / 2
    RCP = RHOCP(TAVE)#, elevation, rh)
    AMOL = -30.0u"cm"
    if zh > 0.0u"m"
        # TODO ustrip to what
        STS = 0.62 / (ustrip(z0) * ustrip(USTAR) / 12)^(9//20)
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
            # TODO ustrip to what
            STS = 0.62 / (ustrip(z0) * ustrip(USTAR) / 12.0)^(9//20)
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
    e = wet_air(T1; rh).P_vap
    RHs .= clamp.(e ./ vapour_pressure.(T) .* 100.0, 0.0, 100.0)

    return (;
        # heights=unique(heights), # There should be not duplicates at this stage?
        heights,
        wind_speeds=VV,      # m/s
        air_temperatures=T,         # deg C
        humidities=RHs,               # %
        qconv=u"W/m^2"(QC),    # W
        ustar=u"m/s"(USTAR)    # m/s
    )
end

@noinline _minimum_heigth_error(heights, z0) =
    error("The minimum height $(minimum(heights)) is not greater than the roughness height ($z0).")

function RHOCP(TAVE)
    # TODO ustrip to what
    return u"(cal*g)/(g*cm^3*K)" * (0.08472 / ustrip(TAVE))
end
function RHOCP(TAVE, elevation, rh)
    dry_air_out = dry_air_properties(u"K"(TAVE); elevation)
    wet_air_out = wet_air_properties(u"K"(TAVE); rh)
    ρ = dry_air_out.ρ_air
    c_p = wet_air_out.c_p
    return u"(cal*g)/(g*cm^3*K)"(ρ * c_p)
end

function PHI(z, GAM, AMOL)
    # TODO ustrip to what
    return (1.0 - min(1.0, GAM * ustrip(z / AMOL)))^(1//4)
end

function PSI1(X)
    return 2.0 * log((1.0 + X) / 2.0) + log((1.0 + X^2) / 2.0) - 2.0 * atan(X) + π / 2.0
end

function PSI2(X)
    return 2.0 * log((1 + X^2.0) / 2.0)
end

function get_Obukhov(TA, TS, V, z, z0, rcptkg, κ)
    AMOL = -30.0u"cm" # initial Monin-Obukhov length cm
    GAM = 16.0 # -
    #RCPTKG = 6.003e-8u"cal/minute/cm/K" #CAL-MIN-CM-C
    z = u"cm"(z)
    z0 = u"cm"(z0)
    ZRATIO = z / z0 + 1.0
    DUM = log(ZRATIO)
    TA = u"K"(TA)
    TS = u"K"(TS)
    DIFFT = TA - TS
    TAVE = (TA + TS) / 2.0
    RCP = RHOCP(TAVE)
    DEL = 1.0
    count = 0
    USTAR = (κ * V / DUM)u"cm/minute"
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
            # TODO ustrip to what
            STS = 0.62 / (ustrip(z0) * ustrip(USTAR) / 12)^(9//20)
            STB = 0.64 / DUM
            QC = RCP * DIFFT * USTAR * STB / (1.0 + STB / STS)
        else
            # TODO ustrip to what
            STS = 0.62 / (ustrip(z0) * ustrip(USTAR) / 12)^(9//20)
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

