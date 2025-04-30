"""
    hour_angle(t::Quantity, lonc::Quantity) -> Quantity

Compute the solar hour angle `h` in radians.

# Arguments
- `t`: Local solar hour (e.g., `14.0`)
- `lonc`: Longitude correction in hours (e.g., `0.5`)

# Returns
- Hour angle `h` as a `Quantity` in radians
- Time at solar noon, `tsn` as a time in hours

# Reference
McCullough & Porter 1971, Eq. 6
"""
function hour_angle(t::Real, lonc::Real=0)
    tsn = 12.0 + lonc                      # solar noon time
    h = (π / 12) * (t - tsn) * u"rad"      # convert hours to radians
    return h, tsn
end

"""
    solar_geometry(d::Real, lat::Quantity, h::Quantity; d0::Real = 80, ω::Real = 2π/365, ϵ::Real = 0.0167, se::Real = 0.39779)

Computes key solar geometry parameters based on McCullough & Porter (1971):

- `ζ`: Auxiliary solar longitude (radians)
- `δ`: Solar declination (radians)
- `z`: Solar zenith angle (radians)
- `AR2`: Square of Earth-to-Sun radius factor (unitless)

# Arguments
- `d`: Day of year (1–365)
- `lat`: Latitude (with angle units, e.g. `u"°"` or `u"rad"`)
- `h`: Hour angle (radians)
- `d0`: Reference day (default: 80)
- `ω`: Angular frequency of Earth’s orbit (default: `2π/365`)
- `ϵ`: Orbital eccentricity (default: `0.0167`)
- `se`: Constant for solar declination amplitude (default: `0.39779`)

# Returns
Tuple: `(ζ, δ, z, AR2)` with angle quantities in radians and AR2 unitless.

# Reference
McCullough & Porter (1971)
"""
function solar_geometry(;
    d::Real=1.0,
    lat::Quantity=83.07305u"°",
    h::Quantity=-2.87979u"rad",
    d0::Real=80,
    ω::Real=2π / 365,
    ϵ::Real=0.0167238,
    se::Real=0.39779
)
    ζ = (ω * (d - d0)) + 2ϵ * (sin(ω * d) - sin(ω * d0))          # Eq.5
    δ = asin(se * sin(ζ))                                         # Eq.4
    cosZ = cos(lat) * cos(δ) * cos(h) + sin(lat) * sin(δ)         # Eq.3
    z = acos(cosZ)u"rad"                                          # Zenith angle
    AR2 = 1 + (2ϵ) * cos(ω * d)                                   # Eq.2
    δ = δ * u"rad"
    ζ = ζ * u"rad"
    return ζ, δ, z, AR2
end

"""
    check_skylight(z, nmax, SRINT, GRINT)

Checks for possible skylight before sunrise or after sunset based on zenith angle.
Modifies SRINT and GRINT at index `nmax` if skylight is present.

# Arguments
- `z::Quantity`: Zenith angle
- `nmax::Int`: Index into result arrays
- `SRINT::Vector{Quantity}`: Scattered radiation array [W/m²]
- `GRINT::Vector{Quantity}`: Global radiation array [W/m²]
"""
function check_skylight(
    z::Quantity,
    nmax::Int,
    SRINT::Vector,
    GRINT::Vector)
    Z = uconvert(°, z).val # convert to degrees
    if Z < 107.0
        if Z > 88.0
            Elog = 41.34615384 - 0.423076923 * Z
            Skylum = (10.0^Elog) * 1.46E-03u"mW * cm^-2"
            SRINT[nmax] = Skylum
            GRINT[nmax] = SRINT[nmax]
        end
    end

    return nothing
end

"""
    elev_corr(elev)

Calculates smooth polynomial approximations of atmospheric constituent correction factors 
as a function of altitude (based on Kearney's modification of the ALTFCT array originally 
from SOLAR.DAT). Input `elev` is the altitude in meters and can include units.

# Description

The array `ELEVFCT(i, j)` represents the **ratio of the total amount of a given 
atmospheric constituent (index j) above the elevation of interest (index i) to that 
above sea level**. The constituent indices are:

- j = 1: Molecular
- j = 2: Aerosol
- j = 3: Ozone
- j = 4: Water vapor

For j = 1–3, values are derived from standard profiles. For water vapor (j = 4), no 
standard profile exists, so only `ELEVFCT(1, 4)` is defined as 1.00.

The elevation index i runs from 1 to 21, corresponding to elevations from sea level to 
20 km in 1 km steps.

This function implements fitted polynomials to reproduce this correction smoothly 
from `elev` (in meters) using continuous approximation.

# Returns

A named tuple with the following keys:
- `ELEVFCT1` for Molecular
- `ELEVFCT2` for Aerosol
- `ELEVFCT3` for Ozone
- `ELEVFCT4` for Water vapor
"""
function elev_corr(elev)
    # Strip units if present and convert to km, then add 1
    elev_km = ustrip(u"m", elev) / 1000 + 1

    ELEVFCT1 = 0.00007277 * elev_km^3 +
               0.00507293 * elev_km^2 -
               0.12482149 * elev_km +
               1.11687469

    ELEVFCT2 = 8.35656e-7 * elev_km^6 -
               6.26384e-5 * elev_km^5 +
               1.86967e-3 * elev_km^4 -
               2.82585e-2 * elev_km^3 +
               2.26739e-1 * elev_km^2 -
               9.25268e-1 * elev_km +
               1.71321

    ELEVFCT3 = 1.07573e-6 * elev_km^5 -
               5.14511e-5 * elev_km^4 +
               7.97960e-4 * elev_km^3 -
               4.90904e-3 * elev_km^2 +
               2.99258e-3 * elev_km +
               1.00238

    ELEVFCT4 = 1.0

    return (ELEVFCT1, ELEVFCT2, ELEVFCT3, ELEVFCT4)
end

function GAMMA(TAU1::Float64)

    CHX = zeros(Float64, 101)
    CHY = zeros(Float64, 101)
    CFA = zeros(Float64, 3)
    AMU = zeros(Float64, 101)
    X1 = zeros(Float64, 101)
    Y1 = zeros(Float64, 101)
    X2 = zeros(Float64, 101)
    Y2 = zeros(Float64, 101)
    AIL = zeros(Float64, 101)
    AI = zeros(Float64, 30)
    XA = zeros(Float64, 4)
    XB = zeros(Float64, 8)
    GAMR = zeros(Float64, 101)
    GAML = zeros(Float64, 101)
    # Set up AMU array
    AMU[1] = 0.0
    for I in 2:101
        AMU[I] = 0.01 * (I - 1)
    end

    # Compute X1, Y1 using dchxy
    CFA[1] = 0.75
    CFA[2] = -0.75
    CFA[3] = 0.0
    NST = 111
    CHX, CHY, NTR = dchxy(TAU1, CFA, NST)
    for I in 1:101
        X1[I] = CHX[I]
        Y1[I] = CHY[I]
    end

    # Compute X2, Y2 using dchxy
    CFA[1] = 0.375
    CFA[2] = -0.375
    NST = 0
    CHX, CHY, NTR = dchxy(TAU1, CFA, NST)
    for I in 1:101
        X2[I] = CHX[I]
        Y2[I] = CHY[I]
    end

    # Compute AIL
    AIL[1] = 0.01 / 3.0
    CNU1 = 4.0 * AIL[1]
    CNU2 = 2.0 * AIL[1]
    for I in 2:2:100
        AIL[I] = CNU1
        AIL[I+1] = CNU2
    end
    AIL[101] = AIL[1]

    # Initialize integrals
    fill!(XA, 0.0)
    fill!(XB, 0.0)

    for I in 1:101
        c1 = AIL[I] * X1[I] * AMU[I]
        XA[1] += c1
        XA[2] += c1 * AMU[I]
        c2 = AIL[I] * Y1[I] * AMU[I]
        XA[3] += c2
        XA[4] += c2 * AMU[I]
        c3 = AIL[I] * X2[I]
        XB[1] += c3
        XB[2] += c3 * AMU[I]
        XB[3] += c3 * AMU[I]^2
        XB[4] += c3 * AMU[I]^3
        c4 = AIL[I] * Y2[I]
        XB[5] += c4
        XB[6] += c4 * AMU[I]
        XB[7] += c4 * AMU[I]^2
        XB[8] += c4 * AMU[I]^3
    end

    AI[1] = XB[1] + XB[5] - 8.0 / 3.0
    AI[2] = XB[2] + XB[6]
    AI[3] = XB[3] + XB[7]
    AI[4] = XB[1] - XB[5] - 8.0 / 3.0
    AI[5] = XB[2] - XB[6]
    AI[6] = XB[3] - XB[7]
    AI[7] = XB[4] - XB[8]
    AI[8] = XA[1] + XA[3]
    AI[9] = XA[2] + XA[4]
    AI[10] = XA[1] - XA[3]
    AI[11] = XA[2] - XA[4]

    AI[12] = (AI[1] - AI[3]) / ((AI[4] - AI[6]) * TAU1 + 2.0 * (AI[5] - AI[7]))
    AI[13] = 1.0 / (AI[4] * AI[10] - AI[5] * AI[11])
    AI[14] = 1.0 / (AI[1] * AI[8] - AI[2] * AI[9] - 2.0 * AI[12] * (AI[5] * AI[8] - AI[4] * AI[9]))
    AI[15] = 2.0 * (AI[8] * AI[10] - AI[9] * AI[11])
    AI[16] = AI[13] * AI[15]
    AI[17] = AI[14] * AI[15]

    CNU1 = 0.5 * (AI[16] - AI[17])
    CNU2 = 0.5 * (AI[16] + AI[17])

    AI[15] = AI[13] * (AI[5] * AI[8] - AI[4] * AI[9])
    AI[16] = AI[14] * (AI[2] * AI[10] - AI[1] * AI[11] - 2.0 * AI[12] * (AI[4] * AI[10] - AI[5] * AI[11]))
    CNU3 = 0.5 * (AI[15] - AI[16])
    CNU4 = 0.5 * (AI[15] + AI[16])

    AI[15] = AI[13] * (AI[2] * AI[10] - AI[1] * AI[11])
    AI[16] = AI[14] * (AI[5] * AI[8] - AI[4] * AI[9])
    CU3 = 0.5 * (AI[15] - AI[16])
    CU4 = 0.5 * (AI[15] + AI[16])

    AI[15] = AI[14] * (AI[1] * AI[8] - AI[2] * AI[9])
    SBAR = 1.0 - 0.375 * AI[12] * (AI[4] - AI[6]) *
                 ((CNU2 - CNU1) * AI[8] + (CU4 - CU3) * AI[2] - AI[15] * AI[6])

    AI[20] = 0.375 * AI[12] * (CNU2 - CNU1) * (AI[4] - AI[6])
    AI[21] = 0.375 * AI[12] * (AI[4] - AI[6])
    AI[22] = AI[21] * (CU4 - CU3)
    AI[23] = AI[21] * AI[15]

    for I in 1:101
        GAML[I] = AI[20] * (X1[I] + Y1[I])
        GAMR[I] = AI[22] * (X2[I] + Y2[I]) - AMU[I] * AI[23] * (X2[I] - Y2[I])
    end
    return GAMR, GAML, SBAR
end

# function dchxy(tau1::Float64, cfa::Vector{Float64}, nst::Int)
#     chx = zeros(Float64, 101)
#     chy = zeros(Float64, 101)

#     if tau1 < 0.0
#         return chx, chy, 1
#     end

#     c = zeros(Float64, 3)
#     u = zeros(Float64, 3)
#     c[1] = cfa[1]
#     c[2] = cfa[2]
#     c[3] = cfa[3]

#     # Calculate initial values for root finding
#     cc = 0.0
#     cd = 0.0
#     for n in 1:3
#         cc += (2n - 1) * c[n]
#         cd += (2n - 1)^2 * c[n]
#     end

#     # Compute roots of characteristic equation (up to 4 roots)
#     q = zeros(ComplexF64, 4)
#     d = zeros(ComplexF64, 4)
#     ntr = 2
#     q[1] = sqrt(Complex(0.25 * cc + 0.5 * sqrt(Complex(cd))))
#     q[2] = sqrt(Complex(0.25 * cc - 0.5 * sqrt(Complex(cd))))
#     d[1] = 1.0 / (2.0 * q[1] * (2.0 * q[1]^2 - cc))
#     d[2] = 1.0 / (2.0 * q[2] * (2.0 * q[2]^2 - cc))

#     # Compute exponential terms using dexpi
#     e1 = dexpi(-real(q[1]) * tau1)
#     e2 = dexpi(-real(q[2]) * tau1)

#     # Initialize mu array (Gauss quadrature nodes)
#     mu = zeros(Float64, 101)
#     mu[1] = 0.0
#     for i in 2:101
#         mu[i] = 0.01 * (i - 1)
#     end

#     # Compute chx and chy using X and Y function definitions
#     for i in 1:101
#         sumx = 0.0 + 0im
#         sumy = 0.0 + 0im
#         for j in 1:ntr
#             qq = q[j]
#             dj = d[j]
#             sumx += dj * mu[i] / (mu[i]^2 - qq^2)
#             sumy += dj * qq / (mu[i]^2 - qq^2)
#         end
#         chx[i] = real(1.0 + 2.0 * mu[i] * sumx)
#         chy[i] = real(2.0 * sumy)
#     end

#     return chx, chy, ntr
# end

function dchxy(TAU1::Float64, CFA::Vector{Float64}, NCASE::Int)
    # Format placeholders (not used directly here, but shown for reference)
    # 5 FORMAT(//T5,'THE PROGRAM IS TERMINATED AS PSI(I) =',D15.6,T50,' FOR I = ',I5//)
    # 20 FORMAT(//T10,'NO COMPUTATIONS CAN BE DONE AS THE COEFFIECIENTS ARE = ',3F15.5//)
    # 30 FORMAT(//T5,' THE PROGRAM IS TERMINATED BECAUSE TAU1 = ',D15.5//)

    # Initialize arrays
    PSI = zeros(Float64, 101)
    AMU = zeros(Float64, 101)
    XA = zeros(Float64, 101)
    XB = zeros(Float64, 101)
    UMA = zeros(Float64, 5)
    ACAP = zeros(Float64, 5)
    TEMX = zeros(Float64, 8)
    TEMY = zeros(Float64, 8)
    RTK = zeros(Float64, 5)
    ALAM = zeros(Float64, 5)
    FNPP = zeros(Float64, 101)
    FNPN = zeros(Float64, 101)
    FNC0 = zeros(Float64, 101)
    FNC1 = zeros(Float64, 101)
    FNX = zeros(Float64, 101)
    FNY = zeros(Float64, 101)
    FNW = zeros(Float64, 101)
    FMC0 = zeros(Float64, 101)
    FMC1 = zeros(Float64, 101)
    XC = zeros(Float64, 101)   # equivalence
    XD = zeros(Float64, 101)   # equivalence
    XE = zeros(Float64, 101)   # equivalence
    CHXA = zeros(Float64, 101)   # equivalence
    CHYA = zeros(Float64, 101)    # equivalence
    CHX = zeros(Float64, 101)
    CHY = zeros(Float64, 101)
    XA = zeros(Float64, 101)
    XB = zeros(Float64, 101)
    # Integer arrays
    NC0 = reshape(Int[
            3, 4, 1, 2,
            2, 4, 1, 3,
            2, 3, 1, 4,
            1, 4, 2, 3,
            1, 3, 2, 4,
            1, 2, 3, 4,
        ], 4, 6)
    NC1 = reshape(Int[
            2, 3, 4, 1,
            1, 3, 4, 2,
            1, 2, 4, 3,
            1, 2, 3, 4,
            4, 1, 2, 3,
            3, 1, 2, 4,
            2, 1, 3, 4,
            1, 2, 3, 4,
        ], 4, 8)

    # Variables
    PERA = 0.0
    PERB = 0.0

    # Terminate if TAU1 is too large or negative
    if TAU1 <= 2.0
        # proceed
    else
        println(" THE PROGRAM IS TERMINATED BECAUSE TAU1 = ", TAU1)
        error("Program terminated due to TAU1 > 2.0")
    end

    if TAU1 < 0.0
        println(" THE PROGRAM IS TERMINATED BECAUSE TAU1 = ", TAU1)
        error("Program terminated due to TAU1 < 0.0")
    end

    PERA = CFA[1] + CFA[2] / 3.0 + 0.2 * CFA[3]
    if NCASE != 0
        PERA = 0.5
    end

    if PERA < 0.0
        println("NO COMPUTATIONS CAN BE DONE AS THE COEFFICIENTS ARE = ", CFA[1], " ", CFA[2], " ", CFA[3])
        error("Program terminated due to PERA < 0.0")
    end

    if PERA > 0.5
        println("NO COMPUTATIONS CAN BE DONE AS THE COEFFICIENTS ARE = ", CFA[1], " ", CFA[2], " ", CFA[3])
        error("Program terminated due to PERA > 0.5")
    end

    # Compute MU, PSI(MU), and weights
    for i in 1:101
        AMU[i] = (i - 1) * 0.01
        TEMA = AMU[i]^2
        PSI[i] = CFA[1] + CFA[2] * TEMA + CFA[3] * (TEMA^2)
        if PSI[i] > -1e-15
            continue
        else
            println("THE PROGRAM IS TERMINATED AS PSI($i) = ", PSI[i])
            println("NO COMPUTATIONS CAN BE DONE AS THE COEFFICIENTS ARE = ", CFA[1], " ", CFA[2], " ", CFA[3])
            error("Program terminated due to PSI(i) < threshold")
        end
    end

    XA[1] = 0.01 / 3.0
    TEMA = 4.0 * XA[1]
    TEMB = 2.0 * XA[1]
    for i in 2:2:100
        XA[i] = TEMA
        XA[i+1] = TEMB
    end
    XA[101] = XA[1]

    # Suppress all intermediate output
    NPRT = 0

    # Compute roots of the characteristic equation
    if NCASE != 0
        KMX = 5
        UMA[1] = 0.97390652851717172
        UMA[2] = 0.86506336668898451
        UMA[3] = 0.67940956829902441
        UMA[4] = 0.43339539412924719
        UMA[5] = 0.14887433898163121

        ACAP[1] = 0.066671344308688138
        ACAP[2] = 0.14945134915058059
        ACAP[3] = 0.21908636251598204
        ACAP[4] = 0.26926671930999636
        ACAP[5] = 0.29552422471475287
    else
        KMX = 4
        N1 = 0
        UMA[1] = 0.96028985649753623
        UMA[2] = 0.79666647741362674
        UMA[3] = 0.52553240991632899
        UMA[4] = 0.18343464249564980

        ACAP[1] = 0.10122853629037626
        ACAP[2] = 0.22238103445337447
        ACAP[3] = 0.31370664587788729
        ACAP[4] = 0.36268378337836198
    end

    for i in 1:KMX
        TEMX[i] = UMA[i]^2
        TEMY[i] = CFA[1] + CFA[2] * TEMX[i] + CFA[3] * TEMX[i]^2
        TEMY[i] = 2.0 * ACAP[i] * TEMY[i]
    end

    if NCASE != 0
        IST = 2
        RTK[1] = 0.0
    else
        IST = 1
    end

    for i in IST:KMX # Fortran line 152
        RTK[i] = (1.0 - TEMY[i]) / TEMX[i]
        if i == 1
            TEMA = 1.0 / UMA[1]^2
            if RTK[1] >= TEMA
                RTK[1] = 0.5 * TEMA
            end
        else
            TEMA = 1.0 / UMA[i-1]^2
            TEMB = 1.0 / UMA[i]^2
            if !(RTK[i] > TEMA && RTK[i] < TEMB)
                RTK[i] = 0.5 * (TEMA + TEMB)
            end
        end
    end

    J = IST # Fortran line 164
    while J <= KMX
        if J == 1
            TEMA = 0.0
            TEMB = 1.0 / UMA[1]^2
            N1 = 0
        else
            TEMA = 1.0 / UMA[J-1]^2
            TEMB = 1.0 / UMA[J]^2
            N1 = 0
        end

        TEMC = 1.0
        for i in 1:KMX
            TEMC -= TEMY[i] / (1.0 - RTK[J] * TEMX[i])
        end
        TEMD = abs(TEMC)
        if TEMD < 1e-14
            J += 1
        else
            N1 += 1
            if N1 > 50
                println("THE PROGRAM IS TERMINATED BECAUSE ROOTS CANNOT BE FOUND TO SATISFY THE CRITERION..")
                println("CFA(1) = $(CFA[1])    CFA(2) = $(CFA[2])    CFA(3) = $(CFA[3])")
                println("THE TROUBLE IS WITH ROOT NUMBER ", J)
                error("Root finding failed.")
            end

            if TEMC > 0.0
                TEMA = RTK[J]
            elseif TEMC < 0.0
                TEMB = RTK[J]
            end

            TEMD = 0.0
            for i in 1:KMX
                TEMD -= (TEMY[i] * TEMX[i]) / (1.0 - RTK[J] * TEMX[i])^2
            end

            TEMC = RTK[J] - TEMC / TEMD

            if TEMC <= TEMA || TEMC >= TEMB
                RTK[J] = 0.5 * (TEMA + TEMB)
            else
                RTK[J] = TEMC
            end
        end
    end

    for i in 1:KMX
        RTK[i] = sqrt(RTK[i])
    end

    if NCASE != 0
        N1 = 11
        KMX = 4
        for j in 1:KMX
            RTK[j] = RTK[j+1]
        end
    end
    UMA[1] = 0.96028985649753623
    UMA[2] = 0.79666647741362674
    UMA[3] = 0.52553240991632899
    UMA[4] = 0.18343464249564980

    ACAP[1] = 0.10122853629037626
    ACAP[2] = 0.22238103445337447
    ACAP[3] = 0.31370664587788729
    ACAP[4] = 0.36268378337836198

    # --- COMPUTE FUNCTIONS LAMDA, P AND W ---
    for j in 1:KMX
        ALAM[j] = 1.0
        for i in 1:KMX
            ALAM[j] *= (RTK[j] * UMA[i] + 1.0) / (RTK[j] * UMA[i] - 1.0)
        end
        ALAM[j] = exp(-RTK[j] * TAU1) / ALAM[j]
    end

    if NPRT != 0
        Printf.printf("%12.5E %12.5E %12.5E\n", CFA[1], CFA[2], CFA[3])
        Printf.printf("%12.5E\n", TAU1)
        Printf.printf("\n")
        for j in 1:KMX
            TEMA = 1.0 / RTK[j]
            # (In the FORTRAN code, TEMA is calculated but not used or printed here)
        end
    end

    for i in 1:101 # Fortran line 225
        FNPP[i] = 1.0
        FNPN[i] = 1.0
        FNW[i] = 1.0
        for j in 1:KMX
            FNPP[i] *= (AMU[i] / UMA[j] - 1.0)
            FNPN[i] *= (-AMU[i] / UMA[j] - 1.0)
            FNW[i] *= (1.0 - RTK[j]^2 * AMU[i]^2)
        end
    end
    # --- COMPUTE C₀ AND C₁ ---

    TEMX[1] = 1.0
    TEMX[8] = 1.0
    for k in 2:7
        TEMX[k] = 1.0
        for i in 1:2
            N1 = NC0[i, k-1]
            for j in 1:2
                N2 = NC0[j+2, k-1]
                TEMX[k] *= (RTK[N1] + RTK[N2]) / (RTK[N1] - RTK[N2])
            end
        end
        TEMX[k] = -TEMX[k]
    end

    for k in 1:4
        TEMY[k] = 1.0
        N2 = NC1[4, k]
        for i in 1:3
            N1 = NC1[i, k]
            TEMY[k] *= (RTK[N1] + RTK[N2]) / (RTK[N1] - RTK[N2])
        end
    end

    for k in 5:8
        TEMY[k] = 1.0
        N1 = NC1[1, k]
        for j in 1:3
            N2 = NC1[j+1, k]
            TEMY[k] *= (RTK[N1] + RTK[N2]) / (RTK[N1] - RTK[N2])
        end
        TEMY[k] = -TEMY[k]
    end

    for i in 1:101 # Fortran line 266
        TEMA = 1.0
        TEMB = 1.0
        for j in 1:4
            TEMA *= (1.0 + RTK[j] * AMU[i])
            TEMB *= (1.0 - RTK[j] * AMU[i])
        end
        FNC0[i] = TEMA
        FMC0[i] = TEMB

        TEMA = 1.0
        TEMB = 1.0
        for j in 1:4
            TEMA *= (1.0 - RTK[j] * AMU[i]) * ALAM[j]
            TEMB *= (1.0 + RTK[j] * AMU[i]) * ALAM[j]
        end
        FNC0[i] += TEMA
        FMC0[i] += TEMB

        IST = 2
        while IST <= 7
            TEMA = 1.0
            TEMB = 1.0
            for k in 1:2
                N2 = NC0[k+2, IST-1]
                TEMA *= (1.0 - RTK[N2] * AMU[i]) * ALAM[N2]
                TEMB *= (1.0 + RTK[N2] * AMU[i]) * ALAM[N2]
            end
            for j in 1:2
                N1 = NC0[j, IST-1]
                TEMA *= (1.0 + RTK[N1] * AMU[i])
                TEMB *= (1.0 - RTK[N1] * AMU[i])
            end
            FNC0[i] += TEMA * TEMX[IST]
            FMC0[i] += TEMB * TEMX[IST]

            IST += 1
        end
    end
    for i in 1:101
        FNC1[i] = 0.0
        FMC1[i] = 0.0
        IST = 1
        while IST <= 4
            N2 = NC1[4, IST]
            TEMA = (1.0 - RTK[N2] * AMU[i]) * ALAM[N2]
            TEMB = (1.0 + RTK[N2] * AMU[i]) * ALAM[N2]
            for j in 1:3
                N1 = NC1[j, IST]
                TEMA *= (1.0 + RTK[N1] * AMU[i])
                TEMB *= (1.0 - RTK[N1] * AMU[i])
            end
            FNC1[i] += TEMY[IST] * TEMA
            FMC1[i] += TEMY[IST] * TEMB
            IST += 1
        end
        while IST <= 8
            N1 = NC1[1, IST]
            TEMA = 1.0 + RTK[N1] * AMU[i]
            TEMB = 1.0 - RTK[N1] * AMU[i]
            for j in 1:3
                N2 = NC1[j+1, IST]
                TEMA *= (1.0 - RTK[N2] * AMU[i]) * ALAM[N2]
                TEMB *= (1.0 + RTK[N2] * AMU[i]) * ALAM[N2]
            end
            FNC1[i] += TEMY[IST] * TEMA
            FMC1[i] += TEMY[IST] * TEMB
            IST += 1
        end
        FNC1[i] = -FNC1[i]
        FMC1[i] = -FMC1[i]
    end

    if NPRT != 0
        # These WRITE statements are formatted outputs; replacing with println for now.
        println()
        println("CFA values: ", CFA[1:3])
        println("TAU1 = ", TAU1)
        println("Table header for output:")

        for i in 1:101
            TEMD = FNC0[i] * FMC0[i] - FNC1[i] * FMC1[i] - (FNC0[1]^2 - FNC1[1]^2) * FNW[i]
            println(AMU[i], " ", FNPP[i], " ", FNPN[i], " ", FNW[i], " ", FNC0[i], " ", FMC0[i], " ", FNC1[i], " ", FMC1[i], " ", TEMD)
        end

        println()
    end
    # COMPUTE THE FOURTH APPROXIMATION OF X AND Y FUNCTIONS
    XB[1] = 0.0 # Fortran line 345
    if TAU1 == 0.0
        XB[1] = 1.0
    end

    for i in 2:101
        XB[i] = exp(-TAU1 / AMU[i])
    end

    TEMA = 1.0 / sqrt(FNC0[1]^2 - FNC1[1]^2)
    for i in 1:101
        TEMC = TEMA / FNW[i]
        FNX[i] = (FNPN[i] * FMC0[i] - XB[i] * FNPP[i] * FNC1[i]) * TEMC
        FNY[i] = (XB[i] * FNPP[i] * FNC0[i] - FNPN[i] * FMC1[i]) * TEMC
    end

    CHXA[1] = 1.0
    CHYA[1] = XB[1]
    nomitr = 1 # Fortran line 362
    converged = false
while converged == false
    for I in 2:101
        for IC in 1:101
            XD[IC] = PSI[IC] * (FNX[I] * FNX[IC] - FNY[I] * FNY[IC]) / (AMU[I] + AMU[IC])
            if I != IC
                XE[IC] = PSI[IC] * (FNY[I] * FNX[IC] - FNX[I] * FNY[IC]) / (AMU[I] - AMU[IC])
            end
            if I <= 3
                XE[I] = 0.5 * (XE[I+1] + XE[I-1])
            else
                if I > 3 && I <= 5
                    # Everett's formula two points on either side
                    XE[I] = 0.0625 * (9.0 * (XE[I+1] + XE[I-1]) - XE[I+3] - XE[I-3])
                else
                    if I > 5 && I <= 96
                        # Everett's formula three points on either side
                        XE[I] = 3.0 * (XE[I+5] + XE[I-5]) + 150.0 * (XE[I+1] + XE[I-1]) - 25.0 * (XE[I+3] + XE[I-3])
                        XE[I] /= 256.0
                    else
                        # Interpolation for I > 96
                        XE[I] = 5.0 * XE[I-1] + 10.0 * XE[I-3] + XE[I-5] - 10.0 * XE[I-2] - 5.0 * XE[I-4]
                    end
                end
            end
            CHXA[I] = 0.0
            CHYA[I] = 0.0
            for IC in 1:101
                CHXA[I] += XA[IC] * XD[IC]
                CHYA[I] += XA[IC] * XE[IC]
            end
            CHXA[I] = 1.0 + AMU[I] * CHXA[I]
            CHYA[I] = XB[I] + AMU[I] * CHYA[I]
        end
    end
    # correction to the approximation

    if nomitr == 1 && TAU1 != 0 # Fortran line 398
        TEMX[1] = -dexpi(-TAU1)
        for n in 2:7
            TEMX[n] = (XB[101] - TAU1 * TEMX[n-1]) / (n - 1)
        end
        PERB = 2.0 * (
            CFA[1] * (0.5 - TEMX[3]) +
            CFA[2] * (0.25 - TEMX[5]) +
            CFA[3] * ((1 / 6) - TEMX[7])
        )
    end
    # accumulate TEMA, TEMB
    if TAU1 != 0.0
        TEMA = 0.0
        TEMB = 0.0
        for i in 1:101
            TEMA += CHXA[i] * PSI[i] * XA[i]  
            TEMB += CHYA[i] * PSI[i] * XA[i]  
        end
        # new TEMC
        c1 = (1 - 2 * PERA) / (1 - TEMA + TEMB)
        TEMC = TAU1 == 0 ? 0.0 : (1 - TEMA - TEMB - c1) / PERB
    end
    if TAU1 == 0.0
        TEMC = 0.0
    end
    for i in 1:101
        TEMD = TEMC * AMU[i] * (1.0 - XB[i])
        CHX[i] = CHXA[i] + TEMD
        CHY[i] = CHYA[i] + TEMD
    end
    if NPRT != 0 # Fortran line 422
        TEMC_out = TEMA^2 - 2.0*TEMA - TEMB^2 + 2.0*PERA
        #@printf(" NOMITR = %d, TEMA = %g, TEMB = %g, TEMC = %g\n", NOMITR, TEMA, TEMB, TEMC_out)
    end

    # Check for convergence
    if nomitr > 1
        #converged = true
        for I in 2:101
            rel_error = abs((CHY[I] - FNY[I]) / CHY[I])
            if rel_error <= 2.0e-4
                converged = true
                break
            end
        end
    end

    if converged
        break
    end

    # Prepare for next iteration
    for I in 1:101
        FNX[I] = CHX[I]
        FNY[I] = CHY[I]
    end
    nomitr += 1
    if nomitr > 15
        break
    end
end

    # if NCASE ≠ 0, generate standard solution (Fortran 975…990)
    if NCASE != 0
        tsumx = 0.0
        tsumb = 0.0
        tsumc = 0.0
        for i in 1:101
            δ = PSI[i] * AMU[i] * XA[i]
            tsumx += δ * CHX[i]
            tsumb += δ * CHY[i]
            tsumc += PSI[i] * CHY[i] * XA[i]
        end
        ratio = tsumc / (tsumx + tsumb)
        for i in 1:101
            Δ = ratio * AMU[i] * (CHX[i] + CHY[i])
            CHX[i] += Δ
            CHY[i] -= Δ
        end
    end

    return CHX, CHY, nomitr

end

function dexpi(x::Float64)
    # Constants
    gamma = 0.57721566490153286

    # Coefficients
    A1 = [0.1193095930415985, 0.3306046932331323, 0.4662347571015760]
    B1 = [0.4679139345726910, 0.3607615730481386, 0.1713244923791703]
    A2 = [0.02823912701457739, 30.52042817823067, 215.8885931211323,
        410.4611319636983, 278.5527592726121, 71.33086969436196, 0.5758931590224375]
    B2 = [10.0, 138.3869728490638, 488.08581830736, 634.8804630786363,
        344.1289899236299, 77.08964199043784, 0.5758934565014882]
    A3 = [0.07630772325814641, 21.23699219410890, 47.45350785776186,
        29.66421696379266, 6.444800036068992, 0.04295808082119383]
    B3 = [10.0, 52.78950049492932, 71.96111390658517, 35.67945294128107,
        6.874380519301884, 0.04295808112146861]
    A4 = [0.1157221173580207, 0.6117574845151307, 1.512610269776419,
        2.833751337743507, 4.599227639418348, 6.844525453115177,
        9.621316842456867, 13.00605499330635, 17.11685518746226,
        22.15109037939701, 28.48796725098400, 37.09912104446692]
    B4 = [0.2647313710554432, 0.3777592758731380, 0.2440820113198776,
        0.09044922221168093, 0.02010238115463410, 0.002663973541865316,
        0.0002032315926629994, 8.365055856819799e-5, 1.668493876540910e-6,
        1.342391030515004e-8, 3.061601635035021e-11, 8.148077467426242e-15]
    A5 = [0.03202844643130281, 0.09555943373680816, 0.1575213398480817,
        0.2168967538130226, 0.2727107356944198, 0.3240468259684878,
        0.3700620957892772, 0.4100009929869515, 0.4432077635022005,
        0.4691372760013664, 0.4873642779856547, 0.4975936099985107]
    B5 = [0.1279381953467522, 0.1258374563468283, 0.1216704729278034,
        0.1155056680537256, 0.1074442701159656, 0.09761865210411389,
        0.08619016153195328, 0.07334648141108031, 0.05929858491543678,
        0.04427743881741981, 0.02853138862893366, 0.01234122979998720]

    if x == 0.0
        error("The argument of DEXPI is very close to zero.")
    elseif x < 0.0
        ax = abs(x)
        if x > -1e-20
            return log(ax) + gamma
        elseif x > -1.5
            yy = exp(-0.5 * ax)
            s = 0.0
            for i in 1:3
                yz = exp(A1[i] * ax)
                s += B1[i] * ((1 - yy / yz) / (0.5 + A1[i]) + (1 - yy * yz) / (0.5 - A1[i]))
            end
            return -0.5 * s + log(ax) + gamma
        elseif x > -4.65
            sumn = evalpoly(ax, reverse(A2))
            sumd = evalpoly(ax, reverse(B2))
            return (sumn / (sumd * x)) * exp(x)
        elseif x > -12.0
            sumn = evalpoly(ax, reverse(A3))
            sumd = evalpoly(ax, reverse(B3))
            return (sumn / (sumd * x)) * exp(x)
        elseif x > -170.0
            dexpi = 0.0
            for j in 1:12
                dexpi += B4[j] / (1 + A4[j] / ax)
            end
            return (exp(x) / ax) * (-dexpi)
        else
            return 0.0
        end
    else
        if x <= 1e-20
            return log(x) + gamma
        elseif x <= 40.0
            yy = exp(0.5 * x)
            dexpi = 0.0
            for j in 1:12
                yz = exp(-A5[j] * x)
                dexpi += ((1 - yy / yz) / (0.5 + A5[j]) + (1 - yy * yz) / (0.5 - A5[j])) * B5[j]
            end
            return -0.5 * dexpi + log(x) + gamma
        elseif x <= 173.0
            dexpi = 0.0
            for j in 1:12
                dexpi += B4[j] / (1 - A4[j] / x)
            end
            return (exp(x) / x) * dexpi
        else
            error("The argument of DEXPI is very large.")
        end
    end
end

"""
    solrad(; days, hours, lat...[, year, lonc, elev, slope, aspect, hori, refl, cmH2O, ϵ,
           ω, se, d0, lamb, iuv, noscat, amr, nmax, Iλ, OZ, τR, τO, τA, τW, Sλ, FD, FDQ,
           S, ER, ERλ]) -> NamedTuple

Compute solar radiation at a given place and time using a detailed atmospheric radiative transfer model.

# Arguments
- `days::Vector{Float64}`: Days of the year (1–365/366) to evaluate.
- `hours::Vector{Float64}`: Decimal hours of the day (0.0–24.0).
- `lat::Quantity`: Latitude in degrees, e.g. `43.0u"°"`.

# Keyword Arguments
- `year::Real=2001`: Year used for ozone table lookup.
- `lonc::Real=0.0`: Longitude correction in hours (positive west of standard meridian).
- `elev::Quantity=0.0u"m"`: Elevation above sea level.
- `slope::Quantity=0u"°"`: Slope angle of the surface.
- `aspect::Quantity=0u"°"`: Azimuth of slope aspect (from north).
- `hori::Vector{Quantity}`: Horizon angles for each of 24 azimuth sectors (default 0°).
- `refl::Real=0.10`: Ground reflectance (albedo), fraction [0, 1].
- `cmH2O::Real=1`: Precipitable water in cm for atmospheric column (e.g. 0.1: dry, 1.0: moist, 2.0: humid).
- `ϵ::Real=0.0167238`: Orbital eccentricity of Earth.
- `ω::Real=2π/365`: Mean angular orbital velocity of Earth (radians/day).
- `se::Real=0.39779`: Precomputed solar elevation constant.
- `d0::Real=80`: Reference day for declination calculations.
- `lamb::Bool=false`: If `true`, returns wavelength-specific irradiance components.
- `iuv::Bool=false`: If `true`, uses the full gamma-function model for diffuse radiation (expensive).
- `noscat::Bool=true`: If `true`, disables scattered light computations (faster).
- `amr::Quantity=25.0u"km"`: Mixing ratio height of the atmosphere.
- `nmax::Integer=111`: Maximum number of wavelength intervals.
- `Iλ::Vector{Quantity}`: Vector of wavelength bins (e.g. in `nm`).
- `OZ::Matrix{Float64}`: Ozone column depth table indexed by latitude band and month (size 19×12).
- `τR`, `τO`, `τA`, `τW`: Vectors of optical depths per wavelength for Rayleigh scattering, ozone, aerosols, and water vapor.
- `Sλ::Vector{Quantity}`: Solar spectral irradiance per wavelength bin (e.g. in `mW * cm^-2 * nm^-1`).
- `FD`, `FDQ`: Auxiliary data vectors for biologically effective radiation models.

# Returns
A named tuple containing:
- `λ::Vector`: Wavelengths (typically in µm).
- `λDirect::Vector`: Spectral direct irradiance [W/m²/µm].
- `λRayleigh::Vector`: Spectral Rayleigh-scattered irradiance [W/m²/µm].
- `λScattered::Vector`: Spectral diffuse (scattered) irradiance [W/m²/µm].
- `λGlobal::Vector`: Spectral global irradiance (direct + diffuse) [W/m²/µm].

- `Direct::Float64`: Broadband direct irradiance [W/m²].
- `Rayleigh::Float64`: Broadband Rayleigh-scattered irradiance [W/m²].
- `Scattered::Float64`: Broadband scattered (diffuse) irradiance [W/m²].
- `Global::Float64`: Broadband global irradiance (direct + diffuse) [W/m²].

- `Zenith::Float64`: Solar zenith angle [degrees].
- `doy::Int`: Day of year.
- `hour::Float64`: Decimal hour of day.

# Notes
- Radiation units are returned in `W/m²`. Internally, units like `mW/cm²` are used and converted as necessary using `Unitful.jl`.
- Topographic shading is included via the `hori` input (horizon angle mask).
- Outputs are computed for each (day, hour) combination in the input vectors.
- The function can be run in a high-resolution spectral mode (`lamb=true`) or in a broadband mode (`lamb=false`).
- In optical air mass 'arims' calculation the difference between apparent and true zenith angle is neglected for z less than 88 degrees
- Variation of airms with altitude is ignored since it is negligible up to at least 6 km above sea level
"""
function solrad(;
    days::Vector{<:Real}=[15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349],
    hours::Vector{<:Real}=collect(0.0:24.0),
    year::Real=2001.0, # needed to determine if a leap year
    lat::Quantity=43.1379u"°",
    lonc::Real=0.0, # longitude correction, hours
    elev::Quantity=276.0u"m", # elevation, m
    slope::Quantity=0u"°",
    aspect::Quantity=0u"°",
    hori::Vector{typeof(0.0u"°")}=fill(0.0, 24) .* u"°",
    refl::Real=0.15, # substrate solar reflectivity (decimal %)
    cmH2O::Real=1, # precipitable cm H2O in air column, 0.1 = VERY DRY; 1 = MOIST AIR CONDITIONS; 2 = HUMID, TROPICAL CONDITIONS (note this is for the whole atmospheric profile, not just near the ground)
    ϵ::Real=0.0167238,
    ω::Real=2π / 365,
    se::Real=0.39779,
    d0::Real=80.0,
    iuv::Bool=false, # Use gamma function for scattered solar radiation? (computationally intensive)
    noscat::Bool=true,
    amr::Quantity=25.0u"km",
    nmax::Integer=111, # maximum number of wavelengths
    Iλ::Vector{typeof(0.0u"nm")}=float.([ # wavelengths across which to integrate
        290, 295, 300, 305, 310, 315, 320, 330, 340, 350, 360, 370, 380, 390,
        400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700,
        720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020,
        1080, 1100, 1120, 1140, 1160, 1180, 1200, 1220, 1240, 1260, 1280, 1300, 1320,
        1380, 1400, 1420, 1440, 1460, 1480, 1500, 1540, 1580, 1600, 1620, 1640, 1660,
        1700, 1720, 1780, 1800, 1860, 1900, 1950, 2000, 2020, 2050, 2100, 2120, 2150,
        2200, 2260, 2300, 2320, 2350, 2380, 2400, 2420, 2450, 2490, 2500, 2600, 2700,
        2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000
    ]) * u"nm",
    OZ::Matrix{<:Real}=reshape([
            0.31, 0.31, 0.32, 0.32, 0.31, 0.3, 0.27, 0.24, 0.23, 0.22, 0.23, 0.24,
            0.27, 0.3, 0.32, 0.33, 0.34, 0.34, 0.33, 0.3, 0.31, 0.31, 0.31, 0.3, 0.29, 0.28,
            0.25, 0.24, 0.22, 0.24, 0.26, 0.28, 0.32, 0.36, 0.39, 0.4, 0.4, 0.39, 0.3, 0.31,
            0.31, 0.3, 0.29, 0.28, 0.26, 0.24, 0.24, 0.23, 0.24, 0.26, 0.29, 0.33, 0.38, 0.42,
            0.45, 0.46, 0.46, 0.27, 0.28, 0.29, 0.3, 0.3, 0.29, 0.27, 0.25, 0.24, 0.23, 0.25,
            0.27, 0.3, 0.34, 0.38, 0.4, 0.42, 0.43, 0.42, 0.34, 0.35, 0.34, 0.33, 0.32, 0.31,
            0.28, 0.25, 0.24, 0.24, 0.26, 0.28, 0.3, 0.34, 0.37, 0.39, 0.4, 0.4, 0.39, 0.38,
            0.4, 0.39, 0.38, 0.36, 0.33, 0.28, 0.25, 0.24, 0.24, 0.25, 0.27, 0.3, 0.33, 0.35,
            0.36, 0.36, 0.36, 0.34, 0.43, 0.44, 0.43, 0.41, 0.39, 0.35, 0.29, 0.25, 0.24, 0.24,
            0.25, 0.26, 0.29, 0.31, 0.33, 0.34, 0.34, 0.33, 0.32, 0.45, 0.46, 0.45, 0.42, 0.4,
            0.37, 0.31, 0.26, 0.24, 0.23, 0.24, 0.26, 0.28, 0.3, 0.31, 0.32, 0.31, 0.3, 0.3,
            0.41, 0.42, 0.43, 0.42, 0.4, 0.38, 0.32, 0.26, 0.24, 0.23, 0.24, 0.26, 0.27, 0.28,
            0.3, 0.3, 0.29, 0.28, 0.27, 0.37, 0.38, 0.4, 0.4, 0.39, 0.37, 0.32, 0.26, 0.24,
            0.22, 0.23, 0.25, 0.26, 0.27, 0.28, 0.28, 0.28, 0.27, 0.26, 0.34, 0.36, 0.38, 0.39,
            0.37, 0.34, 0.29, 0.26, 0.24, 0.22, 0.23, 0.25, 0.26, 0.28, 0.29, 0.3, 0.29, 0.29,
            0.28, 0.31, 0.32, 0.34, 0.35, 0.35, 0.32, 0.29, 0.25, 0.23, 0.22, 0.23, 0.25, 0.27,
            0.29, 0.3, 0.31, 0.31, 0.31, 0.3
        ], (19, 12)),
    τR::Vector{<:Real}=[
        1.41, 1.31, 1.23, 1.14, 1.05, 0.99, 0.92, 0.81, 0.72, 0.63, 0.56, 0.5,
        0.45, 0.4, 0.36, 0.3, 0.25, 0.2, 0.17, 0.15, 0.12, 0.1, 0.09, 0.08, 0.07, 0.06,
        0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01,
        0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ],
    τO::Vector{<:Real}=[
        11.5, 6.3, 3.2, 1.62, 0.83, 0.44, 0.26, 0.03, 0.02, 0.01, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ],
    τA::Vector{<:Real}=[0.269904738, 0.266147825, 0.262442906, 0.258789404, 0.255186744, 0.251634356, 0.248131676, 0.2412732,
        0.234606887, 0.228128378, 0.221833385, 0.215717692, 0.20977715, 0.204007681, 0.198405272, 0.187685927,
        0.177588357, 0.168082846, 0.159140695, 0.150734206, 0.142836655, 0.135422274, 0.128466227, 0.12194459,
        0.115834329, 0.110113284, 0.104760141, 0.099754417, 0.09507644, 0.090707328, 0.086628967, 0.082823998,
        0.07927579, 0.075968428, 0.072886691, 0.070016034, 0.067342571, 0.064853053, 0.062534858, 0.060375964,
        0.058364941, 0.056490925, 0.054743609, 0.053113222, 0.051590514, 0.050166738, 0.046408775, 0.045302803,
        0.044259051, 0.043271471, 0.042334415, 0.041442618, 0.040591184, 0.039775572, 0.038991583, 0.038235345,
        0.037503301, 0.036792197, 0.036099067, 0.034101935, 0.033456388, 0.032817888, 0.032184949, 0.031556287,
        0.030930816, 0.030307633, 0.029065372, 0.027825562, 0.027205981, 0.026586556, 0.025967391, 0.025348692,
        0.024114005, 0.023498886, 0.021669152, 0.021066668, 0.019292088, 0.018144698, 0.016762709, 0.015451481,
        0.014949794, 0.014224263, 0.013093462, 0.012670686, 0.012070223, 0.011164062, 0.010241734, 0.009731103,
        0.009507687, 0.009212683, 0.008965785, 0.008827751, 0.008710756, 0.008574128, 0.008462605, 0.008446967,
        0.008539475, 0.009015237, 0.009748444, 0.010586023, 0.011359647, 0.011901268, 0.012062153, 0.011735443,
        0.010882215, 0.009561062, 0.007961182, 0.006438984, 0.005558204, 0.006133532, 0.009277754
    ],
    τW::Vector{<:Real}=[
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0.123, 0.117, 0.1, 0.23, 0.174, 0.058, 0, 0.024, 0.027,
        0.036, 0.215, 0.25, 0.136, 0.058, 0.047, 0.036, 0.042, 0.098, 0.044, 0, 0.038,
        0.83, 0, 0, 0.38, 0.289, 0.258, 0.173, 0.008, 0, 0, 0, 0, 0, 0, 0, 0, 0.57, 0.76, 0,
        0.185, 0.291, 0.178, 0.196, 0.112, 0.075, 0.074, 0.07, 0.007, 0, 0, 0, 0.086,
        0.122, 0.132, 0.14, 0.207, 0.259, 0, 0, 0, 0.549, 0.297, 0.462, 0.52, 0.374, 0.222,
        0.614, 0.058, 0.038, 0.03, 0.04, 0.16
    ],
    Sλ::Vector{typeof(0.0u"mW * cm^-2 * nm^-1")}=[
        48.2, 58.4, 51.4, 60.2, 68.6, 75.7, 81.9, 103.7, 105, 107.4, 105.5,
        117.3, 111.7, 109.9, 143.3, 175.8, 182.3, 208, 208.5, 194.6, 183.3, 178.3,
        169.5, 170.5, 164.6, 157.6, 151.7, 146.8, 141.8, 136.9, 131.4, 126, 120, 115,
        110.7, 105, 100, 95, 91, 88, 85, 83, 80, 77, 75, 70, 61, 59, 56, 54, 51, 49, 48, 45,
        42, 41, 40, 39, 38, 34, 33, 32, 31, 30, 29, 28, 26, 25, 24, 24, 23, 22, 21, 19, 16, 15,
        12, 11, 10.7, 10.3, 10, 9.7, 9, 8.8, 8.5, 7.9, 7.4, 6.8, 6.7, 6.6, 6.5, 6.4, 6.2,
        5.9, 5.5, 5.4, 4.8, 4.3, 3.9, 3.5, 3.1, 2.6, 2.3, 1.9, 1.7, 1.5, 1.4, 1.2, 1.1, 1, 1
    ] * u"mW * cm^-2 * nm^-1",

    FD::Matrix{<:Real}=reshape([
            8.00e-5, 6.50e-5, 4.00e-5, 2.30e-5, 1.00e-5, 4.50e-6, 1.00e-6, 1.00e-7, 5.50e-9,
            1.00e-9, 3.50e-10, 1.60e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10,
            1.00e-10, 1.00e-10, 1.00e-3, 9.50e-4, 9.00e-4, 8.00e-4, 7.00e-4, 6.00e-4,
            4.50e-4, 3.00e-4, 1.70e-4, 8.00e-5, 3.30e-5, 1.80e-5, 1.00e-5, 8.00e-6,
            6.30e-6, 5.00e-6, 4.10e-6, 3.50e-6, 3.00e-6, 3.50e-2, 3.30e-2, 3.10e-2,
            2.90e-2, 2.50e-2, 2.20e-2, 1.75e-2, 1.30e-2, 7.50e-3, 4.50e-3, 2.50e-3,
            1.30e-3, 5.60e-4, 2.70e-4, 1.40e-4, 7.10e-5, 4.20e-5, 2.60e-5, 1.70e-5,
            1.55e-1, 1.40e-1, 1.40e-1, 1.30e-1, 1.20e-1, 1.10e-1, 1.00e-1, 9.00e-2,
            8.20e-2, 7.00e-2, 5.20e-2, 3.50e-2, 2.20e-2, 1.00e-2, 4.00e-3, 1.20e-3,
            4.50e-4, 1.90e-4, 8.00e-5, 3.70e-1, 3.70e-1, 3.60e-1, 3.50e-1, 3.30e-1,
            3.10e-1, 2.90e-1, 2.60e-1, 2.30e-1, 2.00e-1, 1.70e-1, 1.50e-1, 1.00e-1,
            7.00e-2, 3.80e-2, 1.60e-2, 5.00e-3, 1.20e-3, 3.00e-4, 5.50e-1, 5.50e-1,
            5.50e-1, 5.30e-1, 5.10e-1, 5.00e-1, 4.70e-1, 4.40e-1, 4.00e-1, 3.60e-1,
            3.10e-1, 2.70e-1, 2.25e-1, 1.80e-1, 1.20e-1, 6.00e-2, 2.50e-2, 4.50e-3,
            7.00e-4, 6.50e-1, 6.50e-1, 6.50e-1, 6.50e-1, 6.20e-1, 6.00e-1, 5.70e-1,
            5.50e-1, 5.00e-1, 4.50e-1, 4.20e-1, 3.80e-1, 3.25e-1, 2.70e-1, 1.90e-1,
            1.15e-1, 5.00e-2, 1.10e-2, 1.20e-3, 7.88e-1, 7.80e-1, 8.00e-1, 8.00e-1,
            8.00e-1, 7.60e-1, 7.35e-1, 7.10e-1, 7.00e-1, 6.50e-1, 6.00e-1, 5.50e-1,
            5.14e-1, 4.50e-1, 3.50e-1, 2.68e-1, 1.50e-1, 6.27e-2, 1.20e-2, 7.48e-1,
            7.40e-1, 7.40e-1, 7.30e-1, 7.20e-1, 7.10e-1, 7.04e-1, 6.90e-1, 6.70e-1,
            6.20e-1, 5.70e-1, 5.30e-1, 5.16e-1, 4.80e-1, 3.90e-1, 2.90e-1, 1.70e-1,
            7.62e-2, 3.00e-2, 7.00e-1, 7.00e-1, 7.00e-1, 6.90e-1, 6.80e-1, 6.80e-1,
            6.60e-1, 6.50e-1, 6.30e-1, 6.00e-1, 5.60e-1, 5.21e-1, 5.00e-1, 4.70e-1,
            3.90e-1, 3.00e-1, 1.85e-1, 9.00e-2, 2.60e-2, 6.51e-1, 6.50e-1, 6.50e-1,
            6.40e-1, 6.30e-1, 6.25e-1, 6.22e-1, 6.00e-1, 5.90e-1, 5.70e-1, 5.50e-1,
            5.20e-1, 4.89e-1, 4.60e-1, 3.90e-1, 3.08e-1, 2.00e-1, 9.55e-2, 2.20e-2
        ], (11, 19)),
    FDQ::Matrix{<:Real}=reshape([
            8.00e-6, 7.00e-6, 5.20e-6, 3.50e-6, 1.70e-6, 5.50e-7, 1.00e-7, 2.50e-8, 6.00e-9,
            1.50e-9, 3.00e-10, 6.00e-11, 1.00e-11, 1.00e-11, 1.00e-11, 1.00e-11, 1.00e-11,
            1.00e-11, 1.00e-11, 6.10e-4, 6.00e-4, 5.50e-4, 4.50e-4, 3.40e-4, 2.30e-4,
            1.20e-4, 5.50e-5, 2.60e-5, 1.20e-5, 6.00e-6, 3.50e-6, 2.00e-6, 1.50e-6,
            1.00e-6, 6.00e-7, 4.00e-7, 2.30e-7, 1.00e-7, 2.40e-2, 2.30e-2, 2.20e-2,
            2.10e-2, 1.80e-2, 1.50e-2, 1.20e-2, 7.50e-3, 4.80e-3, 2.50e-3, 1.20e-3,
            5.00e-4, 2.50e-4, 1.00e-4, 4.50e-5, 2.00e-5, 1.00e-5, 5.50e-6, 2.50e-6,
            1.30e-1, 1.20e-1, 1.10e-1, 1.00e-1, 9.10e-2, 8.50e-2, 7.20e-2, 6.50e-2,
            5.20e-2, 4.00e-2, 2.80e-2, 1.70e-2, 1.00e-2, 3.00e-3, 1.00e-3, 4.00e-4,
            1.70e-4, 6.70e-5, 2.50e-5, 3.40e-1, 3.30e-1, 3.20e-1, 3.00e-1, 2.90e-1,
            2.70e-1, 2.50e-1, 2.00e-1, 1.70e-1, 1.50e-1, 1.20e-1, 8.00e-2, 5.80e-2,
            3.00e-2, 1.30e-2, 6.20e-3, 2.00e-3, 6.00e-4, 1.90e-4, 5.40e-1, 5.30e-1,
            5.10e-1, 5.00e-1, 4.80e-1, 4.50e-1, 4.20e-1, 3.70e-1, 3.20e-1, 2.70e-1,
            2.20e-1, 1.70e-1, 1.20e-1, 8.00e-2, 5.00e-2, 2.30e-2, 8.00e-3, 4.70e-3,
            9.00e-4, 6.50e-1, 6.40e-1, 6.20e-1, 6.10e-1, 6.00e-1, 5.60e-1, 5.10e-1,
            4.70e-1, 4.20e-1, 3.60e-1, 3.00e-1, 2.50e-1, 1.80e-1, 1.30e-1, 8.00e-2,
            5.00e-2, 2.00e-2, 4.60e-3, 1.50e-3, 8.20e-1, 8.10e-1, 8.10e-1, 8.00e-1,
            7.90e-1, 7.50e-1, 7.00e-1, 6.50e-1, 6.00e-1, 5.40e-1, 4.60e-1, 3.80e-1,
            3.20e-1, 2.50e-1, 1.80e-1, 1.20e-1, 6.00e-2, 2.50e-2, 8.00e-3, 8.60e-1,
            8.40e-1, 8.20e-1, 8.00e-1, 7.50e-1, 7.00e-1, 6.80e-1, 6.30e-1, 5.80e-1,
            5.00e-1, 4.40e-1, 3.70e-1, 3.20e-1, 2.40e-1, 1.70e-1, 1.20e-1, 6.00e-2,
            2.70e-2, 1.00e-2, 8.00e-1, 7.90e-1, 7.70e-1, 7.60e-1, 7.30e-1, 6.85e-1,
            6.40e-1, 5.90e-1, 5.40e-1, 4.75e-1, 4.20e-1, 3.65e-1, 3.20e-1, 2.50e-1,
            1.80e-1, 1.30e-1, 7.00e-2, 2.90e-2, 1.10e-2, 7.50e-1, 7.40e-1, 7.30e-1,
            7.20e-1, 7.00e-1, 6.70e-1, 6.10e-1, 5.50e-1, 5.00e-1, 4.50e-1, 4.00e-1,
            3.60e-1, 3.20e-1, 2.60e-1, 1.90e-1, 1.40e-1, 8.00e-2, 3.10e-2, 1.20e-2
        ], (11, 19)),
    S::Vector{<:Real}=[0.2, 0.255, 0.315, 0.365, 0.394, 0.405, 0.405, 0.395, 0.37, 0.343, 0.32]
)
    ndays = length(days)    # number of days
    ntimes = length(hours)  # number of times
    nsteps = ndays * ntimes # total time steps

    # arrays to hold every time step's radiation between 300 and 320 nm in 2 nm steps
    GRλs = fill(0.0, nsteps, nmax)u"mW/nm/cm^2" # wavelength-specific global radiation
    DRRλs = fill(0.0, nsteps, nmax)u"mW/nm/cm^2"# wavelength-specific direct Rayleigh radiation
    DRλs = fill(0.0, nsteps, nmax)u"mW/nm/cm^2" # wavelength-specific direct radiation
    SRλs = fill(0.0, nsteps, nmax)u"mW/nm/cm^2" # wavelength-specific scattered radiation
    GRs = fill(0.0, nsteps)u"mW/cm^2"           # total global radiation
    DRRs = fill(0.0, nsteps)u"mW/cm^2"          # total direct Rayleigh radiation
    DRs = fill(0.0, nsteps)u"mW/cm^2"           # total direct radiation
    SRs = fill(0.0, nsteps)u"mW/cm^2"           # total scattered radiation

    # arrays to hold zenith angles each step
    Zs = fill(0.0, nsteps)u"°"                  # zenith angles
    ZSLs = fill(0.0, nsteps)u"°"                # slope zenith angles
    HHs = fill(0.0, ndays)                      # hour angles
    tsns = fill(0.0, ndays)                     # hour at solar noon
    DOYs = Vector{Int}(undef, nsteps)           # day of year
    times = Vector{Real}(undef, nsteps)         # time

    step = 1
    HH = 0.0 # initialise sunrise hour angle
    tsn = 12.0 # initialise time of solar noon
    for i in 1:ndays
        # arrays to hold radiation for a given hour between 300 and 320 nm in 2 nm steps
        GRINT = fill(0.0, nmax)u"mW/cm^2"   # integrated global radiation component (direct + scattered)
        DRRINT = fill(0.0, nmax)u"mW/cm^2"  # integrated direct Rayleigh radiation component
        DRINT = fill(0.0, nmax)u"mW/cm^2"   # integrated direct radiation component
        SRINT = fill(0.0, nmax)u"mW/cm^2"   # integrated scattered radiation component
        AIλ = fill(0.0, nmax)u"nm"
        GRλ = GRINT * u"1/nm"               # wavelength-specific global radiation component (direct + scattered)
        DRRλ = GRINT * u"1/nm"              # wavelength-specific direct Rayleigh radiation component
        DRλ = GRINT * u"1/nm"               # wavelength-specific direct radiation component
        SRλ = GRINT * u"1/nm"               # wavelength-specific scattered radiation component

        for j in 1:ntimes
            d = days[i]
            t = hours[j]
            h, tsn = hour_angle(t, lonc) # hour angle (radians)
            ζ, δ, z, AR2 = solar_geometry(d=d, lat=lat, h=h, d0=d0, ω=ω, ϵ=ϵ, se=se) # compute ecliptic, declination, zenith angle and (a/r)^2
            Z = uconvert(°, z)
            Zsl = Z
            check_skylight(z, nmax, SRINT, GRINT) # checking zenith angle for possible skylight before sunrise or after sunset
            # testing cos(h) to see if it exceeds +1 or -1
            TDTL = -tan(δ) * tan(lat) # from eq.7 McCullough & Porter 1971
            if abs(TDTL) >= 1 # long day or night
                H = π
            else
                H = abs(acos(TDTL))
            end
            # check if sunrise
            HH = 12 * H / π
            ts = t - tsn
            if !((ts <= 0 && abs(ts) - HH > 0) || (ts > 0 && ts - HH >= 0) || (TDTL >= 1)) # sun is up, proceed
                h, tsn = hour_angle(t, lonc) # hour angle (radians)
                ζ, δ, z, AR2 = solar_geometry(d=d, lat=lat, h=h, d0=d0, ω=ω, ϵ=ϵ, se=se) # compute ecliptic, declination, zenith angle and (a/r)^2 - redundant?
                alt = (π / 2 - z)u"rad"
                cazsun = (sin(δ) - sin(lat) * sin(alt)) / (cos(lat) * cos(alt)) # cos(solar azimuth)
                #      Error checking for instability in trig function
                if cazsun < -0.9999999
                    azsun = π
                else
                    if cazsun > 0.9999999
                        azsun = 0
                    else
                        azsun = acos(cazsun)
                    end
                end
                if h <= 0
                    # Morning
                    dazsun = uconvert(°, azsun).val
                else
                    # Afternoon
                    if sign(lat) < 0
                        dazsun = 360° - uconvert(°, azsun).val
                    else
                        dazsun = 180° + (180° - uconvert(°, azsun).val)
                    end
                end

                cz = cos(z)
                intcz = Int(floor(100.0 * cz + 1.0))
                Z = uconvert(°, z)  # zenith angle in degrees

                # horizon angle - check this works when starting at 0 rather than e.g. 15 deg
                azi = range(0°, stop=360° - 360° / length(hori), length=length(hori))
                ahoriz = hori[argmin(abs.(dazsun .- azi))]

                # slope zenith angle calculation (Eq. 3.15 in Sellers 1965. Physical Climatology. U. Chicago Press)
                if slope > 0°
                    czsl = cos(z) * cos(slope) + sin(z) * sin(slope) * cos(dazsun - aspect)
                    zsl = acos(czsl)
                    Zsl = min(uconvert(°, zsl), 90°) # cap at 90 degrees if sun is below slope horizon
                    intczsl = Int(floor(100.0 * czsl + 1.0))
                else
                    czsl = cz
                    zsl = z
                    Zsl = Z
                    intczsl = intcz
                end

                # refraction correction check
                if z < 1.5358896
                    # skip refraction correction
                else
                    refr = 16.0 + ((z - 1.53589) * 15) / (π / 90)
                    refr = (refr / 60) * (π / 180)
                    z -= refr
                end

                # optical air mass (Rozenberg 1966 formula p.159 in book 'Twilight') ---
                airms = 1.0 / (cos(z) + (0.025 * exp(-11.0 * cos(z))))

                # atmospheric ozone lookup
                # convert latitude in degrees to nearest 10-degree index
                tlat = (lat + 100.0°) / 10.0°
                llat = Int(floor(tlat))
                allat = llat
                ala = allat + 0.5
                if tlat > ala
                    llat += 1
                end
                # clamp llat index to valid range
                mon = month(Date(year, 1, 1) + Day(d - 1)) # month from day of year
                llat = clamp(llat, 1, size(OZ, 1))
                ozone = OZ[llat, mon]  # ozone thickness (cm) from lookup table
                ELEVFCT1, ELEVFCT2, ELEVFCT3, ELEVFCT4 = elev_corr(elev)

                # deal with this:
                # c     mutliplier to correct hourly solar data for horizon angle
                #     if(altdeg.lt.ahoriz)then
                # c	   diffuse only - cut down to diffuse fraction      
                #     TDD(111+IT)=TDD(111+IT)* (0.12 + 0.83 * ((CCMINN(IDAY) + 
                #    &  CCMAXX(IDAY))/ 2. / 100.)) ! from Butt et al. 2010 Agricultural and Forest Meteorology 150 (2010) 361–368
                #     endif
                P = get_pressure(elev) # pressure from elevation

                for N in 1:nmax
                    τλ1 = (P / 101300u"Pa") * τR[N] * ELEVFCT1
                    τλ2 = (25u"km" / amr) * τA[N] * ELEVFCT2
                    τλ3 = (ozone / 0.34) * τO[N] * ELEVFCT3
                    τλ4 = τW[N] * sqrt(airms * cmH2O * ELEVFCT4)
                    τλ = ((float(τλ1) + τλ2 + τλ3) * airms) + τλ4

                    if τλ > 80.0 # making sure that at low sun angles air mass doesn't make τλ too large
                        τλ = 80.0
                    end

                    part1 = Sλ[N] * AR2 * cz
                    part2 = τλ > 0.0 ? exp(-τλ) : 0.0
                    if part2 < 1.0e-24
                        DRλ[N] = 0u"mW / cm^2 / nm"
                    else
                        DRλ[N] = ((ustrip(part1) * part2) / 1000) * u"mW / cm^2 / nm"
                    end

                    # so the integrator doesn't get confused at very low sun angles
                    if DRλ[N] < 1.0e-25u"mW / cm^2 / nm"
                        DRλ[N] = 1.0e-25u"mW / cm^2 / nm"
                    end

                    DRRλ[N] = (Sλ[N] * AR2 * cz) * exp(-float(τλ1) * airms) / 1000.0

                    if alt < ahoriz
                        DRλ[N] = 1.0e-24u"mW / cm^2 / nm"
                        DRRλ[N] = 1.0e-24u"mW / cm^2 / nm"
                    end

                    # Sky (SRλ) and Global Radiation (GRλ)
                    if noscat == false
                        SRλ[N] = 0.0u"mW / cm^2 / nm"
                    elseif iuv
                        if τλ1 >= 0.03
                            GAMR, GAML, SBAR = GAMMA(τλ1)
                            SRλ[N] = (
                                         ((float(GAML[intcz]) + float(GAMR[intcz])) / (2.0 * (1.0 - refl * float(SBAR))))
                                         -
                                         exp(-float(τλ1) * airms)
                                     ) * cz * Sλ[N] * AR2 / 1000.0
                        else
                            SRλ[N] = 0.0u"mW / cm^2 / nm"
                        end
                    else
                        if N > 11
                            SRλ[N] = 0.0u"mW / cm^2 / nm"
                        else
                            I = round(Int, (ustrip(Z) + 5) / 5)
                            FDAV = FD[N, I]
                            FDQDAV = FDQ[N, I]
                            SRλ[N] = (Sλ[N] / π) * (FDAV + FDQDAV * (refl / (1.0 - (refl * S[N])))) / 1000.0
                            SRλ[N] *= AR2
                        end
                    end

                    GRλ[N] = SRλ[N] + DRλ[N]

                    if N == 1
                        SRINT[1] = 0.0u"mW / cm^2"
                        DRRINT[1] = 0.0u"mW / cm^2"
                        DRINT[1] = 0.0u"mW / cm^2"
                        GRINT[1] = 0.0u"mW / cm^2"
                    else
                        AIλ[N] = Iλ[N]
                        AIλ[N-1] = Iλ[N-1]

                        Δλ = AIλ[N] - AIλ[N-1]

                        DRINT[N] = DRINT[N-1] + (Δλ * DRλ[N-1]) + (0.5 * Δλ * (DRλ[N] - DRλ[N-1]))
                        DRRINT[N] = DRRINT[N-1] + (Δλ * DRRλ[N-1]) + (0.5 * Δλ * (DRRλ[N] - DRRλ[N-1]))
                        SRINT[N] = SRINT[N-1] + (Δλ * SRλ[N-1]) + (0.5 * Δλ * (SRλ[N] - SRλ[N-1]))
                        GRINT[N] = GRINT[N-1] + (Δλ * GRλ[N-1]) + (0.5 * Δλ * (GRλ[N] - GRλ[N-1]))
                    end
                end
            else # sunrise, sunset or long day

            end
            # Store into row `step`
            GRλs[step, :] .= GRλ
            DRRλs[step, :] .= DRRλ
            DRλs[step, :] .= DRλ
            SRλs[step, :] .= SRλ
            GRs[step] = GRINT[nmax]
            DRRs[step] = DRRINT[nmax]
            DRs[step] = DRINT[nmax]
            SRs[step] = SRINT[nmax]
            Zs[step] = Z
            ZSLs[step] = Zsl
            DOYs[step] = d
            times[step] = t
            Zs[Zs.>90°] .= 90°
            step += 1
        end
        HHs[i] = HH     # save today's sunrise hour angle
        tsns[i] = tsn   # save today's time of sunrise
    end
    return (
        Zenith=Zs,
        ZenithSlope=ZSLs,
        HHsr=HHs,
        tsn=tsns,
        doy=DOYs,
        hour=times,
        Rayleigh=DRRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        Direct=DRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        Scattered=SRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        Global=GRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        λ=Iλ,
        λRayleigh=DRRλs .* (10u"W/m^2" / 1u"mW/cm^2"),
        λDirect=DRλs .* (10u"W/m^2" / 1u"mW/cm^2"),
        λScattered=SRλs .* (10u"W/m^2" / 1u"mW/cm^2"),
        λGlobal=GRλs .* (10u"W/m^2" / 1u"mW/cm^2"),
    )
end

function get_longwave(;
    elev,
    rh,
    tair,
    tsurf,
    slep,
    sle,
    cloud,
    viewf
)
    # Longwave radiation (handle both IR modes)
    # Constants
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ) # Stefan-Boltzmann constant, W/m^2/K^4
    P_atmos = get_pressure(elev)
    wet_air_out = wet_air(u"K"(tair); rh=rh, P_atmos=P_atmos)

    # Atmospheric radiation
    P_vap = wet_air_out.P_vap
    arad = u"W/m^2"(ustrip(1.72 * (ustrip(u"kPa"(P_vap))/ustrip(u"K"(tair))) ^ (1.0/7.0)) * σ * (u"K"(tair)) ^ 4) # Campbell and Norman 1998 eq. 10.10 to get emissivity of sky
    #arad = ((0.0000092 * (u"K"(tair))^2) * σ * (u"K"(tair))^4) / 1u"K^2" # Swinbank, Eq. 10.11 in Campbell and Norman 1998

    # Cloud radiation temperature (shade approximation, TAIR - 2°C)
    crad = σ * slep * (u"K"(tair) - 2u"K")^4

    # Hillshade radiation temperature (approximated as air temperature)
    hrad = σ * slep * (u"K"(tair))^4

    # Ground surface radiation temperature
    srad = σ * sle * (u"K"(tsurf))^4

    # Clear sky fraction
    clr = 1.0 - cloud / 100.0
    clear = arad * clr
    clod = crad * (cloud / 100)
    qradsk = (clear + clod) * ((100 - shade) / 100.0)
    qradvg = (shade / 100.0) * hrad
    qradgr = ((100.0 - shade) / 100.0) * srad + (shade / 100.0) * hrad
    qradhl = hrad
    qrad = (qradsk + qradvg) * viewf + qradhl * (1.0 - viewf) - qradgr
    tsky = (((qradsk + qradvg) * viewf + qradhl * (1.0 - viewf)) / σ)^0.25
    return tsky, qrad, arad, crad, hrad, srad, qradsk, qradvg, qradgr, qradhl
end