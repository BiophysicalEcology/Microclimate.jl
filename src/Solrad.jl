"""
    get_pressure(h::Quantity;
                 h_ref::Quantity = 0u"m",
                 P_ref::Quantity = 101325u"Pa",
                 L_ref::Quantity = -0.0065u"K/m",
                 T_ref::Quantity = 288u"K",
                 g_0::Quantity = 9.80665u"m/s^2",
                 M::Quantity = 0.0289644u"kg/mol") -> Quantity

Computes atmospheric pressure at a given altitude using the barometric formula,
assuming a constant temperature lapse rate (standard tropospheric approximation).

# Arguments
- `h`: Altitude at which to compute pressure (with length units, e.g. `u"m"`).
- `h_ref`: Reference altitude (default: `0u"m"`).
- `P_ref`: Pressure at the reference altitude (default: `101325u"Pa"`, standard sea-level pressure).
- `L_ref`: Temperature lapse rate (default: `-0.0065u"K/m"`).
- `T_ref`: Temperature at the reference altitude (default: `288u"K"`).
- `g_0`: Standard gravitational acceleration (default: `9.80665u"m/s^2"`).
- `M`: Molar mass of dry air (default: `0.0289644u"kg/mol"`).

# Returns
- `P_a`: Atmospheric pressure at altitude `h` (with pressure units, e.g. `u"Pa"`).

# Notes
- Uses the simplified barometric formula assuming a linear lapse rate and ideal gas behavior.
- The universal gas constant `R` is used from the `Unitful` package.

# Example
```julia
using Unitful

P = get_pressure(1500u"m")
"""
function get_pressure(h::Quantity;
    h_ref::Quantity = 0u"m",
    P_ref::Quantity = 101325u"Pa",
    L_ref::Quantity = -0.0065u"K/m",
    T_ref::Quantity = 288u"K",
    g_0::Quantity = 9.80665u"m/s^2",
    M::Quantity = 0.0289644u"kg/mol")

P_a = P_ref * (1 + (L_ref / T_ref) * h) ^ ((-g_0 * M) / (R * L_ref))

return P_a
end

"""
    hour_angle(t::Quantity, lonc::Quantity) -> Quantity

Compute the solar hour angle `h` in radians.

# Arguments
- `t`: Local solar time (e.g., `14.0u"hr"`)
- `lonc`: Longitude correction (e.g., `0.5u"hr"`)

# Returns
- Hour angle `h` as a `Quantity` in radians
- Time at solar noon, `tsn`  as a time in hours

# Reference
McCullough & Porter 1971, Eq. 6
"""
function hour_angle(t::Real, lonc::Real = 0)
    tsn = 12.0 + lonc                        # solar noon time
    h = (2π / 24) * (t - tsn) * u"rad"      # convert hours to radians
    return h, tsn
end

"""
    solar_geometry(d::Real, lat::Quantity, h::Quantity; d0::Real = 80, ω::Real = 2π/365, ϵ::Real = 0.0167, SE::Real = 0.39779)

Computes key solar geometry parameters based on McCullough & Porter (1971):

- `ζ`: Auxiliary solar longitude (radians)
- `δ`: Solar declination (radians)
- `Z`: Solar zenith angle (radians)
- `AR2`: Square of Earth-to-Sun radius factor (unitless)

# Arguments
- `d`: Day of year (1–365)
- `lat`: Latitude (with angle units, e.g. `u"°"` or `u"rad"`)
- `h`: Hour angle (radians)
- `d0`: Reference day (default: 80)
- `ω`: Angular frequency of Earth’s orbit (default: `2π/365`)
- `ϵ`: Orbital eccentricity (default: `0.0167`)
- `SE`: Constant for solar declination amplitude (default: `0.39779`)

# Returns
Tuple: `(ζ, δ, Z, AR2)` with angle quantities in radians and AR2 unitless.

# Reference
McCullough & Porter (1971)
"""
function solar_geometry(;
    d::Real = 1.,
    lat::Quantity = 83.07305u"°",
    h::Quantity = -2.87979u"rad",
    d0::Real = 80,
    ω::Real = 2π / 365,
    ϵ::Real = 0.0167238,
    SE::Real = 0.39779
)
    ζ = (ω * (d - d0)) + 2ϵ * (sin(ω * d) - sin(ω * d0))          # Eq.5
    δ = (asin(SE * sin(ζ)) * sign(lat))                             # Eq.4
    cosZ = cos(lat) * cos(δ) * cos(h) + sin(lat) * sin(δ)           # Eq.3
    Z = acos(cosZ)u"rad"                                            # Zenith angle
    AR2 = 1 + (2ϵ) * cos(ω * d)                                     # Eq.2
    δ = δ*u"rad"
    ζ = ζ*u"rad"
    return ζ, δ, Z, AR2
end

"""
    check_skylight(Z, NMAX, SRINT, GRINT)

Checks for possible skylight before sunrise or after sunset based on zenith angle.
Modifies SRINT and GRINT at index `NMAX` if skylight is present.

# Arguments
- `Z::Quantity`: Zenith angle
- `NMAX::Int`: Index into result arrays
- `SRINT::Vector{Quantity}`: Scattered radiation array [W/m²]
- `GRINT::Vector{Quantity}`: Global radiation array [W/m²]
"""
function check_skylight(
    Z::Quantity,
    NMAX::Int,
    SRINT::Vector,
    GRINT::Vector)
    Zdeg = uconvert(°, Z).val # convert to degrees
    if Zdeg < 107.0
        if Zdeg > 88.0
            Elog = 41.34615384 - 0.423076923 * Zdeg
            Skylum = (10.0 ^ Elog)*1.46E-03u"mW * cm^-2"
            SRINT[NMAX] = Skylum
            GRINT[NMAX] = SRINT[NMAX]
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

    ELEVFCT2 = 8.35656e-7  * elev_km^6 -
              6.26384e-5  * elev_km^5 +
              1.86967e-3  * elev_km^4 -
              2.82585e-2  * elev_km^3 +
              2.26739e-1  * elev_km^2 -
              9.25268e-1  * elev_km +
              1.71321

    ELEVFCT3 = 1.07573e-6  * elev_km^5 -
              5.14511e-5  * elev_km^4 +
              7.97960e-4  * elev_km^3 -
              4.90904e-3  * elev_km^2 +
              2.99258e-3  * elev_km +
              1.00238

    ELEVFCT4 = 1.0

    return (ELEVFCT1, ELEVFCT2, ELEVFCT3, ELEVFCT4)
end

function GAMMA(TAU1::Float64)

    CHX = zeros(Float64, 101)
    CHY = zeros(Float64, 101)
    CFA = zeros(Float64, 3)
    AMU = zeros(Float64, 101)
    X1  = zeros(Float64, 101)
    Y1  = zeros(Float64, 101)
    X2  = zeros(Float64, 101)
    Y2  = zeros(Float64, 101)
    AIL = zeros(Float64, 101)
    AI  = zeros(Float64, 30)
    XA  = zeros(Float64, 4)
    XB  = zeros(Float64, 8)
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
        AIL[I + 1] = CNU2
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

    AI[1]  = XB[1] + XB[5] - 8.0 / 3.0
    AI[2]  = XB[2] + XB[6]
    AI[3]  = XB[3] + XB[7]
    AI[4]  = XB[1] - XB[5] - 8.0 / 3.0
    AI[5]  = XB[2] - XB[6]
    AI[6]  = XB[3] - XB[7]
    AI[7]  = XB[4] - XB[8]
    AI[8]  = XA[1] + XA[3]
    AI[9]  = XA[2] + XA[4]
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
    SBAR[] = 1.0 - 0.375 * AI[12] * (AI[4] - AI[6]) *
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

function dchxy(t::Vector{Float64}, x::Matrix{Float64}, y::Matrix{Float64}, ng::Int)
    # Constants
    PI = π

    # Temporary variables
    a = zeros(Float64, ng)
    b = zeros(Float64, ng)
    c = zeros(Float64, ng)
    d = zeros(Float64, ng)
    h = zeros(Float64, ng)
    k = zeros(Float64, ng)
    p = zeros(Float64, ng)
    q = zeros(Float64, ng)

    for j in 1:ng
        tj = t[j]
        t2 = tj^2
        t3 = tj^3
        t4 = tj^4

        # Compute coefficients for each angular quadrature point
        a[j] = 1.0 + 0.5 * t2
        b[j] = 0.5 * (3.0 + t2)
        c[j] = 0.5 * (5.0 + t2)
        d[j] = (35.0 + 30.0 * t2 + 3.0 * t4) / 8.0
    end

    # Main integration loop over each pair of directions
    for i = 1:ng
        ti = t[i]
        t2i = ti^2
        t3i = ti^3
        t4i = ti^4
        t5i = ti^5
        t6i = ti^6
        t7i = ti^7

        for j = 1:ng
            tj = t[j]
            t2j = tj^2
            t3j = tj^3
            t4j = tj^4
            t5j = tj^5
            t6j = tj^6
            t7j = tj^7

            den = (ti + tj)^3 * (1 + ti * tj)

            # Compute X[i,j]
            xij = ti + tj
            xij += (3.0 * (ti + tj)) / (1.0 + ti * tj)
            xij += 1.5 * (ti * t2j + tj * t2i)
            xij += 0.75 * (ti * t4j + tj * t4i)
            xij += (5.0 / 8.0) * (ti * t6j + tj * t6i)
            xij /= den
            x[i, j] = xij

            # Compute Y[i,j]
            yij = ti + tj
            yij += (3.0 * (ti + tj)) / (1.0 + ti * tj)
            yij += 1.5 * (ti * t2j + tj * t2i)
            yij += 0.75 * (ti * t4j + tj * t4i)
            yij += (5.0 / 8.0) * (ti * t6j + tj * t6i)
            yij += (63.0 / 128.0) * (ti * t8(tj) + tj * t8(ti)) # helper function below
            yij /= den
            y[i, j] = yij
        end
    end

    return x, y
end

"""
    solrad(; days, hours, lat...[, year, lonc, elev, slope, aspect, hori, REFL, cmH2O, ϵ,
           ω, SE, d0, lamb, IUV, NOSCAT, AMR, NMAX, Iλ, OZ, τR, τO, τA, τW, Sλ, FD, FDQ,
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
- `REFL::Real=0.10`: Ground reflectance (albedo), fraction [0, 1].
- `cmH2O::Real=1`: Precipitable water in cm for atmospheric column (e.g. 0.1: dry, 1.0: moist, 2.0: humid).
- `ϵ::Real=0.0167238`: Orbital eccentricity of Earth.
- `ω::Real=2π/365`: Mean angular orbital velocity of Earth (radians/day).
- `SE::Real=0.39779`: Precomputed solar elevation constant.
- `d0::Real=80`: Reference day for declination calculations.
- `lamb::Bool=false`: If `true`, returns wavelength-specific irradiance components.
- `IUV::Bool=false`: If `true`, uses the full gamma-function model for diffuse radiation (expensive).
- `NOSCAT::Bool=true`: If `true`, disables scattered light computations (faster).
- `AMR::Quantity=25.0u"km"`: Mixing ratio height of the atmosphere.
- `NMAX::Integer=111`: Maximum number of wavelength intervals.
- `Iλ::Vector{Quantity}`: Vector of wavelength bins (e.g. in `nm`).
- `OZ::Matrix{Float64}`: Ozone column depth table indexed by latitude band and month (size 19×12).
- `τR`, `τO`, `τA`, `τW`: Vectors of optical depths per wavelength for Rayleigh scattering, ozone, aerosols, and water vapor.
- `Sλ::Vector{Quantity}`: Solar spectral irradiance per wavelength bin (e.g. in `mW * cm^-2 * nm^-1`).
- `FD`, `FDQ`: Auxiliary data vectors for biologically effective radiation models.
- `S::Vector{Float64}`: Effective sensitivity weights for biologically effective radiation.
- `ER::Vector{Float64}`: Spectral effectiveness response data.
- `ERλ::Vector{Quantity}`: Spectral effectiveness wavelengths (e.g. `mW/cm²/nm/s`).

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
"""
function solrad(;
    days::Vector{Float64}=collect(1.0:365.0),
    hours::Vector{Float64}=collect(0.0:24.0),
    year::Real=2001,
    lat::Quantity=43.07305u"°",
    lonc::Real=0.0, # longitude correction, hours
    elev::Quantity=0.0u"m", # elevation, m
    slope::Quantity=0u"°",
    aspect::Quantity=0u"°",
    hori=fill(0.0, 24) * u"°",
    REFL::Real=0.10, # substrate solar reflectivity (decimal %)
    cmH2O::Real=1, # precipitable cm H2O in air column, 0.1 = VERY DRY; 1 = MOIST AIR CONDITIONS; 2 = HUMID, TROPICAL CONDITIONS (note this is for the whole atmospheric profile, not just near the ground)
    ϵ::Real=0.0167238,
    ω::Real=2π / 365,
    SE::Real=0.39779,
    d0::Real=80,
    lamb::Bool=false, # Return wavelength-specific solar radiation output?
    IUV::Bool=false, # Use gamma function for scattered solar radiation? (computationally intensive)
    NOSCAT::Bool=true,
    AMR::Quantity=25.0u"km",
    NMAX::Integer=111, # maximum number of wavelengths
    Iλ=[ # wavelengths across which to integrate
        290, 295, 300, 305, 310, 315, 320, 330, 340, 350, 360, 370, 380, 390,
        400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700,
        720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020,
        1080, 1100, 1120, 1140, 1160, 1180, 1200, 1220, 1240, 1260, 1280, 1300, 1320,
        1380, 1400, 1420, 1440, 1460, 1480, 1500, 1540, 1580, 1600, 1620, 1640, 1660,
        1700, 1720, 1780, 1800, 1860, 1900, 1950, 2000, 2020, 2050, 2100, 2120, 2150,
        2200, 2260, 2300, 2320, 2350, 2380, 2400, 2420, 2450, 2490, 2500, 2600, 2700,
        2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000
    ] * u"nm",
    OZ=reshape([
            0.31, 0.30, 0.30, 0.27, 0.34, 0.38, 0.43, 0.45, 0.41, 0.37, 0.34, 0.31,
            0.31, 0.31, 0.31, 0.28, 0.35, 0.40, 0.44, 0.46, 0.42, 0.38, 0.36, 0.32,
            0.32, 0.31, 0.31, 0.29, 0.34, 0.39, 0.43, 0.45, 0.43, 0.40, 0.38, 0.34,
            0.32, 0.31, 0.30, 0.30, 0.33, 0.38, 0.41, 0.42, 0.42, 0.40, 0.39, 0.35,
            0.31, 0.30, 0.29, 0.30, 0.32, 0.36, 0.39, 0.40, 0.40, 0.39, 0.37, 0.35,
            0.30, 0.29, 0.28, 0.29, 0.31, 0.33, 0.35, 0.37, 0.38, 0.37, 0.34, 0.32,
            0.27, 0.28, 0.26, 0.27, 0.28, 0.28, 0.29, 0.31, 0.32, 0.32, 0.29, 0.29,
            0.24, 0.25, 0.24, 0.25, 0.25, 0.25, 0.25, 0.26, 0.26, 0.26, 0.26, 0.25,
            0.23, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.23,
            0.22, 0.22, 0.23, 0.23, 0.24, 0.24, 0.24, 0.23, 0.23, 0.22, 0.22, 0.22,
            0.23, 0.24, 0.24, 0.25, 0.26, 0.25, 0.25, 0.24, 0.24, 0.23, 0.23, 0.23,
            0.24, 0.26, 0.26, 0.27, 0.28, 0.27, 0.26, 0.26, 0.26, 0.25, 0.25, 0.25,
            0.27, 0.28, 0.29, 0.30, 0.30, 0.30, 0.29, 0.28, 0.27, 0.26, 0.26, 0.27,
            0.30, 0.32, 0.33, 0.34, 0.34, 0.33, 0.31, 0.30, 0.28, 0.27, 0.28, 0.29,
            0.32, 0.36, 0.38, 0.38, 0.37, 0.35, 0.33, 0.31, 0.30, 0.28, 0.29, 0.30,
            0.33, 0.39, 0.42, 0.40, 0.39, 0.36, 0.34, 0.32, 0.30, 0.28, 0.30, 0.31,
            0.34, 0.40, 0.45, 0.42, 0.40, 0.36, 0.34, 0.31, 0.29, 0.28, 0.29, 0.31,
            0.34, 0.40, 0.46, 0.43, 0.40, 0.36, 0.33, 0.30, 0.28, 0.27, 0.29, 0.31,
            0.33, 0.39, 0.46, 0.42, 0.39, 0.34, 0.32, 0.30, 0.27, 0.26, 0.28, 0.30
        ], (19, 12)),
    τR=[
        1.41, 1.31, 1.23, 1.14, 1.05, 0.99, 0.92, 0.81, 0.72, 0.63, 0.56, 0.5,
        0.45, 0.4, 0.36, 0.3, 0.25, 0.2, 0.17, 0.15, 0.12, 0.1, 0.09, 0.08, 0.07, 0.06,
        0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01,
        0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ],
    τO=[
        11.5, 6.3, 3.2, 1.62, 0.83, 0.44, 0.26, 0.03, 0.02, 0.01, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ],
    τA=[0.269904738, 0.266147825, 0.262442906, 0.258789404, 0.255186744, 0.251634356, 0.248131676, 0.2412732,
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
    τW=[
        1.06, 1.03, 1.01, 0.99, 0.97, 0.95, 0.93, 0.91, 0.89, 0.87, 0.85, 0.84, 0.83,
        0.82, 0.81, 0.8, 0.8, 0.79, 0.78, 0.77, 0.76, 0.75, 0.74, 0.73, 0.72, 0.71, 0.7,
        0.7, 0.69, 0.68, 0.67, 0.66, 0.66, 0.65, 0.64, 0.63, 0.63, 0.62, 0.61, 0.6, 0.59,
        0.58, 0.58, 0.57, 0.56, 0.55, 0.54, 0.53, 0.52, 0.51, 0.5, 0.49, 0.48, 0.47,
        0.46, 0.45, 0.44, 0.43, 0.42, 0.41, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34,
        0.33, 0.32, 0.31, 0.3, 0.29, 0.28, 0.27, 0.26, 0.25, 0.24, 0.23, 0.22, 0.21,
        0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13, 0.12, 0.11, 0.1, 0.09, 0.08, 0.07,
        0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ],
    Sλ=[
           48.2, 58.4, 51.4, 60.2, 68.6, 75.7, 81.9, 103.7, 105, 107.4, 105.5,
           117.3, 111.7, 109.9, 143.3, 175.8, 182.3, 208, 208.5, 194.6, 183.3, 178.3,
           169.5, 170.5, 164.6, 157.6, 151.7, 146.8, 141.8, 136.9, 131.4, 126, 120, 115,
           110.7, 105, 100, 95, 91, 88, 85, 83, 80, 77, 75, 70, 61, 59, 56, 54, 51, 49, 48, 45,
           42, 41, 40, 39, 38, 34, 33, 32, 31, 30, 29, 28, 26, 25, 24, 24, 23, 22, 21, 19, 16, 15,
           12, 11, 10.7, 10.3, 10, 9.7, 9, 8.8, 8.5, 7.9, 7.4, 6.8, 6.7, 6.6, 6.5, 6.4, 6.2,
           5.9, 5.5, 5.4, 4.8, 4.3, 3.9, 3.5, 3.1, 2.6, 2.3, 1.9, 1.7, 1.5, 1.4, 1.2, 1.1, 1, 1
       ] * 10 * u"mW * cm^-2 * nm^-1",
    FD=[
        8.00e-05, 6.50e-05, 4.00e-05, 2.30e-05, 1.00e-05, 4.50e-06, 1.00e-06, 1.00e-07, 5.50e-09,
        1.00e-09, 3.50e-10, 1.60e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10,
        1.00e-10, 1.00e-03, 9.50e-04, 9.00e-04, 8.00e-04, 7.00e-04, 6.00e-04, 4.50e-04, 3.00e-04,
        1.70e-04, 8.00e-05, 3.30e-05, 1.80e-05, 1.00e-05, 8.00e-06, 6.30e-06, 5.00e-06, 4.10e-06,
        3.50e-06, 3.00e-06, 3.50e-02, 3.30e-02, 3.10e-02, 2.90e-02, 2.50e-02, 2.20e-02, 1.75e-02,
        1.30e-02, 7.50e-03, 4.50e-03, 2.50e-03, 1.30e-03, 5.60e-04, 2.70e-04, 1.40e-04, 7.10e-05,
        4.20e-05, 2.60e-05, 1.70e-05, 1.55e-01, 1.40e-01, 1.40e-01, 1.30e-01, 1.20e-01, 1.10e-01,
        1.00e-01, 9.00e-02, 8.20e-02, 7.00e-02, 5.20e-02, 3.50e-02, 2.20e-02, 1.00e-02, 4.00e-03,
        1.20e-03, 4.50e-04, 1.90e-04, 8.00e-05, 3.70e-01, 3.70e-01, 3.60e-01, 3.50e-01, 3.30e-01,
        3.10e-01, 2.90e-01, 2.60e-01, 2.30e-01, 2.00e-01, 1.70e-01, 1.50e-01, 1.00e-01, 7.00e-02,
        3.80e-02, 1.60e-02, 5.00e-03, 1.20e-03, 3.00e-04, 5.50e-01, 5.50e-01, 5.50e-01, 5.30e-01,
        5.10e-01, 5.00e-01, 4.70e-01, 4.40e-01, 4.00e-01, 3.60e-01, 3.10e-01, 2.70e-01, 2.25e-01,
        1.80e-01, 1.20e-01, 6.00e-02, 2.50e-02, 4.50e-03, 7.00e-04, 6.50e-01, 6.50e-01, 6.50e-01,
        6.50e-01, 6.20e-01, 6.00e-01, 5.70e-01, 5.50e-01, 5.00e-01, 4.50e-01, 4.20e-01, 3.80e-01,
        3.25e-01, 2.70e-01, 1.90e-01, 1.15e-01, 5.00e-02, 1.10e-02, 1.20e-03, 7.88e-01, 7.80e-01,
        8.00e-01, 8.00e-01, 8.00e-01, 7.60e-01, 7.35e-01, 7.10e-01, 7.00e-01, 6.50e-01, 6.00e-01,
        5.50e-01, 5.14e-01, 4.50e-01, 3.50e-01, 2.68e-01, 1.50e-01, 6.27e-02, 1.20e-02, 7.48e-01,
        7.40e-01, 7.40e-01, 7.30e-01, 7.20e-01, 7.10e-01, 7.04e-01, 6.90e-01, 6.70e-01, 6.20e-01,
        5.70e-01, 5.30e-01, 5.16e-01, 4.80e-01, 3.90e-01, 2.90e-01, 1.70e-01, 7.62e-02, 3.00e-02,
        7.00e-01, 7.00e-01, 7.00e-01, 6.90e-01, 6.80e-01, 6.80e-01, 6.60e-01, 6.50e-01, 6.30e-01,
        6.00e-01, 5.60e-01, 5.21e-01, 5.00e-01, 4.70e-01, 3.90e-01, 3.00e-01, 1.85e-01, 9.00e-02,
        2.60e-02, 6.51e-01, 6.50e-01, 6.50e-01, 6.40e-01, 6.30e-01, 6.25e-01, 6.22e-01, 6.00e-01,
        5.90e-01, 5.70e-01, 5.50e-01, 5.20e-01, 4.89e-01, 4.60e-01, 3.90e-01, 3.08e-01, 2.00e-01,
        9.55e-02, 2.20e-02
    ],
    FDQ=[
        8.00e-06, 7.00e-06, 5.20e-06, 3.50e-06, 1.70e-06, 5.50e-07, 1.0e-07, 2.50e-08,
        6.00e-09, 1.50e-09, 3.00e-10, 6.00e-11, 1.00e-11, 1.00e-11, 1.00e-11, 1.00e-11,
        1.00e-11, 1.00e-11, 1.00e-11, 6.10e-04, 6.00e-04, 5.50e-04, 4.50e-04, 3.40e-04,
        2.30e-04, 1.20e-04, 5.50e-05, 2.60e-05, 1.20e-05, 6.00e-06, 3.50e-06, 2.00e-06,
        1.50e-06, 1.00e-06, 6.00e-07, 4.00e-07, 2.30e-07, 1.00e-07, 2.40e-02, 2.30e-02,
        2.20e-02, 2.10e-02, 1.80e-02, 1.50e-02, 1.20e-02, 7.50e-03, 4.80e-03, 2.50e-03,
        1.20e-03, 5.00e-04, 2.50e-04, 1.00e-04, 4.50e-05, 2.00e-05, 1.00e-05, 5.50e-06,
        2.50e-06, 1.30e-01, 1.20e-01, 1.10e-01, 1.00e-01, 9.10e-02, 8.50e-02, 7.20e-02,
        6.50e-02, 5.20e-02, 4.00e-02, 2.80e-02, 1.70e-02, 1.00e-02, 3.00e-03, 1.00e-03,
        4.00e-04, 1.70e-04, 6.70e-05, 2.50e-05, 3.40e-01, 3.30e-01, 3.20e-01, 3.00e-01,
        2.90e-01, 2.70e-01, 2.50e-01, 2.00e-01, 1.70e-01, 1.50e-01, 1.20e-01, 8.00e-02,
        5.80e-02, 3.00e-02, 1.30e-02, 6.20e-03, 2.00e-03, 6.00e-04, 1.90e-04, 5.40e-01,
        5.30e-01, 5.10e-01, 5.00e-01, 4.80e-01, 4.50e-01, 4.20e-01, 3.70e-01, 3.20e-01,
        2.70e-01, 2.20e-01, 1.70e-01, 1.20e-01, 8.00e-02, 5.00e-02, 2.30e-02, 8.00e-03,
        4.70e-03, 9.00e-04, 6.50e-01, 6.40e-01, 6.20e-01, 6.10e-01, 6.00e-01, 5.60e-01,
        5.10e-01, 4.70e-01, 4.20e-01, 3.60e-01, 3.00e-01, 2.50e-01, 1.80e-01, 1.30e-01,
        8.00e-02, 5.00e-02, 2.00e-02, 4.60e-03, 1.50e-03, 8.20e-01, 8.10e-01, 8.10e-01,
        8.00e-01, 7.90e-01, 7.50e-01, 7.00e-01, 6.50e-01, 6.00e-01, 5.40e-01, 4.60e-01,
        3.80e-01, 3.20e-01, 2.50e-01, 1.80e-01, 1.30e-01, 7.00e-02, 2.90e-02, 1.10e-02,
        7.50e-01, 7.40e-01, 7.30e-01, 7.20e-01, 7.00e-01, 6.70e-01, 6.10e-01, 5.50e-01,
        5.00e-01, 4.50e-01, 4.00e-01, 3.60e-01, 3.20e-01, 2.60e-01, 1.90e-01, 1.40e-01,
        8.00e-02, 3.10e-02, 1.20e-02
    ],
    S=[0.2, 0.255, 0.315, 0.365, 0.394, 0.405, 0.405, 0.395, 0.37, 0.343, 0.32],
    ER=[6.19, 6.86, 11.6, 25.1, 224.0, 560.0, 1160.0],
    ERλ=[12, 15, 20, 34, 88, 224, 300, 430, 600, 830, 1160] * u"mW/cm^2/nm/s"
)

    # ILAM(N)      NTH WAVELENGTH (NM)
    # DRLAM(N)     DIRECT COMPONENT RADIATION-NTH SPECTRAL ESTIMATE
    # DRRLAM(N)    DIRECT RAYLEIGH COMPONENT RADIATION-NTH SPECTRAL ESTIMATE
    # SRLAM(N)     DIFFUSE COMPONENT RADIATION-NTH SPECTRAL ESTIMATE
    # GRLAM(N)     GLOBAL RADIATION-NTH SPECTRAL ESTIMATE (GRLAM=DRLAM+SRLAM)
    # EPLAM(N)     ERYTHYMAL PRODUCT (GRLAM TIMES ACTION SPECTRUM ER)-NTH SPEC-
    #              TRAL ESTIMATE
    # IELAM(N)     WAVELENGTHS-BETWEEN 300 AND 320 NMS IN 2 NM STEPS
    # GPL(N)       GLOBAL RADIATION BETWEEN 300 AND 320 NM IN 2 NM STEPS
    # EPL(N)       ERYTHEMAL PRODUCT BETWEEN 300 AND 320 NM IN 2 NM STEPS
    # LAMAX        WVLTH AT WHICH MAX OF ERYTHEMA PROD. CURVE OCCURS

    # THE FOLLOWING CONTAIN THE INTEGRATED VALUES FOR THE SPECTRAL ESTIMATE
    # DESIGNATED, THAT IS, INTEGRATED OVER WAVELENGTH FROM 290 NM TO THE NTH
    # WAVELENGTH ILAM(N)-------

    # DRINT(N)     DIRECT RADIATION COMPONENT
    # DRRINT(N)    DIRECT RAYLEIGH RADIATION COMPONENT
    # SRINT(N)     SCATTERED RADIATION COMPONENT
    # GRINT(N)     GLOBAL RADIATION COMPONENT
    # EPINT(N)     ERYTHEMA PRODUCT
    # REPINT(N)    RECIPROCAL OF EPINT

    # Define ER as a constant array - action spectrum 


    # Define S as a constant array


    # CRIPPS AND RAMSAY ( BR. JOURNAL OF DERMATOL.-JUNE 1970) HAVE
    # DETERMINED THE MINIMAL ENERGY TO PRODUCE AN ERYTHEMA (SUNBURN) AS
    # A FUNCTION OF WAVELENGTH. THE ARRAY ER(N) ARE THESE VALUES FOR WAVE-
    # LENGTHS FROM 290 TO 320 NANOMETERS (NM) IN  5 NM STEPS. THE ARRAY
    # ERλ(N) ARE VALUES FOR WAVELENGTHS FROM 300 TO 320 NM IN 2 NM STEPS.

    # Arrays with NMAX elements set to 0.0
    GRINT = fill(0.0, NMAX)u"mW/cm^2" # GLOBAL RADIATION COMPONENT (DIRECT + SCATTERED)
    DRRINT = fill(0.0, NMAX)u"mW/cm^2" # DIRECT RAYLEIGH RADIATION COMPONENT
    DRINT = fill(0.0, NMAX)u"mW/cm^2" # DIRECT RADIATION COMPONENT
    SRINT = fill(0.0, NMAX)u"mW/cm^2" # SCATTERED RADIATION COMPONENT
    AIλ = fill(0.0, NMAX)u"nm"
    GRλ = fill(0.0, NMAX)u"mW/cm^2/nm" # integrated global radiation component

    EPINT = fill(0.0, 7) * u"mW * cm^-2"
    REPINT = fill(0.0, 7) * u"cm^2 * mW^-1"
    EPλ = fill(0.0, 7) * u"mW/cm^2/nm"

    DRRλ = GRINT * u"1/nm"
    DRλ = GRINT * u"1/nm"
    SRλ = GRINT * u"1/nm"
    P = get_pressure(elev)

    ndays = length(days)
    ntimes = length(hours)
    nsteps = ndays * ntimes  # total time steps
    # Radiation arrays indexed by [day, hour]
    GRINTs = fill(0.0, nsteps, NMAX)u"mW/cm^2"   # wavelength-specific GLOBAL RADIATION
    DRRINTs = fill(0.0, nsteps, NMAX)u"mW/cm^2"  # wavelength-specific DIRECT RAYLEIGH RADIATION
    DRINTs = fill(0.0, nsteps, NMAX)u"mW/cm^2"   # wavelength-specific DIRECT RADIATION
    SRINTs = fill(0.0, nsteps, NMAX)u"mW/cm^2"   # wavelength-specific SCATTERED RADIATION
    GRs = fill(0.0, nsteps)u"mW/cm^2"   # total GLOBAL RADIATION
    DRRs = fill(0.0, nsteps)u"mW/cm^2"  # total DIRECT RAYLEIGH RADIATION
    DRs = fill(0.0, nsteps)u"mW/cm^2"   # total DIRECT RADIATION
    SRs = fill(0.0, nsteps)u"mW/cm^2"   # total SCATTERED RADIATION
    Zs = fill(0.0, nsteps)u"°"   # zenith angle
    DOYs = Vector{Int}(undef, nsteps)
    TIMEs = Vector{Int}(undef, nsteps)
    step = 1
    for i in 1:ndays
        for j in 1:ntimes
            d = days[i]
            t = hours[j]
            h, tsn = hour_angle(t, lonc) # hour angle (radians)
            ζ, δ, Z, AR2 = solar_geometry(d=d, lat=lat, h=h, d0=d0, ω=ω, ϵ=ϵ, SE=SE) # compute ecliptic, declination, zenith angle and (a/r)^2
            ZZ = uconvert(°, Z)
            check_skylight(Z, NMAX, SRINT, GRINT) # checking zenith angle for possible skylight before sunrise or after sunset
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
                h, tsn = hour_angle(t) # hour angle (radians) - redundant?
                ζ, δ, Z, AR2 = solar_geometry(d=d, lat=lat, h=h, d0=d0, ω=ω, ϵ=ϵ, SE=SE) # compute ecliptic, declination, zenith angle and (a/r)^2 - redundant?
                alt = (π / 2 - Z)rad
                cazsun = (sin(lat) * sin(alt) - sin(δ)) / (cos(lat) * cos(alt)) # cos(solar azimuth)
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

                # horizon angle - check this works when starting at 0 rather than e.g. 15 deg
                azi = range(0°, stop=360° - 360° / length(hori), length=length(hori))
                ahoriz = hori[argmin(abs.(dazsun .- azi))]

                # slope zenith angle calculation (Eq. 3.15 in Sellers 1965. Physical Climatology. U. Chicago Press)
                if slope > 0°
                    czsl = cos(Z) * cos(slope) + sin(Z) * sin(slope) * cos(dazsun - aspect)
                    zsl = acos(czsl)
                    dzsl = min(uconvert(°, zsl), 90°) # cap at 90 degrees if sun is below slope horizon
                else
                    dzsl = 0°
                end

                # computing optical air mass (airms) from knowledge of true zenith angle
                # using the Rozenberg formula (see p.159 in book 'Twilight' by G.V.
                # Rozenberg, Plenum Press, New York, 1966)
                # the difference between apparent and true zenith angle is neglected
                # for z less than 88 degrees
                # variation of airms with altitude is ignored since it is negligible up to
                # at least 6 km above sea level

                # refraction correction check
                if Z < 1.5358896
                    # skip refraction correction
                else
                    refr = 16.0 + ((Z - 1.53589) * 15) / (π / 90)
                    refr = (refr / 60) * (π / 180)
                    Z -= refr
                end

                # optical air mass (Rozenberg formula) ---
                airms = 1.0 / (cos(Z) + (0.025 * exp(-11.0 * cos(Z))))
                cz = cos(Z)
                intcz = Int(floor(100.0 * cz + 1.0))
                Zdeg = uconvert(°, Z)  # zenith angle in degrees

                # ATMOSPHERIC OZONE LOOKUP
                # Convert latitude in degrees to nearest 10-degree index
                tlat = (lat + 100.0) / 10.0
                llat = Int(floor(tlat))
                allat = llat
                ala = allat + 0.5
                if tlat > ala
                    llat += 1
                end

                # Clamp llat index to valid range
                m = month(Date(year, 1, 1) + Day(d - 1)) # month from day of year
                llat = clamp(llat, 1, size(OZ, 1))
                ozone = OZ[llat, m]  # ozone thickness (cm) from lookup table

                #had = 15.0 * (t - tsn) # hour angle, degrees
                ELEVFCT1, ELEVFCT2, ELEVFCT3, ELEVFCT4 = elev_corr(elev)

                # deal with this:
                # c     mutliplier to correct hourly solar data for horizon angle
                #     if(altdeg.lt.ahoriz)then
                # c	   diffuse only - cut down to diffuse fraction      
                #     TDD(111+IT)=TDD(111+IT)* (0.12 + 0.83 * ((CCMINN(IDAY) + 
                #    &  CCMAXX(IDAY))/ 2. / 100.)) ! from Butt et al. 2010 Agricultural and Forest Meteorology 150 (2010) 361–368
                #     endif

                for N in 1:NMAX
                    τλ1 = (P / 101300u"Pa") * τR[N] * ELEVFCT1
                    τλ2 = (25u"km" / AMR) * τA[N] * ELEVFCT2
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
                    if DRλ[N] < 1.0e-24u"mW / cm^2 / nm"
                        DRλ[N] = 1.0e-24u"mW / cm^2 / nm"
                    end

                    DRRλ[N] = (Sλ[N] * AR2 * cz) * exp(-float(τλ1) * airms) / 1000.0

                    if alt < ahoriz
                        DRλ[N] = 1.0e-24u"mW / cm^2"
                        DRRλ[N] = 1.0e-24u"mW / cm^2"
                    end

                    # Sky (SRλ) and Global Radiation (GRλ)
                    if NOSCAT == 0 || N > 11
                        # skip to 400
                        SRλ[N] = 0.0u"mW / cm^2 / nm"
                    elseif IUV
                        if τλ1 >= 0.03
                            GAMR, GAML, SBAR = GAMMA(τλ1)
                            SRλ[N] = (
                                         ((float(GAML[INTCZ]) + float(GAMR[INTCZ])) / (2.0 * (1.0 - REFL * float(SBAR))))
                                         -
                                         exp(-float(τλ1) * airms)
                                     ) * CZ * Sλ[N] * AR2 / 1000.0
                        else
                            SRλ[N] = 0.0u"mW / cm^2 / nm"
                        end
                    else
                        B = Zdeg / 5.0
                        IA = floor(Int, B)
                        C = B - IA
                        I = C > 0.50 ? IA + 2 : IA + 1

                        FDAV = FD[N, I]
                        FDQDAV = FDQ[N, I]

                        SRλ[N] = (Sλ[N] / π) * (FDAV + FDQDAV * (REFL / (1.0 - (REFL * S[N])))) / 1000.0
                        SRλ[N] *= AR2

                        if N <= 7
                            EPλ[N] = (SRλ[N] + DRλ[N]) / ER[N]
                        end
                    end

                    GRλ[N] = SRλ[N] + DRλ[N]

                    if N == 1
                        SRINT[1] = 0.0u"mW / cm^2"
                        DRRINT[1] = 0.0u"mW / cm^2"
                        DRINT[1] = 0.0u"mW / cm^2"
                        GRINT[1] = 0.0u"mW / cm^2"
                        EPINT[1] = 0.0u"mW / cm^2"
                    else
                        AIλ[N] = Iλ[N]
                        AIλ[N-1] = Iλ[N-1]

                        Δλ = AIλ[N] - AIλ[N-1]

                        DRINT[N] = DRINT[N-1] + (Δλ * DRλ[N-1]) + (0.5 * Δλ * (DRλ[N] - DRλ[N-1]))
                        DRRINT[N] = DRRINT[N-1] + (Δλ * DRRλ[N-1]) + (0.5 * Δλ * (DRRλ[N] - DRRλ[N-1]))
                        SRINT[N] = SRINT[N-1] + (Δλ * SRλ[N-1]) + (0.5 * Δλ * (SRλ[N] - SRλ[N-1]))
                        GRINT[N] = GRINT[N-1] + (Δλ * GRλ[N-1]) + (0.5 * Δλ * (GRλ[N] - GRλ[N-1]))

                        if N <= 7
                            EPINT[N] = EPINT[N-1] + (Δλ * EPλ[N-1]) + (0.5 * Δλ * (EPλ[N] - EPλ[N-1]))
                            if EPINT[N] != 0.0
                                REPINT[N] = 1 / EPINT[N]
                            end
                        end
                    end
                end
            else # sunrise, sunset or long day

            end
            # Store into row `step`
            GRINTs[step, :] .= GRINT
            DRRINTs[step, :] .= DRRINT
            DRINTs[step, :] .= DRINT
            SRINTs[step, :] .= SRINT
            GRs[step] = GRINT[NMAX]
            DRRs[step] = DRRINT[NMAX]
            DRs[step] = DRINT[NMAX]
            SRs[step] = SRINT[NMAX]
            Zs[step] = ZZ
            DOYs[step] = d
            TIMEs[step] = t  # optional: start hours from 0
            Zs[Zs .> 90°] .= 90°
            step += 1
        end
    end
    return (
        λ          = Iλ,
        λDirect    = DRINTs .* (10u"W/m^2" / 1u"mW/cm^2"),
        λRayleigh  = DRRINTs .* (10u"W/m^2" / 1u"mW/cm^2"),
        λScattered = SRINTs .* (10u"W/m^2" / 1u"mW/cm^2"),
        λGlobal    = GRINTs .* (10u"W/m^2" / 1u"mW/cm^2"),
        Rayleigh   = DRRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        Direct     = DRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        Scattered  = SRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        Global     = GRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        Zenith     = Zs,
        doy        = DOYs,
        hour       = TIMEs,
    )
end