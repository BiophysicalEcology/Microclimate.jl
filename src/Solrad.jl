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
- `h`: Elevation at which to compute pressure (with length units, e.g. `u"m"`).
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
    R=Unitful.R
P_a = P_ref * (1 + (L_ref / T_ref) * h) ^ ((-g_0 * M) / (R * L_ref))

return P_a
end

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
function hour_angle(t::Real, lonc::Real = 0)
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
    d::Real = 1.,
    lat::Quantity = 83.07305u"°",
    h::Quantity = -2.87979u"rad",
    d0::Real = 80,
    ω::Real = 2π / 365,
    ϵ::Real = 0.0167238,
    se::Real = 0.39779
)
    ζ = (ω * (d - d0)) + 2ϵ * (sin(ω * d) - sin(ω * d0))          # Eq.5
    δ = asin(se * sin(ζ))                                         # Eq.4
    cosZ = cos(lat) * cos(δ) * cos(h) + sin(lat) * sin(δ)         # Eq.3
    z = acos(cosZ)u"rad"                                          # Zenith angle
    AR2 = 1 + (2ϵ) * cos(ω * d)                                   # Eq.2
    δ = δ*u"rad"
    ζ = ζ*u"rad"
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
            Skylum = (10.0 ^ Elog)*1.46E-03u"mW * cm^-2"
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

function dchxy(tau1::Float64, cfa::Vector{Float64}, nst::Int)
    chx = zeros(Float64, 101)
    chy = zeros(Float64, 101)

    if tau1 < 0.0
        return chx, chy, 1
    end

    c = zeros(Float64, 3)
    u = zeros(Float64, 3)
    c[1] = cfa[1]
    c[2] = cfa[2]
    c[3] = cfa[3]

    # Calculate initial values for root finding
    cc = 0.0
    cd = 0.0
    for n in 1:3
        cc += (2n - 1) * c[n]
        cd += (2n - 1)^2 * c[n]
    end

    # Compute roots of characteristic equation (up to 4 roots)
    q = zeros(ComplexF64, 4)
    d = zeros(ComplexF64, 4)
    ntr = 2
    q[1] = sqrt(Complex(0.25 * cc + 0.5 * sqrt(Complex(cd))))
    q[2] = sqrt(Complex(0.25 * cc - 0.5 * sqrt(Complex(cd))))
    d[1] = 1.0 / (2.0 * q[1] * (2.0 * q[1]^2 - cc))
    d[2] = 1.0 / (2.0 * q[2] * (2.0 * q[2]^2 - cc))

    # Compute exponential terms using dexpi
    e1 = dexpi(-real(q[1]) * tau1)
    e2 = dexpi(-real(q[2]) * tau1)

    # Initialize mu array (Gauss quadrature nodes)
    mu = zeros(Float64, 101)
    mu[1] = 0.0
    for i in 2:101
        mu[i] = 0.01 * (i - 1)
    end

    # Compute chx and chy using X and Y function definitions
    for i in 1:101
        sumx = 0.0 + 0im
        sumy = 0.0 + 0im
        for j in 1:ntr
            qq = q[j]
            dj = d[j]
            sumx += dj * mu[i] / (mu[i]^2 - qq^2)
            sumy += dj * qq / (mu[i]^2 - qq^2)
        end
        chx[i] = real(1.0 + 2.0 * mu[i] * sumx)
        chy[i] = real(2.0 * sumy)
    end

    return chx, chy, ntr
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
    year::Real=2001., # needed to determine if a leap year
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
    # GPL(N)       GLOBAL RADIATION BETWEEN 300 AND 320 NM IN 2 NM STEPS
    GRINT = fill(0.0, nmax)u"mW/cm^2" # GLOBAL RADIATION COMPONENT (DIRECT + SCATTERED)
    DRRINT = fill(0.0, nmax)u"mW/cm^2" # DIRECT RAYLEIGH RADIATION COMPONENT
    DRINT = fill(0.0, nmax)u"mW/cm^2" # DIRECT RADIATION COMPONENT
    SRINT = fill(0.0, nmax)u"mW/cm^2" # SCATTERED RADIATION COMPONENT
    AIλ = fill(0.0, nmax)u"nm"
    GRλ = fill(0.0, nmax)u"mW/cm^2/nm" # integrated global radiation component
    DRRλ = GRINT * u"1/nm"
    DRλ = GRINT * u"1/nm"
    SRλ = GRINT * u"1/nm"
    P = get_pressure(elev)

    ndays = length(days)
    ntimes = length(hours)
    nsteps = ndays * ntimes  # total time steps
    GRINTs = fill(0.0, nsteps, nmax)u"mW/cm^2"   # wavelength-specific GLOBAL RADIATION
    DRRINTs = fill(0.0, nsteps, nmax)u"mW/cm^2"  # wavelength-specific DIRECT RAYLEIGH RADIATION
    DRINTs = fill(0.0, nsteps, nmax)u"mW/cm^2"   # wavelength-specific DIRECT RADIATION
    SRINTs = fill(0.0, nsteps, nmax)u"mW/cm^2"   # wavelength-specific SCATTERED RADIATION
    GRs = fill(0.0, nsteps)u"mW/cm^2"   # total GLOBAL RADIATION
    DRRs = fill(0.0, nsteps)u"mW/cm^2"  # total DIRECT RAYLEIGH RADIATION
    DRs = fill(0.0, nsteps)u"mW/cm^2"   # total DIRECT RADIATION
    SRs = fill(0.0, nsteps)u"mW/cm^2"   # total SCATTERED RADIATION
    Zs = fill(0.0, nsteps)u"°"   # zenith angle
    ZSLs = fill(0.0, nsteps)u"°"   # slope zenith angle
    HHs = fill(0.0, ndays)
    tsns = fill(0.0, ndays)
    DOYs = Vector{Int}(undef, nsteps)
    times = Vector{Real}(undef, nsteps)
    step = 1
    HH = 0. # initialise sunrise hour angle
    tsn = 12. # initialise time of solar noon
    for i in 1:ndays
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
                m = month(Date(year, 1, 1) + Day(d - 1)) # month from day of year
                llat = clamp(llat, 1, size(OZ, 1))
                ozone = OZ[llat, m]  # ozone thickness (cm) from lookup table
                ELEVFCT1, ELEVFCT2, ELEVFCT3, ELEVFCT4 = elev_corr(elev)

                # deal with this:
                # c     mutliplier to correct hourly solar data for horizon angle
                #     if(altdeg.lt.ahoriz)then
                # c	   diffuse only - cut down to diffuse fraction      
                #     TDD(111+IT)=TDD(111+IT)* (0.12 + 0.83 * ((CCMINN(IDAY) + 
                #    &  CCMAXX(IDAY))/ 2. / 100.)) ! from Butt et al. 2010 Agricultural and Forest Meteorology 150 (2010) 361–368
                #     endif

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
                    if noscat == false || N > 11
                        # skip to 400
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
                        I = round(Int, (ustrip(Z) + 5) / 5)
                        FDAV = FD[N, I]
                        FDQDAV = FDQ[N, I]
                        SRλ[N] = (Sλ[N] / π) * (FDAV + FDQDAV * (refl / (1.0 - (refl * S[N])))) / 1000.0
                        SRλ[N] *= AR2
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
            GRINTs[step, :] .= GRINT
            DRRINTs[step, :] .= DRRINT
            DRINTs[step, :] .= DRINT
            SRINTs[step, :] .= SRINT
            GRs[step] = GRINT[nmax]
            DRRs[step] = DRRINT[nmax]
            DRs[step] = DRINT[nmax]
            SRs[step] = SRINT[nmax]
            Zs[step] = Z
            ZSLs[step] = Zsl
            DOYs[step] = d
            times[step] = t  # optional: start hours from 0
            Zs[Zs .> 90°] .= 90°
            step += 1
        end
        HHs[i] = HH # save today's sunrise hour angle
        tsns[i] = tsn # save today's time of sunrise
    end
    return (
        Zenith      = Zs,
        ZenithSlope = ZSLs,
        HHsr        = HHs,
        tsn         = tsns,
        doy         = DOYs,
        hour        = times,
        Rayleigh    = DRRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        Direct      = DRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        Scattered   = SRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        Global      = GRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        λ           = Iλ,
        λDirect     = DRINTs .* (10u"W/m^2" / 1u"mW/cm^2"),
        λRayleigh   = DRRINTs .* (10u"W/m^2" / 1u"mW/cm^2"),
        λScattered  = SRINTs .* (10u"W/m^2" / 1u"mW/cm^2"),
        λGlobal     = GRINTs .* (10u"W/m^2" / 1u"mW/cm^2"),
    )
end

"""
    vapour_pressure(T)

Calculates saturation vapour pressure (Pa) for a given air temperature.

# Arguments
- `T`: air temperature in K.
"""
function vapour_pressure(T)
    T = Unitful.ustrip(T)# + 273.15
    logP_vap = T
    if T <= 273.15
        logP_vap = -9.09718 * (273.16 / T - 1) - 3.56654 * log10(273.16 / T) + 0.876793 * (1 - T / 273.16) + log10(6.1071)
    else
        logP_vap = -7.90298 * (373.16 / T - 1) + 5.02808 * log10(373.16 / T) - 1.3816E-07 * (10^(11.344 * (1 - T / 373.16)) - 1) + 8.1328E-03 * (10^(-3.49149 * (373.16 / T - 1)) - 1) + log10(1013.246)
    end
    (10^logP_vap) * 100u"Pa"
end

"""
    wet_air(T_drybulb, T_wetbulb, rh, T_dew, P_atmos, fO2, fCO2, fN2)
    wet_air(T_drybulb; kw...)

Calculates several properties of humid air as output variables below. The program
is based on equations from List, R. J. 1971. Smithsonian Meteorological Tables. Smithsonian
Institution Press. Washington, DC. wet_air must be used in conjunction with function vapour_pressure.

Input variables are shown below. The user must supply known values for T_drybulb and P (P at one standard
atmosphere is 101 325 pascals). Values for the remaining variables are determined by whether the user has
either (1) psychrometric data (T_wetbulb or rh), or (2) hygrometric data (T_dew)

# TODO fix this desctiption
(1) Psychrometric data:
If T_wetbulb is known but not rh, then set rh=-1 and dp=999
If rh is known but not T_wetbulb then set T_wetbulb=0 and dp=999

(2) Hygrometric data:
If T_dew is known then set T_wetublb = 0 and rh = 0.

# Arguments

- `T_drybulb`: Dry bulb temperature (K)
- `T_wetbulb`: Wet bulb temperature (K)
- `rh`: Relative humidity (%)
- `T_dew`: Dew point temperature (K)
- `P`: Barometric pressure (Pa)
- `fO2`; fractional O2 concentration in atmosphere, -
- `fCO2`; fractional CO2 concentration in atmosphere, -
- `fN2`; fractional N2 concentration in atmosphere, -
# - `P_vap`: Vapour pressure (Pa)
# - `P_vap_sat`: Saturation vapour pressure (Pa)
# - `ρ_vap`: Vapour density (kg m-3)
# - `r_w Mixing`: ratio (kg kg-1)
# - `T_vir`: Virtual temperature (K)
# - `T_vinc`: Virtual temperature increment (K)
# - `ρ_air`: Density of the air (kg m-3)
# - `cp`: Specific heat of air at constant pressure (J kg-1 K-1)
# - `ψ`: Water potential (Pa)
# - `rh`: Relative humidity (%)

"""
function wet_air(T_drybulb; 
    T_wetbulb=T_drybulb, 
    rh=0, 
    T_dew=nothing, 
    P_atmos=101325u"Pa",
    fO2=0.2095,
    fCO2=0.0004,
    fN2=0.79
)
    return wet_air(T_drybulb, T_wetbulb, rh, T_dew, P_atmos, fO2, fCO2, fN2)
end
function wet_air(T_drybulb, T_wetbulb=T_drybulb, rh=0, T_dew=nothing, P_atmos=101325u"Pa", fO2 = 0.2095, fCO2 = 0.0004, fN2 = 0.79)
    f_w = 1.0053 # (-) correction factor for the departure of the mixture of air and water vapour from ideal gas laws
    M_w = (1molH₂O |> u"kg")/1u"mol" # molar mass of water
    M_a = ((fO2*molO₂ + fCO2*molCO₂ + fN2*molN₂) |> u"kg")/1u"mol" # molar mass of air
    P_vap_sat = vapour_pressure(T_drybulb)
    if T_dew !== nothing
        P_vap = vapour_pressure(T_dew)
        rh = (P_vap / P_vap_sat) * 100
    else
        if rh !== nothing
            P_vap = P_vap_sat * rh / 100
        else
            δ_bulb = T_drybulb - T_wetbulb
            δ_P_vap = (0.000660 * (1 + 0.00115 * (Unitful.ustrip(T_wetbulb)-273.15)) * Unitful.ustrip(P) * Unitful.ustrip(δ_bulb))u"Pa"
            P_vap = vapour_pressure(T_wetbulb) - δ_P_vap
            relhumid = (P_vap / P_vap_sat) * 100
        end
    end
    r_w = ((0.62197 * f_w * P_vap) / (P_atmos - f_w * P_vap))u"kg/kg"
    ρ_vap = P_vap * M_w / (0.998 * Unitful.R * T_drybulb)
    ρ_vap = Unitful.uconvert(u"kg/m^3",ρ_vap) # simplify units
    T_vir = T_drybulb * ((1.0 + r_w / (18.016 / 28.966)) / (1 + r_w))
    T_vinc = T_vir - T_drybulb
    ρ_air = (M_a / Unitful.R) * P_atmos / (0.999 * T_vir)
    ρ_air = Unitful.uconvert(u"kg/m^3",ρ_air) # simplify units
    cp = ((1004.84 + (r_w * 1846.40)) / (1 + r_w))u"J/K/kg"
    ψ = if min(rh) <= 0
        -999u"Pa"
    else
        (4.615e+5 * Unitful.ustrip(T_drybulb) * log(rh / 100))u"Pa"
    end

    return (;P_vap, P_vap_sat, ρ_vap, r_w, T_vinc, ρ_air, cp, ψ, rh)
end

"""
    dry_air(T_drybulb; kw...)
    dry_air(T_drybulb, P_atmos, elev, fO2, fCO2, fN2)

"""
dry_air(T_drybulb; P_atmos=nothing, elev=0m, fO2 = 0.2095, fCO2 = 0.0004, fN2 = 0.79) = dry_air(T_drybulb, P_atmos, elev, fO2, fCO2, fN2)
function dry_air(T_drybulb, P_atmos, elev, fO2, fCO2, fN2)
    σ = Unitful.uconvert(u"W/m^2/K^4",Unitful.σ) # Stefan-Boltzmann constant, W/m^2/K^4, extract σ when calling Unitful when units issue is fixed in Unitful
    M_a = ((fO2*molO₂ + fCO2*molCO₂ + fN2*molN₂) |> u"kg")/1u"mol" # molar mass of air
    if isnothing(P_atmos)
        P_atmos = get_pressure(elev)#P_std * ((1 - (0.0065 * elev / 288m))^(1 / 0.190284))
    end
    ρ_air = (M_a / Unitful.R) * P_atmos / (T_drybulb)
    ρ_air = Unitful.uconvert(u"kg/m^3",ρ_air) # simplify units
    vis_not = 1.8325e-5u"kg/m/s"
    T_not = 296.16u"K"
    c = 120u"K"
    μ = (vis_not * (T_not + c) / (T_drybulb + c)) * (T_drybulb / T_not)^1.5 # kg / m.s
    ν = μ / ρ_air # m2 / s or J.s/kg
    dif_vpr = 2.26e-5u"m^2/s" * ((T_drybulb / 273.15u"K")^1.81) * (1.e5u"Pa" / P_atmos) # m2 / s
    k_fluid = (0.02425 + (7.038e-5 * (Unitful.ustrip(T_drybulb) - 273.15)))u"W/m/K"
    L_v = (2.5012e6 - 2.3787e3 * (Unitful.ustrip(T_drybulb) - 273.15))u"J/kg"

    tcoeff = 1 / T_drybulb
    ggroup = 0.0980616u"m/s^2" * tcoeff / (ν^2) # 1 / m3.K
    bbemit = σ * ((T_drybulb)^4) # W/m2
    emtmax = 2.897e-3u"K*m" / (T_drybulb) # m

    return (;P_atmos, ρ_air, μ, ν, dif_vpr, k_fluid, L_v, tcoeff, ggroup, bbemit, emtmax)
end

function get_profile(;
    refhyt=1.2u"m",
    ruf=0.004u"m",
    zh=0.0u"m",
    d0=0.0u"m",
    κ = 0.4, # Kármán constant
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
        return (1 - min(1, GAM * z / AMOL))^0.25
    end

    function PSI1(X)
        return 2 * log((1 + X) / 2) + log((1 + X^2) / 2) - 2 * atan(X) + π / 2
    end

    function PSI2(X)
        return 2 * log((1 + X^2) / 2)
    end

    function get_Obukhov(TA, TS, VEL, z, z0, rcptkg, κ)

        AMOL = -30.0u"cm" # initial Monin-Obukhov length cm
        GAM = 16.0 # -
        #RCPTKG = 6.003e-8u"cal/minute/cm/K" #CAL-MIN-CM-C
        z = u"cm"(z)
        z0 = u"cm"(z0)
        ZRATIO = z / z0 + 1
        DUM = log(ZRATIO)
        DIFFT = TA - TS
        TAVE = (TA + TS) / 2.0
        RCP = RHOCP(TAVE)
        DEL = 1.0
        count = 0
        USTAR = 0.0
        QC = 0.0
        STO = 0.0
        STB = 0.0
        STS = 0.0

        while DEL > 1e-2 && count < 500
            count += 1
            X = PHI(z, GAM, AMOL)
            Y = PSI1(X)
            YY = PSI2(X)
            USTAR = κ * (VEL * 100 * 60) / (log(z / z0) - Y)

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
    rcptkg = u"cal*minute^2/m^4"(ρ * cp * TREF / (κ * g))

    GAM = 16
    ZRATIO = z / z0 + 1.0
    DUM = log(ZRATIO)
    USTAR = κ * V / DUM
    DIFFT = T1 - T3
    TAVE = (T3 + T1) / 2
    RCP = RHOCP(TAVE, elev, rh)
    AMOL = -30.0
    if zh > 0.0u"m"
        STS = 0.62 / (ustrip(z0) * ustrip(USTAR) / 12)^0.45
        STB = 0.64 / DUM
        QC = RCP * DIFFT * USTAR * STB / (1.0 + STB / STS)

        for i in 2:NAIR
            if T1 ≥ T3 || T3 ≤ maxsurf || ZEN ≥ 90
                VV[i] = (USTAR / κ) * log(ZZ[i] / z0 + 1)
            else
                X1 = PHI(ZZ[i])
                Y1 = PSI1(X1)
                ADUM = ZZ[i] / z0 - Y1
                VV[i] = (USTAR / κ) * log(ADUM)
            end

            A = (T1 - T3) / (1 - log((z - d0_cm) / zh_cm))
            T0 = T1 + A * log((z - d0_cm) / zh_cm)
            T[i] = T0 - A * log((ZZ[i] - d0_cm) / zh_cm)
        end
    else
        if T1 ≥ T3 || T3 ≤ u"K"(maxsurf) || ZEN ≥ 90
            STS = 0.62 / (ustrip(z0) * ustrip(USTAR) / 12)^0.45
            STB = 0.64 / DUM
            QC = RCP * DIFFT * USTAR * STB / (1.0 + STB / STS)

            for i in 2:NAIR
                VV[i] = (USTAR / κ) * log(ZZ[i] / z0 + 1.0)
                TZO = (T1 * STB + T3 * STS) / (STB + STS)
                T[i] = TZO + (T1 - TZO) * log(ZZ[i] / z0 + 1.0) / DUM
            end
        else
            for i in 2:NAIR
                X1 = PHI(ZZ[i])
                Y1 = PSI1(X1)
                YY2 = PSI2(X1)
                X = PHI(z)
                Y = PSI1(X)
                YY = PSI2(X)
                ADUM = ZZ[i] / z0 - Y1
                VV[i] = (USTAR / κ) * log(ADUM)

                Obukhov_out = get_Obukhov(T1, T3, V, ZZ[i], z0, RCP, κ)
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
        heights=heights,
        VELs=u"m/s".(VV),      # m/s
        TAs=u"°C".(T),         # deg C
        RHs=RHs,               # %
        QCONV=u"W/m^2"(QC),    # W
        USTAR=u"m/s"(USTAR)    # m/s
    )
end

function sinec!(
    ITAIR, TIMARY, TAIRRY, TMIN, TMAX, TMIN2, TMAX2, TIMSR, TIMSS, TIMTMX, daily, iday
)
    TMIN = u"K"(TMIN)
    TMAX = u"K"(TMAX)
    TMIN2 = u"K"(TMIN2)
    TMAX2 = u"K"(TMAX2)
    T=TMIN[1] # initialise T
    A = (TMAX - TMIN) / 2
    TSR = TMIN
    TREF = (TIMTMX - TIMSR) / 2 + TIMSR
    SS = 360.0 * (TIMSS - TREF) / (2.0 * (TIMTMX - TIMSR))
    SY = SS / 57.29577
    ZS = sin(SY)
    TSS = A * ZS + TMIN + A
    TAU = 3.0 / ((2400.0 - TIMSS) + TIMSR)

    for I in 1:24
        J = I + 1
        TIME = I * 100.0
        if TIME <= TIMSR
            T = TMIN
        elseif TIME >= TIMSS
            TI = (2400.0 - TIMSS) - (2400.0 - TIME)
            E = TI * TAU
            if daily
                T = (TSS - TMIN2) / exp(E) + TMIN2
            else
                T = (TSS - TSR) / exp(E) + TSR
            end
        elseif TIME > TIMSR && TIME < TIMSS
            X = 360.0 * (TIME - TREF) / (2.0 * (TIMTMX - TIMSR))
            Y = X / 57.29577
            Z = sin(Y)
            T = A * Z + TMIN + A
        else
            TI = (2400.0 - TIMSS) + TIME
            if daily && iday > 1
                T = ((TMIN - 273.0 - ITAIR) / TIMSR) * TIME + ITAIR
            else
                E = TI * TAU
                T = (TSS - TSR) / exp(E) + TSR
            end
        end

        T = u"°C"(T)
        ITIME = Int(floor(TIME / 100.0))
        FRMIN = TIME / 100.0 - ITIME
        ITIME *= 60
        FRMIN *= 100.0
        TIMEC = ITIME + FRMIN

        TIMARY[J] = TIMEC
        TAIRRY[J] = T
    end

    # Set first time step
    TIMARY[1] = 0.0
    if daily == 1 && iday > 1
        TAIRRY[1] = ITAIR
    else
        TAIRRY[1] = T
    end
    ITAIR = TAIRRY[25]
end

function vsine(VMIN, VMAX, TIMSR, TIMSS, TIMIN, TIMAX, daily, iday, IVAR)
    vinit = 0.0 * VMIN[1]
    XA = fill(0.0, 25)
    YA = fill(vinit, 25)
    TIMARY = fill(0.0, 25)  
    vave = (VMAX + VMIN) / 2.0

    if daily == 1 && iday > 1
        vinit = IVAR
    end

    vsmIN = VMIN + 0.01 * vave
    vsmAX = VMAX - 0.01 * vave

    if TIMIN < TIMAX
        # morning min, afternoon max
        ITEST1 = Int(floor(TIMIN / 100))
        ITEST2 = Int(floor(TIMAX / 100))
        slope1 = (vave - vsmIN) / (100.0 - TIMIN)
        slope2 = (vsmIN - vsmAX) / (TIMIN - TIMAX)
        slope3 = (vsmAX - vave) / (TIMAX - 2400.0)
    else
        # morning max, afternoon min
        ITEST1 = Int(floor(TIMAX / 100))
        ITEST2 = Int(floor(TIMIN / 100))
        slope1 = (vave - vsmAX) / (100.0 - TIMAX)
        slope2 = (vsmAX - vsmIN) / (TIMAX - TIMIN)
        slope3 = (vsmIN - vave) / (TIMIN - 2400.0)
    end


    for I in 1:25
        XA[I] = I * 100.0 - 100.0
        TIME = XA[I]

        ITIME = Int(floor(TIME / 100.0))
        FRMIN = (TIME / 100.0) - ITIME
        ITIME *= 60
        FRMIN *= 100.0
        TIMEC = ITIME + FRMIN
        TIMARY[I] = TIMEC

        if I == 1
            YA[I] = (daily == 1 && iday > 1) ? vinit : vave
            continue
        end

        if I < ITEST1
            YA[I] = YA[1] - slope1 * (XA[1] - XA[I])
            continue
        end

        if I == ITEST1
            YA[I] = (TIMIN < TIMAX) ? vsmIN : vsmAX
            continue
        end

        if I < ITEST2
            if TIMIN < TIMAX
                YA[I] = vsmIN - slope2 * (TIMIN - XA[I])
            else
                YA[I] = vsmAX - slope2 * (TIMAX - XA[I])
                if YA[I] > vsmAX
                    YA[I] = vsmAX
                end
            end
            continue
        end

        if I == ITEST2
            YA[I] = (TIMIN < TIMAX) ? vsmAX : vsmIN
            continue
        end

        if I > ITEST2
            if TIMIN < TIMAX
                YA[I] = vsmAX - slope3 * (TIMAX - XA[I])
            else
                YA[I] = vsmIN - slope3 * (TIMIN - XA[I])
            end
        end
    end

    # Fix any negative values
    for JCT in 1:25
        if YA[JCT] < 0.0 * VMIN[1]
            if JCT < 25 && YA[JCT + 1] > 0.0 * VMIN[1]
                YA[JCT] = (YA[JCT - 1] + YA[JCT + 1]) / 2.0
            else
                YA[JCT] = abs(YA[JCT])
            end
        end
    end

    return YA
end

function hourly_vars(;
    TMINN::Vector,
    TMAXX::Vector,
    WNMINN::Vector,
    WNMAXX::Vector,
    RHMINN::Vector,
    RHMAXX::Vector,
    CCMINN::Vector,
    CCMAXX::Vector,
    solrad_out::Any,
    TIMINS::Vector=[0, 0, 1, 1],
    TIMAXS::Vector=[1, 1, 0, 0],
    daily::Bool=false)

    ndays = length(TMINN)
    Tairs = fill(TMINN[1], 25 * ndays)
    Clds = fill(CCMINN[1], 25 * ndays)
    RHs = fill(RHMINN[1], 25 * ndays)
    WNs = fill(WNMINN[1], 25 * ndays)

    ITAIR = TMINN[1] # initial air temperature for daily
    IWN = WNMINN[1]
    IRH = RHMAXX[1]
    ICLD = CCMINN[1]

    for iday in 1:ndays
        TIMARY = fill(0.0, 25)
        TAIRRY = fill(0.0u"°C", 25)
        WNARRY = fill(IWN, 25)
        RHARRY = fill(IRH, 25)
        CLDARRY = fill(ICLD, 25)

        HH = solrad_out.HHsr[iday]
        tsn = solrad_out.tsn[iday]

        #     Air temperature calculations
        #     SUNSET IN MILITARY TIME
        HHINT = trunc(HH)
        FRACT = (HH - HHINT) * 60.0
        HSINT = trunc(tsn)
        FRACTS = (tsn - HSINT) * 60.0
        TIMSS = (HSINT * 100.0 + FRACTS) + (HHINT * 100.0 + FRACT)
        #     SUNRISE IN MILITARY TIME
        DELT = tsn - HH
        HRINT = trunc(DELT)
        FRACTR = (DELT - HRINT) * 60.0
        #     TIME OF AIR TEMPERATURE MINIMUM (NOTE: A KLUGE OF 200, I.E. 2
        #     HOURS IS BEING USED TO MAKE SINEC WORK WHEN MORNING MINIMUM
        #     IS NOT AT TRUE SUNRISE, THE ALGORITHM IS 2 HOURS OFF FOR SOME
        #     UNKNOWN (6/7/89) REASON)
        if TIMINS[1] > 0
            TSRHR = TIMINS[1] + 2.0
        else
            TSRHR = TIMINS[1]
        end
        TIMSR = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        #     TIME OF AIR TEMPERATURE MAXIMUM
        TSNHR = TIMAXS[1]# + TIMCOR (TIMCOR never used in Fortran version - uninitialised?)
        TIMTMX = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        #     SETTING TMIN, TMAX FROM ARRAYS OBTAINED FROM SUB. IOSOLR
        TMIN = TMINN[iday]
        TMAX = TMAXX[iday]
        if iday < ndays & ndays > 1
            TMIN2 = TMINN[iday+1]
            TMAX2 = TMAXX[iday+1]
        else
            TMIN2 = TMIN
            TMAX2 = TMAX
        end

        #     SETTING TIME OF MIN & MAX (HOURS BEFORE SUNRISE OR
        #     AFTER SOLAR NOON)
        TIMIN = TIMSR
        TIMAX = TIMTMX
        sinec!(ITAIR, TIMARY, TAIRRY, TMIN, TMAX, TMIN2, TMAX2, TIMSR, TIMSS, TIMTMX, daily, iday)

        # wind speed
        VMIN = WNMINN[iday]
        VMAX = WNMAXX[iday]
        #      SETTING MAX & MIN TIMES RELATIVE TO SUNRISE & SOLAR NOON
        #     TIME OF MINIMUM
        TSRHR = TIMINS[2]
        TIMSR = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        #      TIME OF MAXIMUM at sunrise for relative humidity
        TSNHR = TIMAXS[2] #+ TIMCOR
        TIMTMX = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        TIMIN = TIMTMX
        TIMAX = TIMSR
        IVAR = IWN
        WNARRY = vsine(VMIN, VMAX, TIMSR, TIMSS, TIMIN, TIMAX, daily, iday, IVAR)

        # relative humidity
        VMIN = RHMINN[iday]
        VMAX = RHMAXX[iday]
        #      SETTING MAX & MIN TIMES RELATIVE TO SUNRISE & SOLAR NOON
        #     TIME OF MINIMUM
        TSRHR = TIMINS[3]
        TIMSR = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        #      TIME OF MAXIMUM at sunrise for relative humidity
        TSNHR = TIMAXS[3] #+ TIMCOR
        TIMTMX = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        TIMIN = TIMTMX
        TIMAX = TIMSR
        IVAR = IRH
        RHARRY = vsine(VMIN, VMAX, TIMSR, TIMSS, TIMIN, TIMAX, daily, iday, IVAR)

        # cloud cover
        VMIN = CCMINN[iday]
        VMAX = CCMAXX[iday]
        #      SETTING MAX & MIN TIMES RELATIVE TO SUNRISE & SOLAR NOON
        #     TIME OF MINIMUM
        TSRHR = TIMINS[4]
        TIMSR = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        #      TIME OF MAXIMUM at sunrise for relative humidity
        TSNHR = TIMAXS[4] #+ TIMCOR
        TIMTMX = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        TIMIN = TIMTMX
        TIMAX = TIMSR
        IVAR = ICLD
        CLDARRY = vsine(VMIN, VMAX, TIMSR, TIMSS, TIMIN, TIMAX, daily, iday, IVAR)

        Tairs[(iday*25-24):(iday*25)] = TAIRRY[1:25]
        WNs[(iday*25-24):(iday*25)] = WNARRY[1:25]
        RHs[(iday*25-24):(iday*25)] = RHARRY[1:25]
        Clds[(iday*25-24):(iday*25)] = CLDARRY[1:25]
    end
    return Tairs, WNs, RHs, Clds
end

# Julia translation of FORTRAN code that calculates variable thermal conductivity and specific heat
# for soil and snow layers, based on Campbell et al. (1994) and Campbell & Norman (1998)

function soil_properties(T_soil, θ_soil, node, soilprops, numtyps, elev)
    NON = length(node)
    runmoist = false # to do
    runsnow = false # to do
    θ_sat = soilprops[:, 2] # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    λ_m = soilprops[:, 3] # mineral thermal conductivity, W/m/K
    cp_m = soilprops[:, 4] # mineral specific heat capacity, J/kg/K
    ρ_m = u"kg/m^3".(soilprops[:, 5]) # mineral density, kg/m^3
    ρ_dry = u"kg/m^3".(soilprops[:, 1]) # dry density, kg/m^3
    λ_b = fill(λ_m[1], NON) # bulk thermal conductivity, W/m/K
    cp_b = fill(cp_m[1], NON) # bulk specific heat capacity, J/kg/K
    ρ_b = fill(ρ_m[1], NON) # bulk density, kg/m^3

    ii, ij = runsnow ? (9, 8) : (1, 0)
    j = 1
    for i in ii:NON # ii ensures that it starts at soil depth if snow on top 
        if i >= node[j] + ij
            j = min(j + 1, numtyps)
        end
        # don't make it volumetric, but rather mass-specific, so don't multiply by kg/m3 converters (constants from Campbell and Norman 1998, Table 8.2)
        if runmoist
            m = runsnow ? θ_soil[i-8] : θ_soil[i]
            cp_b[i] = ρ_dry[j] / ρ_m[j] * cp_m[j] + m * 4184.0u"J/kg/K"
        else
            cp_b[i] = ρ_dry[j] / ρ_m[j] * cp_m[j] + θ_sat[j] * θ_soil[j] * 4184.0u"J/kg/K"
        end
        # constants from Campbell and Norman 1998, Table 8.2
        if runmoist
            m = runsnow ? θ_soil[i-8] : θ_soil[i]
            ρ_b[i] = m * 1000.0u"kg/m^3" + ρ_dry[j]
        else
            ρ_b[i] = θ_sat[j] * θ_soil[j] * 1000.0u"kg/m^3" + ρ_dry[j]
        end

        p_a0 = 101325.0u"Pa" # standard sea level air pressure, Pa
        p_a = get_pressure(elev)

        T_K = T_soil[i]
        T_C = u"°C"(T_K)
        λ_a = (0.024 + 7.73e-5 * ustrip(T_C) - 2.6e-8 * ustrip(T_C)^2)u"W/m/K" # thermal conductivity of dry air (equation 9 in Campbell et al. 1994)
        λ_w = (0.554 + 2.24e-3 * ustrip(T_C) - 9.87e-6 * ustrip(T_C)^2)u"W/m/K" # thermal conductivity of water (equation 8 in Campbell et al. 1994)

        hr = 1.0 # relative humidity
        D_v0 = 2.12e-5u"m^2/s" # vapour diffusivity in air (m2/s), standard value at 0 deg C and sea level pressure (Campbell et al. 1994)
        ρ_hat0 = 44.65u"mol/m^3" # molar density of air (mol/m3), standard value at 0 deg C and sea level pressure (Campbell et al. 1994)

        D_v = D_v0 * (p_a0 / p_a) * (T_K / 273.15u"K")^1.75 # temperature/pressure-corrected vapour diffusivity in air (m2/s) (p. 309 in Campbell et al. 1994)
        ρ_hat = ρ_hat0 * (p_a / p_a0) * (273.15 / T_K) # temperature/pressure-corrected molar density of air (mol/m3) (p. 309 in Campbell et al. 1994)
        λ_vap = (45144.0 - 48.0 * ustrip(T_C))u"J/mol" # J/mol latent heat of vaporization (Cambell et al. 1994, p. 309)

        rh = 99.0
        e_a = wet_air(T_K; rh=rh, P_atmos = p_a).P_vap
        e_a1 = wet_air(T_K-1u"K"; rh=rh, P_atmos = p_a).P_vap
        e_a2 = wet_air(T_K+1u"K"; rh=rh, P_atmos = p_a).P_vap

        ∇x = (e_a2 - e_a1) / 2.0 # slope of the vapour pressure function centred at focal temperature

        # these could vary with soil texture but the relationship isn't strong
        # power for liquid recirculation, mean in Table 2 of Campell et al. 1994, excluding peat moss value
        # q_0=4
        # q=q_0*(T_L/303.0u"K")^2
        q = 4.0 # using a typical value for 'q' of 4 - program becomes unstable if this is temperature dependent
        θ_0 = 0.162 # mean in Table 2 of Campell et al. 1994, excluding peat moss value

        ϕ_m = ρ_dry[j] / ρ_m[j] # volume fraction of minerals
        θ = runmoist == 1 ? (runsnow == 1 ? θ_soil[i-8] : θ_soil[i]) : θ_soil[j] * θ_sat[j] # volume fraction of water

        ϕ_g = max(0.0, 1.0 - θ - ϕ_m) # volume fraction of gas

        f_w = 1.0 / (1.0 + (θ / θ_0)^(-q)) # eq 8.17 Campbell and Norman 1988, using temperature-specific q

        λ_g = λ_a + λ_vap * ∇x * hr * f_w * ρ_hat * D_v / (p_a - e_a) # eq 8.18 Campbell and Norman 1988
        λ_f = λ_g + f_w * (λ_w - λ_g) # eq 8.19 Campbell and Norman 1988

        g_a = 0.1 # 0.1 for mineral soils, 0.33 for organic, p 125, Campbell and Norman 1988
        g_c = 1.0 - 2.0 * g_a # p 125, Campbell and Norman

        # equation 8.20 in Campbell and Norman 1988
        ϵ_g = 2.0 / (3.0 * (1.0 + g_a * (λ_g / λ_f - 1.0))) + 1.0 / (3.0 * (1.0 + g_c * (λ_g / λ_f - 1.0))) 
        ϵ_w = 2.0 / (3.0 * (1.0 + g_a * (λ_w / λ_f - 1.0))) + 1.0 / (3.0 * (1.0 + g_c * (λ_w / λ_f - 1.0)))
        ϵ_m = 2.0 / (3.0 * (1.0 + g_a * (λ_m[j] / λ_f - 1.0))) + 1.0 / (3.0 * (1.0 + g_c * (λ_m[j] / λ_f - 1.0)))
        
        # equation 8.13 in Campbell and Norman 1988
        λ_b[i] = (θ * ϵ_w * λ_w + ϕ_m * ϵ_m * λ_m[j] + ϕ_g * ϵ_g * λ_g) /
                       (θ * ϵ_w + ϕ_m * ϵ_m + ϕ_g * ϵ_g)

        #convert to cal/min/g
        #λ_b[i] = λ_b[i] / 418.6 * 60.0
        #cp_b[i] /= 4186.0
        #ρ_b[i] /= 1000.0
    end

    #HTOFN = 333500.0  # J/kg
    # if runsnow
    #     snowdens = 0.0
    #     if densfun[1] > 0
    #         if densfun[3] > 0
    #             snowdens = (densfun[1] - densfun[2]) * (1.0 - exp(-densfun[3] * cursnow - densfun[4] * snowage)) + densfun[2]
    #         else
    #             snowdens = min(0.9167, densfun[1] * snowage + densfun[2])
    #         end
    #     end

    #     if cursnow >= minsnow
    #         e_a = wet_air(T_L; rh=100, P_atmos = p_a).r_w
    #         cpsnow = (2100.0 * snowdens + (1.005 + 1.82 * (RW / (1.0 + RW))) * 1000.0 * (1.0 - snowdens))
    #         snowcond2 = (0.00395 + 0.00084 * snowdens * 1000.0 - 0.0000017756 * (snowdens * 1000.0)^2 +
    #                      0.00000000380635 * (snowdens * 1000.0)^3) / 418.6 * 60.0

    #         for i in 1:8
    #             λ_b[i] = snowcond > 0.0 ? snowcond : snowcond2
    #             cp_b[i] = (T_C > -0.45 && T_C <= 0.4) ? (cpsnow + HTOFN) / 4186.0 : cpsnow / 4186.0
    #             ρ_b[i] = (T_C > 0.0 || (snode[i] < 1e-8 && snode[min(8, i + 1)] > 0.0)) ? snowdens : 0.0
    #         end
    #     else
    #         for i in 1:8
    #             λ_b[i] = λ_b[9]
    #             cp_b[i] = cp_b[9]
    #             ρ_b[i] = 0.0
    #         end
    #     end
    # end

    return λ_b, cp_b, ρ_b
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

Base.@kwdef struct MicroParams
    soilprops::Matrix{Union{Unitful.AbstractQuantity, Float64}}
    dep::Vector{<:Unitful.AbstractQuantity}
    refhyt::Quantity
    ruf::Quantity
    d0::Quantity
    zh::Quantity
    slope::Quantity
    shade::Float64
    viewf::Float64
    elev::Quantity
    refl::Float64
    sle::Float64
    slep::Float64 # fix at 1?
    pctwet::Float64
    tdeep::Quantity
    θ_soil::Vector{Float64}
    runmoist::Bool
    runsnow::Bool
end
function soil_energy_balance!(dT, T, p::MicroParams, t)
    # extract parameters
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
    θ_soil = p.θ_soil # parameter for now

    N = length(dep)
    #dT = fill(0.0u"K/minute", N)
    #dT .= (0.0u"K/minute")

    # get soil properties and convert to cal/cm/g/C
    λ_b, cp_b, ρ_b = soil_properties(T, θ_soil, node, soilprops, numtyps, elev)
    λ_b = u"cal/cm/K/minute".(λ_b)
    cp_b = u"cal/g/K".(cp_b)
    ρ_b = u"g/cm^3".(ρ_b)

    # Get environmental data at time t
    tair = TAIRt(ustrip(t))
    vel = VELt(ustrip(t))
    zenr = ZENRt(ustrip(t))
    solr = u"cal/cm^2/minute"(max(0.0u"W/m^2", SOLRt(ustrip(t))))
    cloud = CLDt(ustrip(t))
    rh = RHt(ustrip(t))
    zslr = ZSLt(ustrip(t))

    T[N] = tdeep # set boundary condition of deep soil temperature

    depp = fill(0.0u"cm", N + 1)
    depp[1:10] = dep
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
    arad = ustrip(1.72 * (ustrip(u"kPa"(P_vap))/ustrip(u"K"(tair))) ^ (1.0/7.0)) * σ * (u"K"(tair)) ^ 4 # Campbell and Norman 1998 eq. 10.10 to get emissivity of sky
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

    # Convection
    profile_out = get_profile(
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
    dT[1] = (qsolar + qrad + qconv - qevap) / wc[1]

    # Soil conduction for internal nodes
    for i in 2:N-1
        dT[i] = (c[i-1] * (T[i-1] - T[i]) + c[i] * (T[i+1] - T[i])) / wc[i]
    end

    # Lower boundary condition
    dT[N] = 0.0u"K/minute"  # or set T[N] = T_surface from data
end