# Integer arrays TODO give them real names
const NC0 = reshape(Int[
        3, 4, 1, 2,
        2, 4, 1, 3,
        2, 3, 1, 4,
        1, 4, 2, 3,
        1, 3, 2, 4,
        1, 2, 3, 4,
    ], 4, 6)

const NC1 = reshape(Int[
        2, 3, 4, 1,
        1, 3, 4, 2,
        1, 2, 4, 3,
        1, 2, 3, 4,
        4, 1, 2, 3,
        3, 1, 2, 4,
        2, 1, 3, 4,
        1, 2, 3, 4,
    ], 4, 8)

"""
    hour_angle(t::Quantity, longitude_correction::Quantity) -> Quantity

Compute the solar hour angle `h` in radians.

# Arguments
- `t`: Local solar hour (e.g., `14.0`)
- `longitude_correction`: Longitude correction in hours (e.g., `0.5`)

# Returns

- Hour angle `h` as a `Quantity` in radians
- Time at solar noon, `tsn` as a time in hours

# Reference
McCullough & Porter 1971, Eq. 6
"""
function hour_angle(t::Real, longitude_correction::Real=0)
    tsn = 12.0 + longitude_correction                      # solar noon time
    h = (π / 12) * (t - tsn) * u"rad"      # convert hours to radians
    return h, tsn
end

abstract type AbstractSolarGeometryModel end

"""
    solar_geometry(d::Real, latitude::Quantity, h::Quantity; d0::Real = 80, ω::Real = 2π/365, ϵ::Real = 0.0167, se::Real = 0.39779)

Computes key solar geometry parameters based on McCullough & Porter (1971):

- `ζ`: Auxiliary solar longitude (radians)
- `δ`: Solar declination (radians)
- `z`: Solar zenith angle (radians)
- `AR2`: Square of Earth-to-Sun radius factor (unitless)

# Arguments
- `d`: Day of year (1–365)
- `latitude`: Latitude (with angle units, e.g. `u"°"` or `u"rad"`)
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
@kwdef struct McCulloughPorterSolarGeometry <: AbstractSolarGeometryModel
    d0::Real = 80
    ω::Real = 2π / 365
    ϵ::Real = 0.0167238
    se::Real = 0.39784993 #0.39779
end

solar_geometry(::McCulloughPorterSolarGeometry, ::Missing, ; kwargs...) = missing
function solar_geometry(sm::McCulloughPorterSolarGeometry, latitude::Quantity; # =83.07305u"°",
    d::Real, # =1.0,
    h::Quantity, # =-2.87979u"rad",
)
    (; d0, ω, ϵ, se) = sm

    ζ = (ω * (d - d0)) + 2.0ϵ * (sin(ω * d) - sin(ω * d0))          # Eq.5
    δ = asin(se * sin(ζ))                                         # Eq.4
    cosZ = cos(latitude) * cos(δ) * cos(h) + sin(latitude) * sin(δ)         # Eq.3
    z = acos(cosZ)u"rad"                                          # Zenith angle
    AR2 = 1.0 + (2.0ϵ) * cos(ω * d)                                   # Eq.2
    δ = δ * u"rad"
    ζ = ζ * u"rad"
    return(; ζ, δ, z, AR2)
end

"""
    elevation_correction(elevation)

Calculates smooth polynomial approximations of atmospheric constituent correction factors 
as a function of elevation (based on Kearney's modification of the ALTFCT array originally 
from SOLAR.DAT). Input `elevation` is the elevation in meters and can include units.

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
from `elevation` (in meters) using continuous approximation.

# Returns

A `NamedTuple` with the fields:
- `molecular_corr`
- `aerosol_corr`
- `ozone_corr`
- `water_vapor_corr`
"""
function elevation_correction(elevation)
    elev_km = ustrip(u"km", elevation + 1.0u"km")

    molecular_corr =0.00007277 * elev_km^3 +
               0.00507293 * elev_km^2 -
               0.12482149 * elev_km +
               1.11687469

    aerosol_corr =  8.35656e-7 * elev_km^6 -
               6.26384e-5 * elev_km^5 +
               1.86967e-3 * elev_km^4 -
               2.82585e-2 * elev_km^3 +
               2.26739e-1 * elev_km^2 -
               9.25268e-1 * elev_km +
               1.71321

    ozone_corr =    1.07573e-6 * elev_km^5 -
               5.14511e-5 * elev_km^4 +
               7.97960e-4 * elev_km^3 -
               4.90904e-3 * elev_km^2 +
               2.99258e-3 * elev_km +
               1.00238

    water_vapour_corr = 1.0

    return (; molecular_corr, aerosol_corr, ozone_corr, water_vapour_corr)
end

# function GAMMA(TAU1::Float64)

#     CHX = zeros(Float64, 101)
#     CHY = zeros(Float64, 101)
#     CFA = zeros(Float64, 3)
#     AMU = zeros(Float64, 101)
#     X1 = zeros(Float64, 101)
#     Y1 = zeros(Float64, 101)
#     X2 = zeros(Float64, 101)
#     Y2 = zeros(Float64, 101)
#     AIL = zeros(Float64, 101)
#     AI = zeros(Float64, 30)
#     XA = zeros(Float64, 4)
#     XB = zeros(Float64, 8)
#     GAMR = zeros(Float64, 101)
#     GAML = zeros(Float64, 101)
#     # Set up AMU array
#     AMU[1] = 0.0
#     for I in 2:101
#         AMU[I] = 0.01 * (I - 1)
#     end

#     # Compute X1, Y1 using dchxy
#     CFA[1] = 0.75
#     CFA[2] = -0.75
#     CFA[3] = 0.0
#     NST = 111
#     CHX, CHY, NTR = dchxy(TAU1, CFA, NST)
#     for I in 1:101
#         X1[I] = CHX[I]
#         Y1[I] = CHY[I]
#     end

#     # Compute X2, Y2 using dchxy
#     CFA[1] = 0.375
#     CFA[2] = -0.375
#     NST = 0
#     CHX, CHY, NTR = dchxy(TAU1, CFA, NST)
#     for I in 1:101
#         X2[I] = CHX[I]
#         Y2[I] = CHY[I]
#     end

#     # Compute AIL
#     AIL[1] = 0.01 / 3.0
#     CNU1 = 4.0 * AIL[1]
#     CNU2 = 2.0 * AIL[1]
#     for I in 2:2:100
#         AIL[I] = CNU1
#         AIL[I+1] = CNU2
#     end
#     AIL[101] = AIL[1]

#     # Initialize integrals
#     fill!(XA, 0.0)
#     fill!(XB, 0.0)

#     for I in 1:101
#         c1 = AIL[I] * X1[I] * AMU[I]
#         XA[1] += c1
#         XA[2] += c1 * AMU[I]
#         c2 = AIL[I] * Y1[I] * AMU[I]
#         XA[3] += c2
#         XA[4] += c2 * AMU[I]
#         c3 = AIL[I] * X2[I]
#         XB[1] += c3
#         XB[2] += c3 * AMU[I]
#         XB[3] += c3 * AMU[I]^2
#         XB[4] += c3 * AMU[I]^3
#         c4 = AIL[I] * Y2[I]
#         XB[5] += c4
#         XB[6] += c4 * AMU[I]
#         XB[7] += c4 * AMU[I]^2
#         XB[8] += c4 * AMU[I]^3
#     end

#     AI[1] = XB[1] + XB[5] - 8.0 / 3.0
#     AI[2] = XB[2] + XB[6]
#     AI[3] = XB[3] + XB[7]
#     AI[4] = XB[1] - XB[5] - 8.0 / 3.0
#     AI[5] = XB[2] - XB[6]
#     AI[6] = XB[3] - XB[7]
#     AI[7] = XB[4] - XB[8]
#     AI[8] = XA[1] + XA[3]
#     AI[9] = XA[2] + XA[4]
#     AI[10] = XA[1] - XA[3]
#     AI[11] = XA[2] - XA[4]

#     AI[12] = (AI[1] - AI[3]) / ((AI[4] - AI[6]) * TAU1 + 2.0 * (AI[5] - AI[7]))
#     AI[13] = 1.0 / (AI[4] * AI[10] - AI[5] * AI[11])
#     AI[14] = 1.0 / (AI[1] * AI[8] - AI[2] * AI[9] - 2.0 * AI[12] * (AI[5] * AI[8] - AI[4] * AI[9]))
#     AI[15] = 2.0 * (AI[8] * AI[10] - AI[9] * AI[11])
#     AI[16] = AI[13] * AI[15]
#     AI[17] = AI[14] * AI[15]

#     CNU1 = 0.5 * (AI[16] - AI[17])
#     CNU2 = 0.5 * (AI[16] + AI[17])

#     AI[15] = AI[13] * (AI[5] * AI[8] - AI[4] * AI[9])
#     AI[16] = AI[14] * (AI[2] * AI[10] - AI[1] * AI[11] - 2.0 * AI[12] * (AI[4] * AI[10] - AI[5] * AI[11]))
#     CNU3 = 0.5 * (AI[15] - AI[16])
#     CNU4 = 0.5 * (AI[15] + AI[16])

#     AI[15] = AI[13] * (AI[2] * AI[10] - AI[1] * AI[11])
#     AI[16] = AI[14] * (AI[5] * AI[8] - AI[4] * AI[9])
#     CU3 = 0.5 * (AI[15] - AI[16])
#     CU4 = 0.5 * (AI[15] + AI[16])

#     AI[15] = AI[14] * (AI[1] * AI[8] - AI[2] * AI[9])
#     SBAR = 1.0 - 0.375 * AI[12] * (AI[4] - AI[6]) *
#                  ((CNU2 - CNU1) * AI[8] + (CU4 - CU3) * AI[2] - AI[15] * AI[6])

#     AI[20] = 0.375 * AI[12] * (CNU2 - CNU1) * (AI[4] - AI[6])
#     AI[21] = 0.375 * AI[12] * (AI[4] - AI[6])
#     AI[22] = AI[21] * (CU4 - CU3)
#     AI[23] = AI[21] * AI[15]

#     for I in 1:101
#         GAML[I] = AI[20] * (X1[I] + Y1[I])
#         GAMR[I] = AI[22] * (X2[I] + Y2[I]) - AMU[I] * AI[23] * (X2[I] - Y2[I])
#     end
#     return GAMR, GAML, SBAR
# end

"""
    GAMMA(TAU1::Float64) -> (GAMR::Vector{Float64}, GAML::Vector{Float64}, SBAR::Float64)

Compute radiation scattered by a plane parallel homogeneous atmosphere
using the method of the X and Y functions.

# Arguments
- `TAU1::Float64`: Optical depth of the atmosphere.

# Returns
- `GAMR::Vector{Float64}`: Right-hand scattering function (length 101).
- `GAML::Vector{Float64}`: Left-hand scattering function (length 101).
- `SBAR::Float64`: Mean scattering function.
"""
# function GAMMA(TAU1::Float64)

#     # Preallocate
#     CHX  = zeros(Float64, 101)
#     CHY  = zeros(Float64, 101)
#     CFA  = zeros(Float64, 3)
#     AMU  = zeros(Float64, 101)
#     X1   = zeros(Float64, 101)
#     Y1   = zeros(Float64, 101)
#     X2   = zeros(Float64, 101)
#     Y2   = zeros(Float64, 101)
#     AIL  = zeros(Float64, 101)
#     AI   = zeros(Float64, 30)
#     GAMR = zeros(Float64, 101)
#     GAML = zeros(Float64, 101)

#     # Set up AMU array
#     AMU[1] = 0.0
#     for I in 2:101
#         AMU[I] = 0.01 * (I - 1)
#     end

#     # Compute X1, Y1 using dchxy
#     CFA[1] = 0.75
#     CFA[2] = -0.75
#     CFA[3] = 0.0
#     CHX, CHY, _ = dchxy(TAU1, CFA, 111)
#     X1 .= CHX
#     Y1 .= CHY
#     # Compute X2, Y2 using dchxy
#     CFA[1] = 0.375
#     CFA[2] = -0.375
#     CHX, CHY, _ = dchxy(TAU1, CFA, 0)
#     for I in 1:101
#         X2[I] = CHX[I]
#         Y2[I] = CHY[I]
#     end

#     # Compute AIL (quadrature weights)
#     AIL[1] = 0.01 / 3.0
#     CNU1   = 4.0 * AIL[1]
#     CNU2   = 2.0 * AIL[1]
#     for I in 2:2:100
#         AIL[I]   = CNU1
#         AIL[I+1] = CNU2
#     end
#     AIL[101] = AIL[1]

#     # Scalar accumulators instead of XA[1:4], XB[1:8]
#     xa1 = xa2 = xa3 = xa4 = 0.0
#     xb1 = xb2 = xb3 = xb4 = xb5 = xb6 = xb7 = xb8 = 0.0

#     for I in 1:101
#         a  = AMU[I]
#         a2 = a * a
#         a3 = a2 * a

#         c1 = AIL[I] * X1[I] * a
#         xa1 += c1
#         xa2 += c1 * a

#         c2 = AIL[I] * Y1[I] * a
#         xa3 += c2
#         xa4 += c2 * a

#         c3 = AIL[I] * X2[I]
#         xb1 += c3
#         xb2 += c3 * a
#         xb3 += c3 * a2
#         xb4 += c3 * a3

#         c4 = AIL[I] * Y2[I]
#         xb5 += c4
#         xb6 += c4 * a
#         xb7 += c4 * a2
#         xb8 += c4 * a3
#     end

#     # Fill AI vector
#     AI[1]  = xb1 + xb5 - 8.0 / 3.0
#     AI[2]  = xb2 + xb6
#     AI[3]  = xb3 + xb7
#     AI[4]  = xb1 - xb5 - 8.0 / 3.0
#     AI[5]  = xb2 - xb6
#     AI[6]  = xb3 - xb7
#     AI[7]  = xb4 - xb8
#     AI[8]  = xa1 + xa3
#     AI[9]  = xa2 + xa4
#     AI[10] = xa1 - xa3
#     AI[11] = xa2 - xa4

#     AI[12] = (AI[1] - AI[3]) / ((AI[4] - AI[6]) * TAU1 + 2.0 * (AI[5] - AI[7]))
#     AI[13] = 1.0 / (AI[4] * AI[10] - AI[5] * AI[11])
#     AI[14] = 1.0 / (AI[1] * AI[8] - AI[2] * AI[9] -
#                     2.0 * AI[12] * (AI[5] * AI[8] - AI[4] * AI[9]))
#     AI[15] = 2.0 * (AI[8] * AI[10] - AI[9] * AI[11])
#     AI[16] = AI[13] * AI[15]
#     AI[17] = AI[14] * AI[15]

#     CNU1 = 0.5 * (AI[16] - AI[17])
#     CNU2 = 0.5 * (AI[16] + AI[17])

#     AI[15] = AI[13] * (AI[5] * AI[8] - AI[4] * AI[9])
#     AI[16] = AI[14] * (AI[2] * AI[10] - AI[1] * AI[11] -
#                        2.0 * AI[12] * (AI[4] * AI[10] - AI[5] * AI[11]))
#     #CNU3 = 0.5 * (AI[15] - AI[16])
#     #CNU4 = 0.5 * (AI[15] + AI[16])

#     AI[15] = AI[13] * (AI[2] * AI[10] - AI[1] * AI[11])
#     AI[16] = AI[14] * (AI[5] * AI[8] - AI[4] * AI[9])
#     CU3   = 0.5 * (AI[15] - AI[16])
#     CU4   = 0.5 * (AI[15] + AI[16])

#     AI[15] = AI[14] * (AI[1] * AI[8] - AI[2] * AI[9])
#     SBAR  = 1.0 - 0.375 * AI[12] * (AI[4] - AI[6]) *
#                    ((CNU2 - CNU1) * AI[8] + (CU4 - CU3) * AI[2] - AI[15] * AI[6])

#     AI[20] = 0.375 * AI[12] * (CNU2 - CNU1) * (AI[4] - AI[6])
#     AI[21] = 0.375 * AI[12] * (AI[4] - AI[6])
#     AI[22] = AI[21] * (CU4 - CU3)
#     AI[23] = AI[21] * AI[15]

#     for I in 1:101
#         GAML[I] = AI[20] * (X1[I] + Y1[I])
#         GAMR[I] = AI[22] * (X2[I] + Y2[I]) - AMU[I] * AI[23] * (X2[I] - Y2[I])
#     end

#     return GAMR, GAML, SBAR
# end
using StaticArrays

function allocate_gamma()
    (;
        CHX  = zeros(101),
        CHY  = zeros(101),
        AMU  = zeros(101),
        X1   = zeros(101),
        Y1   = zeros(101),
        X2   = zeros(101),
        Y2   = zeros(101),
        AIL  = zeros(101),
        GAMR = zeros(101),
        GAML = zeros(101),
        dchxy_buffers = init_dchxy_buffers(),
    )
end

GAMMA(TAU1::Float64) = GAMMA!(allocate_gamma(), TAU1)

function GAMMA!(buffers, TAU1::Float64)
    # Large arrays (mutable, normal)
    (; CHX, CHY, AMU, X1, Y1, X2, Y2, AIL, GAMR, GAML, dchxy_buffers) = buffers

    # Small fixed-size arrays (use StaticArrays)
    AI  = @MVector zeros(30)

    # Set up AMU array
    AMU[1] = 0.0
    for I in 2:101
        AMU[I] = 0.01 * (I - 1)
    end


    # Compute X1, Y1 using dchxy
    CFA = SVector((0.75, -0.75, 0.0))
    # CHX_, CHY_, _ = dchxy!(dchxy_buffers, TAU1, collect(CFA), 111)
    CHX_, CHY_, _ = dchxy!(dchxy_buffers, TAU1, CFA, 111)
    X1 .= CHX_
    Y1 .= CHY_

    # Compute X2, Y2 using dchxy
    CFA = SVector((0.375, -0.375, 0.0))
    # CHX_, CHY_, _ = dchxy!(dchxy_buffers, TAU1, collect(CFA), 0)
    CHX_, CHY_, _ = dchxy!(dchxy_buffers, TAU1, CFA, 0)
    X2 .= CHX_
    Y2 .= CHY_

    # Compute AIL (quadrature weights)
    AIL[1] = 0.01 / 3.0
    CNU1 = 4.0 * AIL[1]
    CNU2 = 2.0 * AIL[1]
    for I in 2:2:100
        AIL[I] = CNU1
        AIL[I+1] = CNU2
    end
    AIL[101] = AIL[1]

    # Scalar accumulators
    xa1 = xa2 = xa3 = xa4 = 0.0
    xb1 = xb2 = xb3 = xb4 = xb5 = xb6 = xb7 = xb8 = 0.0

    for I in 1:101
        a  = AMU[I]
        a2 = a * a
        a3 = a2 * a

        c1 = AIL[I] * X1[I] * a
        xa1 += c1
        xa2 += c1 * a

        c2 = AIL[I] * Y1[I] * a
        xa3 += c2
        xa4 += c2 * a

        c3 = AIL[I] * X2[I]
        xb1 += c3
        xb2 += c3 * a
        xb3 += c3 * a2
        xb4 += c3 * a3

        c4 = AIL[I] * Y2[I]
        xb5 += c4
        xb6 += c4 * a
        xb7 += c4 * a2
        xb8 += c4 * a3
    end

    # Fill AI vector
    AI[1]  = xb1 + xb5 - 8.0 / 3.0
    AI[2]  = xb2 + xb6
    AI[3]  = xb3 + xb7
    AI[4]  = xb1 - xb5 - 8.0 / 3.0
    AI[5]  = xb2 - xb6
    AI[6]  = xb3 - xb7
    AI[7]  = xb4 - xb8
    AI[8]  = xa1 + xa3
    AI[9]  = xa2 + xa4
    AI[10] = xa1 - xa3
    AI[11] = xa2 - xa4

    AI[12] = (AI[1] - AI[3]) / ((AI[4] - AI[6]) * TAU1 + 2.0 * (AI[5] - AI[7]))
    AI[13] = 1.0 / (AI[4] * AI[10] - AI[5] * AI[11])
    AI[14] = 1.0 / (AI[1] * AI[8] - AI[2] * AI[9] -
                    2.0 * AI[12] * (AI[5] * AI[8] - AI[4] * AI[9]))
    AI[15] = 2.0 * (AI[8] * AI[10] - AI[9] * AI[11])
    AI[16] = AI[13] * AI[15]
    AI[17] = AI[14] * AI[15]

    CNU1 = 0.5 * (AI[16] - AI[17])
    CNU2 = 0.5 * (AI[16] + AI[17])

    AI[15] = AI[13] * (AI[5] * AI[8] - AI[4] * AI[9])
    AI[16] = AI[14] * (AI[2] * AI[10] - AI[1] * AI[11] -
                       2.0 * AI[12] * (AI[4] * AI[10] - AI[5] * AI[11]))

    AI[15] = AI[13] * (AI[2] * AI[10] - AI[1] * AI[11])
    AI[16] = AI[14] * (AI[5] * AI[8] - AI[4] * AI[9])
    CU3   = 0.5 * (AI[15] - AI[16])
    CU4   = 0.5 * (AI[15] + AI[16])

    AI[15] = AI[14] * (AI[1] * AI[8] - AI[2] * AI[9])
    SBAR  = 1.0 - 0.375 * AI[12] * (AI[4] - AI[6]) *
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

function init_dchxy_buffers()
    arrays = (;
        PSI = zeros(101),
        AMU = zeros(101),
        XA = zeros(101),
        XB = zeros(101),
        FNPP = zeros(101),
        FNPN = zeros(101),
        FNC0 = zeros(101),
        FNC1 = zeros(101),
        FNX = zeros(101),
        FNY = zeros(101),
        FNW = zeros(101),
        FMC0 = zeros(101),
        FMC1 = zeros(101),
        XD = zeros(101),
        XE = zeros(101),
        CHXA = zeros(101),
        CHYA = zeros(101),
        CHX = zeros(101),
        CHY = zeros(101),
    )
    return arrays
end

"""
    dchxy(TAU1::Float64, CFA::Vector{Float64}, NCASE::Int) 
        -> (CHX::Vector{Float64}, CHY::Vector{Float64}, nomitr::Int)

Compute Chandrasekhar's X and Y functions for radiative transfer.

# Description
This routine evaluates the X- and Y-functions of Chandrasekhar using
double precision arithmetic. The method starts with the fourth
approximation given in Sec. 59 of Chandrasekhar’s *Radiative Transfer*
(Dover Publications, 1960), and iteratively refines the values
according to the procedure in Sec. 60. Iteration terminates when
successive corrected values of the Y-function agree to four significant
figures.

# Inputs
- `TAU1::Float64`:  
  Normal optical thickness of the atmosphere.  
  Must be ≤ 2.0.

- `CFA::NTuple{3,Float64}`:  
  Coefficients of the characteristic function in polynomial form:  
  ```math
  C(μ) = Σⱼ Aⱼ * μ^(2(j-1)),   j = 1,2,3

Outputs

CHX::Vector{Float64}
Values of the X-function at 101 evenly spaced μ values from 0.00 to 1.00 in steps of 0.01.

CHY::Vector{Float64}
Values of the Y-function at the same μ grid.

nomitr::Int
Number of iterations performed before convergence.

Notes

If ncase != 0, a conservative case is assumed and a standard solution is returned.
The program terminates with an error if:
- tau1 > 2.0
- the characteristic function is negative for any μ
- the integral of the characteristic function exceeds 0.5

References

https://en.wikipedia.org/wiki/Chandrasekhar%27s_X-_and_Y-function

McCullough, E. C., & Porter, W. P. (1971). Computing clear day solar radiation 
spectra for the terrestrial ecological environment. Ecology, 52(6), 1008–1015.
     https://doi.org/10.2307/1933806

"""
dchxy(TAU1::Float64, CFA::Vector{Float64}, NCASE::Int) =
    dchxy!(init_dchxy_buffers(), TAU1, CFA, NCASE)

function dchxy!(buffers, TAU1::Float64, CFA::AbstractVector{Float64}, NCASE::Int)
    PSI   = buffers.PSI
    AMU   = buffers.AMU
    XA    = buffers.XA
    XB    = buffers.XB
    UMA  = @MVector zeros(5)
    ACAP = @MVector zeros(5)
    TEMX = @MVector zeros(8)
    TEMY = @MVector zeros(8)
    RTK  = @MVector zeros(5)
    ALAM = @MVector zeros(5)
    FNPP  = buffers.FNPP
    FNPN  = buffers.FNPN
    FNC0  = buffers.FNC0
    FNC1  = buffers.FNC1
    FNX   = buffers.FNX
    FNY   = buffers.FNY
    FNW   = buffers.FNW
    FMC0  = buffers.FMC0
    FMC1  = buffers.FMC1
    XD    = buffers.XD    # equivalence
    XE    = buffers.XE    # equivalence
    CHXA  = buffers.CHXA  # equivalence
    CHYA  = buffers.CHYA  # equivalence
    CHX   = buffers.CHX
    CHY   = buffers.CHY
    XA    = buffers.XA
    XB    = buffers.XB

    # Variables
    PERA = 0.0

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
        #Printf.printf("%12.5E %12.5E %12.5E\n", CFA[1], CFA[2], CFA[3])
        #Printf.printf("%12.5E\n", TAU1)
        #Printf.printf("\n")
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
    XB[1] = TAU1 == 0.0 ? 1.0 : 0.0 # Fortran line 345

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

    converged, nomitr = _dchxy_converge!(FNX, FNY, AMU, PSI, XA, XB, XD, XE, CHX, CHY, CHXA, CHYA)

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



# """
#     dexpi(x::Float64) -> Float64

# Compute the exponential integral with ~15 significant figure accuracy.  

# Implements the algorithm described in the original FORTRAN code, which switches 
# between polynomial ratios and numerical quadrature depending on the range of `x`.

# # Inputs
# - `x::Float64`  
#   Argument of the exponential integral.

# # Output
# - `E1::Float64`  
#   The exponential integral evaluated at `x`.

# # Method
# Different computational strategies are applied depending on the sign and magnitude of `x`:

# - **For negative `x`:**
#   - `x > -1.0e-20` → `γ + log(|x|)`  
#   - `-1.0e-20 ≥ x > -1.5` → 3-point Gaussian quadrature  
#   - `-1.5 ≥ x > -4.65` → ratio of two 7-term polynomials  
#   - `-4.65 ≥ x > -12.0` → ratio of two 6-term polynomials  
#   - `-12.0 ≥ x > -170.0` → 12-point Gauss–Laguerre quadrature  

# - **For positive `x`:**
#   - `x < 1.0e-20` → `γ + log(x)`  
#   - `1.0e-20 ≤ x ≤ 40.0` → 12-point Gaussian quadrature  
#   - `40.0 < x ≤ 173.0` → 12-point Gauss–Laguerre quadrature  

# # Notes
# - `γ` denotes the Euler–Mascheroni constant.  
# - Accuracy is approximately 15 significant figures across the supported domain.  
# - Outside the ranges listed above, behavior is not guaranteed.
# """
# function dexpi(x::Float64)
#     # Constants
#     gamma = 0.57721566490153286

#     # Coefficients
#     A1 = [0.1193095930415985, 0.3306046932331323, 0.4662347571015760]
#     B1 = [0.4679139345726910, 0.3607615730481386, 0.1713244923791703]
#     A2 = [0.02823912701457739, 30.52042817823067, 215.8885931211323,
#         410.4611319636983, 278.5527592726121, 71.33086969436196, 0.5758931590224375]
#     B2 = [10.0, 138.3869728490638, 488.08581830736, 634.8804630786363,
#         344.1289899236299, 77.08964199043784, 0.5758934565014882]
#     A3 = [0.07630772325814641, 21.23699219410890, 47.45350785776186,
#         29.66421696379266, 6.444800036068992, 0.04295808082119383]
#     B3 = [10.0, 52.78950049492932, 71.96111390658517, 35.67945294128107,
#         6.874380519301884, 0.04295808112146861]
#     A4 = [0.1157221173580207, 0.6117574845151307, 1.512610269776419,
#         2.833751337743507, 4.599227639418348, 6.844525453115177,
#         9.621316842456867, 13.00605499330635, 17.11685518746226,
#         22.15109037939701, 28.48796725098400, 37.09912104446692]
#     B4 = [0.2647313710554432, 0.3777592758731380, 0.2440820113198776,
#         0.09044922221168093, 0.02010238115463410, 0.002663973541865316,
#         0.0002032315926629994, 8.365055856819799e-5, 1.668493876540910e-6,
#         1.342391030515004e-8, 3.061601635035021e-11, 8.148077467426242e-15]
#     A5 = [0.03202844643130281, 0.09555943373680816, 0.1575213398480817,
#         0.2168967538130226, 0.2727107356944198, 0.3240468259684878,
#         0.3700620957892772, 0.4100009929869515, 0.4432077635022005,
#         0.4691372760013664, 0.4873642779856547, 0.4975936099985107]
#     B5 = [0.1279381953467522, 0.1258374563468283, 0.1216704729278034,
#         0.1155056680537256, 0.1074442701159656, 0.09761865210411389,
#         0.08619016153195328, 0.07334648141108031, 0.05929858491543678,
#         0.04427743881741981, 0.02853138862893366, 0.01234122979998720]

#     if x == 0.0
#         error("The argument of DEXPI is very close to zero.")
#     elseif x < 0.0
#         ax = abs(x)
#         if x > -1e-20
#             return log(ax) + gamma
#         elseif x > -1.5
#             yy = exp(-0.5 * ax)
#             s = 0.0
#             for i in 1:3
#                 yz = exp(A1[i] * ax)
#                 s += B1[i] * ((1 - yy / yz) / (0.5 + A1[i]) + (1 - yy * yz) / (0.5 - A1[i]))
#             end
#             return -0.5 * s + log(ax) + gamma
#         elseif x > -4.65
#             sumn = evalpoly(ax, reverse(A2))
#             sumd = evalpoly(ax, reverse(B2))
#             return (sumn / (sumd * x)) * exp(x)
#         elseif x > -12.0
#             sumn = evalpoly(ax, reverse(A3))
#             sumd = evalpoly(ax, reverse(B3))
#             return (sumn / (sumd * x)) * exp(x)
#         elseif x > -170.0
#             dexpi = 0.0
#             for j in 1:12
#                 dexpi += B4[j] / (1 + A4[j] / ax)
#             end
#             return (exp(x) / ax) * (-dexpi)
#         else
#             return 0.0
#         end
#     else
#         if x <= 1e-20
#             return log(x) + gamma
#         elseif x <= 40.0
#             yy = exp(0.5 * x)
#             dexpi = 0.0
#             for j in 1:12
#                 yz = exp(-A5[j] * x)
#                 dexpi += ((1 - yy / yz) / (0.5 + A5[j]) + (1 - yy * yz) / (0.5 - A5[j])) * B5[j]
#             end
#             return -0.5 * dexpi + log(x) + gamma
#         elseif x <= 173.0
#             dexpi = 0.0
#             for j in 1:12
#                 dexpi += B4[j] / (1 - A4[j] / x)
#             end
#             return (exp(x) / x) * dexpi
#         else
#             error("The argument of DEXPI is very large.")
#         end
#     end
# end

"""
    solrad(solar_radiaion_model; kw...)

Compute clear sky solar radiation at a given place and time using a detailed atmospheric radiative transfer model.

# Arguments

- `solar_radiaion_model`:

# Keyword Arguments

- `days::Vector{Float64}`: Days of the year (1–365/366) to evaluate.
- `hours::Vector{Float64}`: Decimal hours of the day (0.0–23.0).
- `latitude::Quantity`: Latitude in degrees, e.g. `43.0u"°"`.
- `longitude_correction::Real=0.0`: Longitude correction in hours (positive west of standard meridian).
- `year::Real`: Year used for ozone table lookup.
- `terrain`
- `albedo::Vector{<:Real}=fill(0.15, length(days))`: Daily ground albedo, fraction [0, 1].

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
- `ZenithSlope::Float64`: Slope solar zenith angle [degrees].
- `Azimuth::Float64`: Solar azimuth angle [degrees].
- `doy::Int`: Day of year.
- `hour::Float64`: Decimal hour of day.

# Notes
- Radiation units are returned in `W/m²`. Internally, units like `mW/cm²` are used and converted as necessary using `Unitful.jl`.
- Topographic shading is included via the `horizon_angles` input (horizon angle mask) but cloud effects on scattered solar should be added later.
- Outputs are computed for each (day, hour) combination in the input vectors.
- In optical air mass 'arims' calculation the difference between apparent and true zenith angle is neglected for z less than 88 degrees
- Variation of airms with altitude is ignored since it is negligible up to at least 6 km above sea level

# References

FD and FDQ derived from tables in Dave and Furukawa (1967)
Dave, J. V., & Furukawa, P. M. (1966). Scattered radiation in the ozone 
 absorption bands at selected levels of a terrestrial, Rayleigh atmosphere (Vol. 7).
 Americal Meteorological Society.
"""
#solrad(::McCulloughPorterSolarGeometry, args::Vararg{Missing}; kwargs...) = missing
function solrad(model::SolarRadiation, args...; kwargs...)
    any(ismissing, args) && return missing
    return solrad_core(model, args...; kwargs...)
end
function solrad(solar_model::SolarRadiation, latitude::Quantity, elevation::Quantity, 
    slope::Quantity, aspect::Quantity, P_atmos::Quantity, albedo::Number;
    days::Real,#Vector{<:Real}=[15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349],
    year::Real=2001, # TODO: this shouldn't have a default
    #terrain::Terrain,
    longitude_correction::Real=0.0, # longitude correction, hours
    hours::Number,#AbstractVector{<:Real}=0:1:23,
    horizon_angles::AbstractVector{<:Quantity}, #albedo::Vector{<:Real}, # substrate albedo (decimal %)
)
    (; solar_geometry_model, cmH2O, iuv, scattered, amr, nmax, Iλ, OZ, τR, τO, τA, τW, Sλ, FD, FDQ, s̄) = solar_model
    #(; elevation, horizon_angles, slope, P_atmos) = terrain

    ndays = length(days)    # number of days
    ntimes = length(hours)  # number of times
    nsteps = ndays * ntimes # total time steps

    # arrays to hold every time step's radiation between 300 and 320 nm in 2 nm steps
    GRλs = fill(0.0u"mW/nm/cm^2", nsteps, nmax) # wavelength-specific global radiation
    DRRλs = fill(0.0u"mW/nm/cm^2", nsteps, nmax)# wavelength-specific direct Rayleigh radiation
    DRλs = fill(0.0u"mW/nm/cm^2", nsteps, nmax) # wavelength-specific direct radiation
    SRλs = fill(0.0u"mW/nm/cm^2", nsteps, nmax) # wavelength-specific scattered radiation
    GRs = fill(0.0u"mW/cm^2", nsteps)           # total global radiation
    DRRs = fill(0.0u"mW/cm^2", nsteps)          # total direct Rayleigh radiation
    DRs = fill(0.0u"mW/cm^2", nsteps)           # total direct radiation
    SRs = fill(0.0u"mW/cm^2", nsteps)           # total scattered radiation
    gamma_buffers = allocate_gamma()

    # arrays to hold zenith and azimuth angles each step
    Zs = fill(90.0u"°", nsteps)                 # zenith angles
    ZSLs = fill(90.0u"°", nsteps)               # slope zenith angles
    AZIs = Vector{Union{Missing,typeof(0.0u"°")}}(undef, nsteps)
    fill!(AZIs, 90.0u"°")   
    HHs = fill(0.0, ndays)                      # hour angles
    tsns = fill(0.0, ndays)                     # hour at solar noon
    DOYs = Vector{Int}(undef, nsteps)           # day of year
    times = Vector{Real}(undef, nsteps)         # time
    step = 1
    HH = 0.0 # initialise sunrise hour angle
    tsn = 12.0 # initialise time of solar noon
    for i in 1:ndays
        # arrays to hold radiation for a given hour between 300 and 320 nm in 2 nm steps
        GRINT = fill(0.0u"mW/cm^2", nmax)   # integrated global radiation component (direct + scattered)
        DRRINT = fill(0.0u"mW/cm^2", nmax)  # integrated direct Rayleigh radiation component
        DRINT = fill(0.0u"mW/cm^2", nmax)   # integrated direct radiation component
        SRINT = fill(0.0u"mW/cm^2", nmax)   # integrated scattered radiation component
        AIλ = fill(0.0u"nm", nmax)
        GRλ = GRINT * u"1/nm"               # wavelength-specific global radiation component (direct + scattered)
        DRRλ = GRINT * u"1/nm"              # wavelength-specific direct Rayleigh radiation component
        DRλ = GRINT * u"1/nm"               # wavelength-specific direct radiation component
        SRλ = GRINT * u"1/nm"               # wavelength-specific scattered radiation component
        alb = albedo#[i]
        for j in 1:ntimes
            d = days[i]
            t = hours[j]
            h, tsn = hour_angle(t, longitude_correction) # hour angle (radians)
            (; ζ, δ, z, AR2) = solar_geometry(solar_geometry_model, latitude; d, h) # compute ecliptic, declination, zenith angle and (a/r)^2
            Z = uconvert(u"°", z)
            Zsl = Z
            amult = 1.0
            if sign(latitude) < 0
                amult = -1.0
            end
            if Z < 107.0u"°"
                if Z > 88.0u"°"
                    # Compute skylight based on G.V. Rozenberg. 1966. Twilight. Plenum Press.
                    # p. 18,19.  First computing lumens: y = b - mx. Regression of data is:
                    Elog = 41.34615384 - 0.423076923 * ustrip(u"°", Z)
                    # Converting lux (lumen/m2) to W/m2 on horizontal surface -
                    # Twilight - scattered skylight before sunrise or after sunset
                    # From p. 239 Documenta Geigy Scientific Tables. 1966. 6th ed. K. Diem, ed.
                    # Mech./elect equiv. of light = 1.46*10^-3 kW/lumen
                    Skylum = (10.0^Elog) * 1.46E-03u"mW * cm^-2"
                    SRINT[nmax] = Skylum
                    GRINT[nmax] = SRINT[nmax]
                    GRs[step] = GRINT[nmax]
                    SRs[step] = SRINT[nmax]
                end
            end
            # testing cos(h) to see if it exceeds +1 or -1
            TDTL = -tan(δ) * tan(latitude) # from eq.7 McCullough & Porter 1971
            if abs(TDTL) >= 1 # long day or night
                H = π
            else
                H = abs(acos(TDTL))
            end
            # check if sunrise
            HH = 12.0 * H / π
            ts = t - tsn

            sun_up = true

            if ts <= 0.0 && abs(ts) > HH
                sun_up = false
            elseif ts > 0.0 && ts >= HH
                sun_up = false
            end

            if sun_up || TDTL == 1 # sun is up, proceed
                alt = (π / 2 - z)u"rad"
                altdeg = uconvert(u"°", alt).val
                # tazsun corresponds to tangent of azimuth
                tazsun = sin(h) / (cos(latitude) * tan(δ) - sin(latitude) * cos(h))
                # sun azimuth in radians
                azsun = atan(tazsun) * amult
                # azimuth in degrees
                dazsun = uconvert(u"°", azsun)
                # correcting for hemisphere/quadrant
                if h <= 0.0
                    # Morning - east of reference
                    if dazsun <= 0.0u"°"
                        # 1st Quadrant (0–90°)
                        dazsun = -1.0 * dazsun
                    else
                        # 2nd Quadrant (90–180°)
                        dazsun = 180.0u"°" - dazsun
                    end
                else
                    # Afternoon - west of reference
                    if dazsun < 0.0u"°"
                        # 3rd Quadrant (180–270°)
                        dazsun = 180.0u"°" - dazsun
                    else
                        # 4th Quadrant (270–360°)
                        dazsun = 360.0u"°" - dazsun
                    end
                end
                # Special case: hour angle = 0
                if h == 0.0
                    dazsun = 180.0u"°"
                end

                cz = cos(z)
                intcz = Int(floor(100.0 * cz + 1.0))
                Z = uconvert(u"°", z)  # zenith angle in degrees

                # horizon angle - check this works when starting at 0 rather than e.g. 15 deg
                azi = range(0u"°", stop=360u"°" - 360u"°" / length(horizon_angles), length=length(horizon_angles))
                ahoriz = horizon_angles[argmin(abs.(dazsun .- azi))]

                # slope zenith angle calculation (Eq. 3.15 in Sellers 1965. Physical Climatology. U. Chicago Press)
                if slope > 0u"°"
                    czsl = cos(z) * cos(slope) + sin(z) * sin(slope) * cos(dazsun - aspect)
                    zsl = acos(czsl)
                    Zsl = min(uconvert(u"°", zsl), 90u"°") # cap at 90 degrees if sun is below slope horizon
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
                cz = cos(z)
                intcz = Int(floor(100.0 * cz + 1.0))
                Z = uconvert(u"°", z)  # zenith angle in degrees

                # atmospheric ozone lookup
                # convert latitude in degrees to nearest 10-degree index
                tlat = (latitude + 100.0u"°") / 10.0u"°"
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

                (; molecular_corr, aerosol_corr, ozone_corr, water_vapour_corr) = elevation_correction(elevation)

                P = P_atmos

                for N in 1:nmax
                    τλ1 = (P / 101300u"Pa") * τR[N] * molecular_corr
                    τλ2 = (25.0u"km" / amr) * τA[N] * aerosol_corr
                    τλ3 = (ozone / 0.34) * τO[N] * ozone_corr
                    τλ4 = τW[N] * sqrt(airms * cmH2O * water_vapour_corr)
                    τλ = ((float(τλ1) + τλ2 + τλ3) * airms) + τλ4

                    if τλ > 80.0 # making sure that at low sun angles air mass doesn't make τλ too large
                        τλ = 80.0
                    end

                    part1 = Sλ[N] * AR2 * cz
                    part2 = τλ > 0.0 ? exp(-τλ) : 0.0
                    if part2 < 1.0e-24
                        DRλ[N] = 0.0u"mW / cm^2 / nm"
                    else
                        # TODO: ustrip to what
                        DRλ[N] = ((ustrip(part1) * part2) / 1000.0) * u"mW / cm^2 / nm"
                    end

                    # so the integrator doesn't get confused at very low sun angles
                    if DRλ[N] < 1.0e-25u"mW / cm^2 / nm"
                        DRλ[N] = 1.0e-25u"mW / cm^2 / nm"
                    end

                    DRRλ[N] = (Sλ[N] * AR2 * cz) * exp(-float(τλ1) * airms) / 1000.0

                    if altdeg < ahoriz
                        DRλ[N] = 1.0e-25u"mW / cm^2 / nm"
                        DRRλ[N] = 1.0e-25u"mW / cm^2 / nm"
                    end

                    # Sky (SRλ) and Global Radiation (GRλ)
                    if scattered == false
                        SRλ[N] = 0.0u"mW / cm^2 / nm"
                    elseif iuv
                        if τλ1 >= 0.03
                            GAMR, GAML, SBAR = GAMMA!(gamma_buffers, τλ1)
                            SRλ[N] = (
                                         ((float(GAML[intcz]) + float(GAMR[intcz])) / (2.0 * (1.0 - alb * float(SBAR))))
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
                            # The option iuv = false has caused the program to enter this section which
                            # computes scattered radiation (SRλ) for 290 nm to 360 nm using a theory
                            # of radiation scattered from a Rayleigh (molecular) atmosphere with
                            # ozone absorption. The functions needed for the computation are stored
                            # as FD(N,I) and FDQ(N,I) where N is the wavelength index and I is
                            # (zenith angle + 5)/5 rounded off to the nearest integer value.
                            # The arrays FD and FDQ are for sea level (P = 1013 mb).
                            # TODO: ustrip to what
                            B = ustrip(Z) / 5
                            IA = trunc(Int, B)
                            C = B - IA
                            if C > 0.5
                                I = IA + 2
                            else
                                I = IA + 1
                            end
                            FDAV = FD[N, I]
                            FDQDAV = FDQ[N, I]
                            SRλ[N] = (Sλ[N] / π) * (FDAV + FDQDAV * (alb / (1.0 - (alb * s̄[N])))) / 1000.0
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
                GRλs[step, :] .= GRλ
                DRRλs[step, :] .= DRRλ
                DRλs[step, :] .= DRλ
                SRλs[step, :] .= SRλ
                GRs[step] = GRINT[nmax]
                DRRs[step] = DRRINT[nmax]
                DRs[step] = DRINT[nmax]
                SRs[step] = SRINT[nmax]
            else # sunrise, sunset or long day
                dazsun = missing
            end
            # Store into row `step`
            Zs[step] = Z
            ZSLs[step] = Zsl
            AZIs[step] = dazsun
            DOYs[step] = d
            times[step] = t
            step += 1
        end
        HHs[i] = HH     # save today's sunrise hour angle
        tsns[i] = tsn   # save today's time of sunrise
    end

    return (
        zenith_angle = Zs,
        zenith_slope_angle = ZSLs,
        azimuth_angle = AZIs,
        hour_angle_sunrise = HHs,
        hour_solar_noon = tsns,
        day_of_year = DOYs,
        hour = times,
        # TODO remove all this allocation from broadcasts
        # why is this conversion needed, what is the 10 about
        rayleigh_total = DRRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        direct_total = DRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        diffuse_total = SRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        global_total = GRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        wavelength = Iλ,
        rayleigh_spectra = DRRλs .* (10u"W/m^2" / 1u"mW/cm^2"),
        direct_spectra = DRλs .* (10u"W/m^2" / 1u"mW/cm^2"),
        diffuse_spectra = SRλs .* (10u"W/m^2" / 1u"mW/cm^2"),
        global_spectra = GRλs .* (10u"W/m^2" / 1u"mW/cm^2"),
    )
end

abstract type AbstractAtmosphericRadiationModel end
struct SwinbankAtmosphericRadiation <: AbstractAtmosphericRadiationModel end
struct CampbellNormanAtmosphericRadiation <: AbstractAtmosphericRadiationModel end

function atmospheric_radiation(::SwinbankAtmosphericRadiation, P_vap, tair)
    # Swinbank, Eq. 10.11 in Campbell and Norman 1998
    arad = uconvert(u"W*m^-2", ((9.2e-6 * (u"K"(tair))^2) * σ * (u"K"(tair))^4) / 1u"K^2")
    return P_vap, arad
end
function atmospheric_radiation(::CampbellNormanAtmosphericRadiation, P_vap, tair)
    # Campbell and Norman 1998 eq. 10.10 to get emissivity of sky
    arad = u"W/m^2"((1.72 * (ustrip(u"kPa", P_vap) / ustrip(u"K", tair + 0.01u"K"))^(1//7)) * σ * (u"K"(tair) + 0.01u"K")^4) 
    return P_vap, arad
end

function longwave_radiation(radiation_model=CampbellNormanAtmosphericRadiation(); 
    terrain, 
    environment_instant,
    surface_temperature,
)
    # TODO these are not the real names
    (; elevation, P_atmos, viewfactor) = terrain
    (; reference_humidity, reference_temperature, surface_emissivity, cloud_emissivity, cloud_cover, shade) = environment_instant

    # Short names, hardly worth it
    tsurf = surface_temperature
    tair = reference_temperature
    rh = reference_humidity
    slep = surface_emissivity
    sle = cloud_emissivity
    cloud = cloud_cover

    # Longwave radiation (handle both IR modes)
    wet_air_out = wet_air_properties(u"K"(tair); rh, P_atmos)

    # Atmospheric radiation
    P_vap, arad = atmospheric_radiation(radiation_model, wet_air_out.P_vap, tair)

    # Cloud radiation temperature (shade approximation, TAIR - 2°C)
    crad = σ * slep * (u"K"(tair) - 2.0u"K")^4

    # Hillshade radiation temperature (approximated as air temperature)
    hrad = σ * slep * (u"K"(tair))^4

    # Ground surface radiation temperature
    srad = σ * sle * (u"K"(tsurf))^4

    # Clear sky fraction
    clr = 1.0 - cloud / 100.0
    clear = arad * clr
    clod = crad * (cloud / 100.0)
    qradsk = (clear + clod) * ((100 - shade) / 100.0)
    qradvg = (shade / 100.0) * hrad
    qradgr = ((100.0 - shade) / 100.0) * srad + (shade / 100.0) * hrad
    qradhl = hrad
    qrad = (qradsk + qradvg) * viewfactor + qradhl * (1.0 - viewfactor) - qradgr
    tsky = (((qradsk + qradvg) * viewfactor + qradhl * (1.0 - viewfactor)) / σ)^(1//4)

    return (;
        # TODO standardise these names with their target uses
        Tsky=u"K"(tsky),
        Qrad=qrad, # e.g. this is Q_infrared in `solar_radiation`
        Qrad_sky=qradsk,
        Qrad_veg=qradvg,
        Qrad_ground=qradgr,
        Qrad_hill=qradhl
    )
end

"""
    cloud_adjust_radiation(cloud, D_cs, B_cs, zenith, doy; a=0.25, b=0.5, gamma=1.0)

Compute global (G), diffuse (D), and direct-beam (B) solar on a horizontal surface
given cloud cover fraction `cloud` (0–1), clear-sky diffuse `D_cs` and direct `B_cs`,
solar zenith angle `zenith` (radians), and day-of-year `doy`.

- Ångström scaling: G = (a + b*S) * (D_cs + B_cs), with S ≈ (1 - cloud)^gamma
- Diffuse fraction via Erbs (uses extraterrestrial horizontal irradiance) via
    a clearness index (Maxwwell 1987) which is the ratio of global to extraterrestrial
    irradiance on a horizontal plane

Returns `(G, D, B)`; works with arrays but needs to not use 'similar' if to work with
    scalars.

Reference
Maxwell, E. L., "A Quasi-Physical Model for Converting Hourly
           Global Horizontal to Direct Normal Insolation", Technical
           Report No. SERI/TR-215-3087, Golden, CO: Solar Energy Research
           Institute, 1987.
"""
function cloud_adjust_radiation!(output, cloud::AbstractArray, D_cs, B_cs, zenith::AbstractArray, doy; 
    a=0.36, b=0.64, gamma=1.0,
)
    (; global_solar, diffuse_solar, direct_solar) = output
    G, D, B = (global_solar, diffuse_solar, direct_solar)
    # Solar geometry
    cosz     = cos.(zenith)
    cosz_pos = max.(cosz, 0.0)

    # 1) Extraterrestrial horizontal irradiance (W/m²)
    I_sc  = 1367.0u"W/m^2"
    E0    = 1.00011 .+ 0.034221*cosd.(360.0 .* (doy .- 1) ./ 365.0) .+
                     0.00128*sind.(360.0 .* (doy .- 1) ./ 365.0) .+
                     0.000719*cosd.(2 .* 360.0 .* (doy .- 1) ./ 365.0) .+
                     0.000077*sind.(2 .* 360.0 .* (doy .- 1) ./ 365.0)
    G0h   = I_sc .* E0 #.* cosz_pos  # on horizontal; 0 at night

    # 2) Ångström–Prescott scaling of clear-sky global by cloud cover
    S     = (1 .- cloud).^gamma                 # approx. sunshine fraction
    T     = a .+ b .* S                         # transmittance
    G_cs  = D_cs .+ B_cs                        # clear-sky global
    G     .= max.(T .* G_cs, 0.0u"W/m^2") #.* (cosz_pos .> 0)  # zero at night

    # 3) Split G into diffuse/direct using Erbs diffuse fraction vs clearness index K_t
    ϵ     = 1e-9u"W/m^2"
    Kt    = G ./ max.(G0h, ϵ)
    Kt    = clamp.(Kt, 0.0, 1.2)
    Fd = similar(Kt) # diffuse fraction
    for i in eachindex(Kt)
        if Kt[i] <= 0.22
            Fd[i] = 1 - 0.09*Kt[i]
        elseif Kt[i] <= 0.80
            Fd[i] = 0.9511 - 0.1604*Kt[i] + 4.388*Kt[i]^2 - 16.638*Kt[i]^3 + 12.336*Kt[i]^4
        else
            Fd[i] = 0.165
        end
    end
    # TODO probably this still allocates because of aliasing 
    Fd .= clamp.(Fd, zero(eltype(Fd)), oneunit(eltype(Fd)))

    D .= Fd .* G
    B .= G .- D

    # Zero everything at night
    # night = (cosz_pos .== 0)
    # D[night] .= 0.0u"W/m^2"
    # B[night] .= 0.0u"W/m^2"
    # G[night] .= 0.0u"W/m^2"

    return output
end

function cloud_adjust_radiation(cloud::AbstractVector, args...; kw...)
    n = length(cloud)
    global_solar = fill(0.0u"W/m^2", n)
    diffuse_solar = fill(0.0u"W/m^2", n)
    direct_solar = fill(0.0u"W/m^2", n)
    output = (; global_solar, diffuse_solar, direct_solar)
    cloud_adjust_radiation!(output, cloud, args...; kw...)

    return output
end

# Separated out from dchxy for easier optimisation
# This algorithm is very expensive
# TODO these argument names are nightmare fuel
@noinline function _dchxy_converge!(FNX, FNY, AMU, PSI, XA, XB, XD, XE, CHX, CHY, CHXA, CHYA)
    nomitr = 1 # Fortran line 362
    TEMC = 0.0 # Initialize before convergence loop
    converged = false

    while !converged
        for I in 2:101
            fnx_i = FNX[I] 
            fny_i = FNY[I] 
            amu_i = AMU[I]

            #######################################################################################################
            # Compute XD and XE for this I
            # The most performance-intensive code of the package: loop inside loop inside while, called from another loop
            # Possibly there is a faster algorithm?
            # works marginally better when each line is separate
            for IC in 1:101
                XD[IC] = PSI[IC] * (fnx_i * FNX[IC] - fny_i * FNY[IC]) / (amu_i + AMU[IC])
            end
            for IC in 1:101
                XE[IC] = PSI[IC] * (fny_i * FNX[IC] - fnx_i * FNY[IC]) / (amu_i - AMU[IC])
            end
            #######################################################################################################

            # Everett's formula / interpolation for XE[I]
            XE[I] = if I <= 3
                0.5 * (XE[I+1] + XE[I-1])
            elseif I <= 5
                0.0625 * (9.0*(XE[I+1] + XE[I-1]) - XE[I+3] - XE[I-3])
            elseif I <= 96
                (3.0*(XE[I+5] + XE[I-5]) + 150.0*(XE[I+1] + XE[I-1]) - 25.0*(XE[I+3] + XE[I-3])) / 256.0
            else
                5.0*XE[I-1] + 10.0*XE[I-3] + XE[I-5] - 10.0*XE[I-2] - 5.0*XE[I-4]
            end


            #########################################################
            # Second most expensive code in the package
            # is a huge performance gain
            sxd = 0.0
            sxe = 0.0
            for ic in 1:101
                sxd += XA[ic] * XD[ic]
                sxe += XA[ic] * XE[ic]
            end
            #########################################################

            CHXA[I] = 1.0 + amu_i * sxd
            CHYA[I] = XB[I] + amu_i * sxe
        end

        # Correction to CHX and CHY
        for i in 1:101
            TEMD = TEMC * AMU[i] * (1.0 - XB[i])
            CHX[i] = CHXA[i] + TEMD
            CHY[i] = CHYA[i] + TEMD
        end

        # Check convergence (same as before)
        if nomitr > 1
            for I in 2:101
                rel_error = abs((CHY[I] - FNY[I]) / CHY[I])
                # TODO this seems wrong? shouldnt it only break if errors are <= 2.0e-4 for all I ?
                if rel_error <= 2.0e-4
                    converged = true
                    break
                end
            end
        end

        # Prepare for next iteration
        for I in 1:101
            FNX[I] = CHX[I]
            FNY[I] = CHY[I]
        end

        nomitr += 1
        nomitr > 15 && break
    end 

    return converged, nomitr
end
