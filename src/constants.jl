const water_heat_capacity = 4184.0u"J/kg/K" # heat capacity of pure water
const water_density = 1000.0u"kg/m^3" # density of pure water
const water_vapour_molar_density_stp = 44.65u"mol/m^3" # molar density of water vapour at STP, p. 309 Campbell et al 1994
const water_vapour_diffusivity_stp = 2.12e-5u"m^2/s" # diffusivity of water vapour at STP, p. 309 Campbell et al 1994


"""
    DEFAULT_HEIGHTS

Default air nodes (m) for the temperature/wind/humidity profile.
"""
const DEFAULT_HEIGHTS = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.2] .* u"m"

"""
    DEFAULT_DEPTHS

Default 19 soil depth nodes (cm), spaced closer near the surface where temperature
and moisture change fastest (see [Soil hydraulics algorithms](@ref)).
"""
const DEFAULT_DEPTHS = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"

"""
    DEFAULT_DAYS

Default day-of-year for one representative day per month (the 15th of each,
approximately), for use with [`MonthlyMinMaxEnvironment`](@ref) and
[`NonConsecutiveDayMode`](@ref).
"""
const DEFAULT_DAYS = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]

"""
    DEFAULT_HOURS

Default hour-of-day grid for solar radiation.
"""
const DEFAULT_HOURS = collect(0.0:1:23.0)

const LATENT_HEAT_FUSION = 333550.0u"J/kg"

"""
    DEFAULT_SNOW_NODE_THRESHOLDS

Default depths (cm) of the up-to-8 snow nodes [`SnowModel`](@ref) activates as the
pack grows.
"""
const DEFAULT_SNOW_NODE_THRESHOLDS = (2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 300.0)

# Bloomsburg & Wang (1969, Soil Sci. Soc. Am. J. 33:686-691), via SHAW's SOILHK:
# conductivity blocked below this ice-free porosity (m3/m3)
const ICE_IMPEDANCE_MIN_POROSITY = 0.13
# Numerical regularization, not part of Bloomsburg & Wang: floor so a fully
# frozen layer's conductivity stays nonzero, avoiding 0/0 in root-uptake resistance
const ICE_CONDUCTIVITY_FLOOR_FACTOR = 1e-6

