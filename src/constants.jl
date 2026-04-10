const water_heat_capacity = 4184.0u"J/kg/K" # heat capacity of pure water
const water_density = 1000.0u"kg/m^3" # density of pure water
const water_vapour_molar_density_stp = 44.65u"mol/m^3" # molar density of water vapour at STP, p. 309 Campbell et al 1994
const water_vapour_diffusivity_stp = 2.12e-5u"m^2/s" # diffusivity of water vapour at STP, p. 309 Campbell et al 1994


const DEFAULT_HEIGHTS = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.2] .* u"m"
const DEFAULT_DEPTHS = [0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]u"cm"

const LATENT_HEAT_FUSION = 333550.0u"J/kg"
const DEFAULT_SNOW_NODE_THRESHOLDS = (2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 300.0)

