const cp_water = 4184.0u"J/kg/K" # heat capacity of pure water
const ρ_water = 1000.0u"kg/m^3" # density of pure water
const ρ_hat0 = 44.65u"mol/m^3" # density of water vapour at STP, p. 309 Campbell et al 1994
const D_v0 = 2.12e-5u"m^2/s" # diffusivity of water vapour at STP, p. 309 Campbell et al 1994


const DEFAULT_HEIGHTS = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.2] .* u"m"
const DEFAULT_DEPTHS = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 100.0, 200.0]u"cm"

