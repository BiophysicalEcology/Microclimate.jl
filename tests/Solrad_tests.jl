using Unitful
using Unitful: °, rad, R, kg
using ModelParameters
using Dates
using Plots
#@unit pct "pct" Percent 0.01 false

#days = collect(1.:365.) # day of year
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
#days = [196]
hours = collect(0.:1:24.) # hour of day
lat = -40° # latitude
elev = 0u"m" # elevation (m)
hori = fill(0.0°, 24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
slope = 90.0° # slope (degrees, range 0-90)
aspect = 0.0° # aspect (degrees, 0 = North, range 0-360)
refl = 0.10 # substrate solar reflectivity (decimal %)
iuv = false
#solrad_out = solrad()
#plot(solrad_out.Global)
#maximum(solrad_out.Global)
solrad_out = solrad(days = days, hours = hours, lat = lat, elev = elev, hori = hori, slope = slope, aspect = aspect, refl = refl, iuv = iuv)
#plot(solrad_out.λ, solrad_out.λGlobal[12, :])
plot(solrad_out.Global)
#plot!(solrad_out.Global)
#plot(solrad_out.Zenith)

profile = get_profile()
scatter(profile.heights, profile.VELs)
scatter(profile.heights, profile.TAs)
scatter(profile.heights, profile.RHs)
