using Microclimate
using Unitful
using Unitful: °, rad, R, kg, m
using Plots
using CSV, DataFrames

λDirect_NMR = DataFrame(CSV.File("data/drlam.csv"))
λRayleigh_NMR = DataFrame(CSV.File("data/drrlam.csv"))
λScattered_NMR = DataFrame(CSV.File("data/srlam.csv"))
metout_NMR = DataFrame(CSV.File("data/metout.csv"))

hours = collect(0.:1:24.)
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
solrad_out = solrad(
    days = days,               # days of year
    hours = hours,              # hours of day
    lat = 43.1379°,             # latitude (degrees)
    elev = 226m,                  # elevation (m)
    hori = fill(0.0°, 24),      # horizon angles 0 degrees azimuth (north) clockwise in 15 degree intervals
    slope = 0.0°,               # slope (degrees, range 0-90)
    aspect = 0.0°,              # aspect (degrees, 0 = North, range 0-360)
    refl = 0.10,                # substrate solar reflectivity (decimal %)
    iuv = false                 # use Dave_Furkawa theory for UV radiation (290-360 nm)?
    )

# extract output
Zenith = solrad_out.Zenith
HHsr = solrad_out.HHsr
tsn = solrad_out.tsn
Global = solrad_out.Global
Scattered = solrad_out.Scattered
Direct = solrad_out.Direct
Rayleigh = solrad_out.Rayleigh
λ = solrad_out.λ
λGlobal = solrad_out.λGlobal
λScattered = solrad_out.λScattered
λDirect = solrad_out.λDirect
λRayleigh = solrad_out.λRayleigh

hours25 = repeat(hours, outer = length(days))
hours24 = repeat(0:1:23, outer = length(days))
hours_remove = findall(x->x==24.0, hours25)
Zenith = solrad_out.Zenith
Global = solrad_out.Global
deleteat!(Zenith, hours_remove)
deleteat!(Global, hours_remove)
plot(Zenith, ylabel="Zenith angle", legend=false)
plot!(metout_NMR.ZEN, linestyle = :dash)
plot(Global, ylabel="Radiation", label="solrad.jl")
plot!(metout_NMR.SOLR, linestyle = :dash, label="NMR")

λ = solrad_out.λ
λGlobal = solrad_out.λGlobal
λDirect = solrad_out.λDirect
λRayleigh = solrad_out.λRayleigh
λScattered = solrad_out.λScattered

hour2do = 9.0
hour = findfirst(x -> x == hour2do, hours25)
hourNMR = findfirst(x -> x == hour2do+0.0, hours24)
λDirect_NMR_hr = collect(λDirect_NMR[hourNMR, 4:114])u"W/m^2/nm"
λRayleigh_NMR_hr = collect(λRayleigh_NMR[hourNMR, 4:114])u"W/m^2/nm"
λScattered_NMR_hr = collect(λScattered_NMR[hourNMR, 4:114])u"W/m^2/nm"

plot(λ, [λDirect[hour, :] λScattered[hour, :] λRayleigh[hour, :]], xlabel="Wavelength", ylabel="Spectral Irradiance", label=["Direct" "Scattered" "Rayleigh"])
plot!(λ, [λDirect_NMR_hr λScattered_NMR_hr λRayleigh_NMR_hr], xlabel="Wavelength", ylabel="Spectral Irradiance", label=["Direct" "Scattered" "Rayleigh"], linestyle=[:dash :dash :dash])

# plot(hours, Zenith, xlabel="hour of day", ylabel="Zenith angle", legend=false)
# plot!([tsn - HHsr, tsn + HHsr], seriestype="vline", color="red", linestyle = [:dash, :dash])

# plot(hours, [Global Direct Scattered Rayleigh], xlabel="hour of day", ylabel="Radiation", label=["Global" "Direct" "Scattered" "Rayleigh"])

# hour = findfirst(x -> x == 6.0, hours)
# #plot(λ, [λGlobal[hour, :] λDirect[hour, :] λScattered[hour, :] λRayleigh[hour, :]], xlabel="Wavelength", ylabel="Spectral Irradiance", label=["Global" "Direct" "Scattered" "Rayleigh"])

# plot(λ, λRayleigh[hour, :], xlabel="Wavelength", ylabel="Spectral Irradiance", label="Rayleigh")
# plot!(λ, λGlobal[hour, :], xlabel="Wavelength", ylabel="Spectral Irradiance", label="Global")
# plot!(λ, λDirect[hour, :], xlabel="Wavelength", ylabel="Spectral Irradiance", label="Direct")
# plot!(λ, λScattered[hour, :], xlabel="Wavelength", ylabel="Spectral Irradiance", label="Scattered")