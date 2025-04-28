using Microclimate
using Unitful
using Unitful: °, rad, R, kg, m
using Plots

hours = collect(0.:0.1:24.)

solrad_out = solrad(
    days = [180],               # days of year
    hours = hours,              # hours of day
    lat = 43.1379°,             # latitude (degrees)
    elev = 0m,                  # elevation (m)
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

plot(hours, Zenith, xlabel="hour of day", ylabel="Zenith angle", legend=false)
plot!([tsn - HHsr, tsn + HHsr], seriestype="vline", color="red", linestyle = [:dash, :dash])

plot(hours, [Global Direct Scattered Rayleigh], xlabel="hour of day", ylabel="Radiation", label=["Global" "Direct" "Scattered" "Rayleigh"])

hour = findfirst(x -> x == 6.0, hours)
#plot(λ, [λGlobal[hour, :] λDirect[hour, :] λScattered[hour, :] λRayleigh[hour, :]], xlabel="Wavelength", ylabel="Spectral Irradiance", label=["Global" "Direct" "Scattered" "Rayleigh"])

plot(λ, λRayleigh[hour, :], xlabel="Wavelength", ylabel="Spectral Irradiance", label="Rayleigh")
plot!(λ, λGlobal[hour, :], xlabel="Wavelength", ylabel="Spectral Irradiance", label="Global")
plot!(λ, λDirect[hour, :], xlabel="Wavelength", ylabel="Spectral Irradiance", label="Direct")
plot!(λ, λScattered[hour, :], xlabel="Wavelength", ylabel="Spectral Irradiance", label="Scattered")