using Microclimate
using Unitful
using Unitful: °, rad, R, kg, m
using Plots
using CSV, DataFrames
using Test

# NicheMapR simulation results
λDirect_NMR = DataFrame(CSV.File("tests/data/drlam_monthly.csv"))
λRayleigh_NMR = DataFrame(CSV.File("tests/data/drrlam_monthly.csv"))
λScattered_NMR = DataFrame(CSV.File("tests/data/srlam_monthly.csv"))
metout_NMR = DataFrame(CSV.File("tests/data/metout_monthly.csv"))

# NicheMapR simulation parameters
microinput_vec = DataFrame(CSV.File("tests/data/init_monthly/microinput.csv"))[:, 2]
names = [
    :doynum, :RUF, :ERR, :Usrhyt, :Refhyt, :Numtyps, :Z01, :Z02, :ZH1, :ZH2,
    :idayst, :ida, :HEMIS, :ALAT, :AMINUT, :ALONG, :ALMINT, :ALREF, :slope,
    :azmuth, :ALTT, :CMH2O, :microdaily, :tannul, :EC, :VIEWF, :snowtemp,
    :snowdens, :snowmelt, :undercatch, :rainmult, :runshade, :runmoist,
    :maxpool, :evenrain, :snowmodel, :rainmelt, :writecsv, :densfun1, :densfun2,
    :densfun3, :densfun4, :hourly, :rainhourly, :lamb, :IUV, :RW, :PC, :RL, :SP, :R1, 
    :IM, :MAXCOUNT, :IR, :message, :fail, :snowcond, :intercept, :grasshade, :solonly, 
    :ZH, :D0, :TIMAXS1, :TIMAXS2, :TIMAXS3, :TIMAXS4, :TIMINS1, :TIMINS2, :TIMINS3, :TIMINS4,
    :spinup, :dewrain, :moiststep, :maxsurf, :ndmax
]
microinput = (; zip(names, microinput_vec)...)

iuv = Bool(Int(microinput[:IUV])) # this makes it take ages if true!
lat = (microinput[:ALAT]+microinput[:AMINUT]/60)*1.0u"°" # latitude
slope = (microinput[:slope])*1.0u"°" # slope
aspect = (microinput[:azmuth])*1.0u"°" # aspect
elev = (microinput[:ALTT])*1.0u"m" # elevation
hori = (DataFrame(CSV.File("tests/data/init_monthly/hori.csv"))[:, 2])*1.0u"°"#fill(0.0u"°", 24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
REFLS = (DataFrame(CSV.File("tests/data/init_monthly/REFLS.csv"))[:, 2]*1.0) # set up vector of soil reflectances for each day (decimal %)
refl = REFLS[1]
CCMINN = (DataFrame(CSV.File("tests/data/init_monthly/CCMINN.csv"))[:, 2]/100.) # min cloud cover (dec%)
CCMAXX = (DataFrame(CSV.File("tests/data/init_monthly/CCMAXX.csv"))[:, 2]/100.) # max cloud cover (dec%)
cloud = CCMINN .+ CCMINN ./ 2

hours = collect(0.:1:24.)
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]*1.0

solrad_out = @inferred solrad(
    days = days,     # days of year
    hours = hours,   # hours of day
    lat = lat,       # latitude (degrees)
    elev = elev,     # elevation (m)
    hori = hori,     # horizon angles 0 degrees azimuth (north) clockwise in 15 degree intervals
    slope = slope,   # slope (degrees, range 0-90)
    aspect = aspect, # aspect (degrees, 0 = North, range 0-360)
    refl = refl,     # substrate solar reflectivity (decimal %)
    iuv = true        # use Dave_Furkawa theory for UV radiation (290-360 nm)?
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

kλ = solrad_out.λ
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


@testset "solar radiation comparisons" begin
    @test ustrip.(u"°", Zenith) ≈ metout_NMR.ZEN atol=1e-4
    @test all(isapprox.(ustrip.(u"W/m^2", Global), metout_NMR.SOLR; atol=1e-1))
    @test ustrip.(u"W/m^2", Global) ≈ metout_NMR.SOLR atol=1e-1 rtol=0.0
    @test λDirect[hour, :] ≈ λDirect_NMR_hr atol=1e-3u"W/nm/m^2"
    @test ustrip.(λScattered[hour, :]) ≈ ustrip.(λScattered_NMR_hr) atol=1e-7#u"W/nm/m^2"
    @test λRayleigh[hour, :] ≈ λRayleigh_NMR_hr atol=1e-3u"W/nm/m^2"
end  
    
diffglob = abs.(ustrip.(u"W/m^2", Global) .- metout_NMR.SOLR)
maxloc = findmax(diffglob)[2]
Global[maxloc]
metout_NMR.SOLR[maxloc]

# plot(hours, Zenith, xlabel="hour of day", ylabel="Zenith angle", legend=false)
# plot!([tsn - HHsr, tsn + HHsr], seriestype="vline", color="red", linestyle = [:dash, :dash])

# plot(hours, [Global Direct Scattered Rayleigh], xlabel="hour of day", ylabel="Radiation", label=["Global" "Direct" "Scattered" "Rayleigh"])

# hour = findfirst(x -> x == 6.0, hours)
# #plot(λ, [λGlobal[hour, :] λDirect[hour, :] λScattered[hour, :] λRayleigh[hour, :]], xlabel="Wavelength", ylabel="Spectral Irradiance", label=["Global" "Direct" "Scattered" "Rayleigh"])

# plot(λ, λRayleigh[hour, :], xlabel="Wavelength", ylabel="Spectral Irradiance", label="Rayleigh")
# plot!(λ, λGlobal[hour, :], xlabel="Wavelength", ylabel="Spectral Irradiance", label="Global")
# plot!(λ, λDirect[hour, :], xlabel="Wavelength", ylabel="Spectral Irradiance", label="Direct")
# plot!(λ, λScattered[hour, :], xlabel="Wavelength", ylabel="Spectral Irradiance", label="Scattered")