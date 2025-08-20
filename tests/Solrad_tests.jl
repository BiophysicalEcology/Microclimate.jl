using Microclimate
using Unitful
using Unitful: °, rad, R, kg, m
using Plots
using CSV, DataFrames
using Test

# to do - compare with McCullough and Porter plots and Insolation.jl and CloudScat.jl? 
# check how CliMA do atmospheric response - looks like it's in RRTMGP.jl

# NicheMapR simulation results
λDirect_NMR = Matrix(DataFrame(CSV.File("tests/data/drlam_monthly.csv"))[:, 4:114])
λRayleigh_NMR = Matrix(DataFrame(CSV.File("tests/data/drrlam_monthly.csv"))[:, 4:114])
λScattered_NMR = Matrix(DataFrame(CSV.File("tests/data/srlam_monthly.csv"))[:, 4:114])
metout_NMR = DataFrame(CSV.File("tests/data/metout_monthly.csv"))
λDirect_NMR_units = λDirect_NMR*u"W/m^2/nm"
λRayleigh_NMR_units = λRayleigh_NMR*u"W/m^2/nm"
λScattered_NMR_units = λScattered_NMR*u"W/m^2/nm"

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
TIMINS = [microinput[:TIMINS1], microinput[:TIMINS2], microinput[:TIMINS3], microinput[:TIMINS4]] # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
TIMAXS = [microinput[:TIMAXS1], microinput[:TIMAXS2], microinput[:TIMAXS3], microinput[:TIMAXS4]] # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
CCMINN = (DataFrame(CSV.File("tests/data/init_monthly/CCMINN.csv"))[:, 2] * 1.0) # min cloud cover (%)
CCMAXX = (DataFrame(CSV.File("tests/data/init_monthly/CCMAXX.csv"))[:, 2] * 1.0) # max cloud cover (c%)
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
    iuv = iuv        # use Dave_Furkawa theory for UV radiation (290-360 nm)?
    )
CLDs = hourly_vars(
    CCMINN,
    CCMAXX,
    solrad_out,
    TIMINS,
    TIMAXS,
    false
)
    
# extract output
Zenith = solrad_out.Zenith
Zenith[Zenith.>90u"°"] .= 90u"°"
Azimuth = solrad_out.Azimuth
plot(Azimuth)
plot(Zenith)
HHsr = solrad_out.HHsr
tsn = solrad_out.tsn
Global = solrad_out.Global
# Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
Global = Global .* (0.36 .+ 0.64 * (1.0 .- (CLDs / 100.0))) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
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
deleteat!(Zenith, hours_remove)
deleteat!(Global, hours_remove)
deleteat!(Direct, hours_remove)
deleteat!(Scattered, hours_remove)
rows_to_remove = 25:25:300
λGlobal = λGlobal[setdiff(1:end, rows_to_remove), :]
λScattered = λScattered[setdiff(1:end, rows_to_remove), :]
λDirect = λDirect[setdiff(1:end, rows_to_remove), :]
λRayleigh = λRayleigh[setdiff(1:end, rows_to_remove), :]
plot(Zenith, ylabel="Zenith angle", legend=false)
plot!(metout_NMR.ZEN, linestyle = :dash)
plot(Global, ylabel="Radiation", label="solrad.jl")
plot!(metout_NMR.SOLR, linestyle = :dash, label="NMR")

month2do = 6
hour2do = 12
i = (month2do - 1) * 24 + hour2do

plot(λ, [λDirect[i, :] λScattered[i, :] λRayleigh[i, :]], xlabel="Wavelength", ylabel="Spectral Irradiance", label=["Direct" "Scattered" "Rayleigh"])
plot!(λ, [λDirect_NMR_units[i, :] λScattered_NMR_units[i, :] λRayleigh_NMR_units[i, :]], xlabel="Wavelength", ylabel="Spectral Irradiance", label=["Direct" "Scattered" "Rayleigh"], linestyle=[:dash :dash :dash])

@testset "solar radiation comparisons" begin
    @test ustrip.(u"°", Zenith) ≈ metout_NMR.ZEN atol=1e-4
    @test all(isapprox.(ustrip.(u"W/m^2", Global), metout_NMR.SOLR; atol=0.5))
    @test λDirect ≈ λDirect_NMR_units atol=1e-4u"W/nm/m^2"
    @test λScattered ≈ λScattered_NMR_units atol=1e-6u"W/nm/m^2"
    @test λRayleigh ≈ λRayleigh_NMR_units atol=1e-4u"W/nm/m^2"
end  
