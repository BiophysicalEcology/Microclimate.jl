using Microclimate
using Unitful
using Unitful: °, rad, R, kg, m
using Plots
using CSV, DataFrames, Polynomials, Statistics
using Test

# to do - compare with McCullough and Porter plots and Insolation.jl and CloudScat.jl? 
# check how CliMA do atmospheric response - looks like it's in RRTMGP.jl

testdir = realpath(joinpath(dirname(pathof(Microclimate)), "../test"))

# NicheMapR simulation results
λDirect_NMR = Matrix(DataFrame(CSV.File("$testdir/data/drlam_monthly.csv"))[:, 4:114])
λRayleigh_NMR = Matrix(DataFrame(CSV.File("$testdir/data/drrlam_monthly.csv"))[:, 4:114])
λScattered_NMR = Matrix(DataFrame(CSV.File("$testdir/data/srlam_monthly.csv"))[:, 4:114])
metout_NMR = DataFrame(CSV.File("$testdir/data/metout_monthly.csv"))
λDirect_NMR_units = λDirect_NMR*u"W/m^2/nm"
λRayleigh_NMR_units = λRayleigh_NMR*u"W/m^2/nm"
λScattered_NMR_units = λScattered_NMR*u"W/m^2/nm"

# NicheMapR simulation parameters
microinput_vec = DataFrame(CSV.File("$testdir/data/init_monthly/microinput.csv"))[:, 2]
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
longlat = (DataFrame(CSV.File("$testdir/data/init_monthly/longlat.csv"))[:, 2] * 1.0)
lat = longlat[2]*1.0u"°" # latitude
lon =  longlat[1]*1.0u"°" # longitude
slope = (microinput[:slope])*1.0u"°" # slope
aspect = (microinput[:azmuth])*1.0u"°" # aspect
elevation = (microinput[:ALTT])*1.0u"m" # elevation
horizon_angles = (DataFrame(CSV.File("$testdir/data/init_monthly/horizon_angles.csv"))[:, 2])*1.0u"°"#fill(0.0u"°", 24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
albedos = (DataFrame(CSV.File("$testdir/data/init_monthly/REFLS.csv"))[:, 2]*1.0) # set up vector of soil albedoectances for each day (decimal %)
τA_NMR = (DataFrame(CSV.File("$testdir/data/init_monthly/TAI.csv"))[:, 2]*1.0)
TIMINS = [microinput[:TIMINS1], microinput[:TIMINS2], microinput[:TIMINS3], microinput[:TIMINS4]] # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
TIMAXS = [microinput[:TIMAXS1], microinput[:TIMAXS2], microinput[:TIMAXS3], microinput[:TIMAXS4]] # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
CCMINN = (DataFrame(CSV.File("$testdir/data/init_monthly/CCMINN.csv"))[:, 2] * 1.0) # min cloud cover (%)
CCMAXX = (DataFrame(CSV.File("$testdir/data/init_monthly/CCMAXX.csv"))[:, 2] * 1.0) # max cloud cover (c%)
cloud = CCMINN .+ CCMINN ./ 2

hours = collect(0.:1:24.)
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]*1.0

# test gads function
# run gads and get mean interpolated output over seasons
optdep_summer = gads(ustrip(lat), ustrip(lon), 1, 0)
optdep_winter = gads(ustrip(lat), ustrip(lon), 1, 1)
optdep_array = hcat(
    optdep_winter[:, 1], 
    mean.(eachrow(hcat(optdep_summer[:, 2], optdep_winter[:, 2])))
)
optdep = DataFrame(LAMBDA = optdep_array[:, 1], OPTDEPTH = optdep_array[:, 2])
xs = optdep.LAMBDA
ys = optdep.OPTDEPTH
xmin, xmax = extrema(xs)  # get min and max
# Scale xs to [-1,1] (can't fit higher order polynomials than 4 otherwise)
scale_xs(x) = 2 * (x - xmin) / (xmax - xmin) - 1
xscaled = scale_xs.(xs)
# fit polynomial in scaled space
p_scaled = Polynomials.fit(xscaled, ys, 6)
# function to evaluate the fit at original coordinates
p_eval(x) = p_scaled(scale_xs(x))
λ = float.([ # wavelengths across which to integrate
        290, 295, 300, 305, 310, 315, 320, 330, 340, 350, 360, 370, 380, 390,
        400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700,
        720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020,
        1080, 1100, 1120, 1140, 1160, 1180, 1200, 1220, 1240, 1260, 1280, 1300, 1320,
        1380, 1400, 1420, 1440, 1460, 1480, 1500, 1540, 1580, 1600, 1620, 1640, 1660,
        1700, 1720, 1780, 1800, 1860, 1900, 1950, 2000, 2020, 2050, 2100, 2120, 2150,
        2200, 2260, 2300, 2320, 2350, 2380, 2400, 2420, 2450, 2490, 2500, 2600, 2700,
        2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000
    ])
τA = p_eval.(λ)
plot(λ, τA)
plot!(λ, τA_NMR, linecolor="grey")

solrad_out = @inferred solrad(;
    days,       # days of year
    hours,      # hours of day
    lat,        # latitude (degrees)
    elevation,       # elevation (m)
    horizon_angles,       # horizon angles 0 degrees azimuth (north) clockwise in 15 degree intervals
    slope,      # slope (degrees, range 0-90)
    aspect,     # aspect (degrees, 0 = North, range 0-360)
    albedos,      # substrate solar albedoectivity (decimal %)
    iuv,        # use Dave_Furkawa theory for UV radiation (290-360 nm)?
    τA          # aerosol profile from gads (global aerosol database)
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
Direct = solrad_out.Direct
Diffuse = solrad_out.Scattered
# Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
doy  = repeat(days, inner=length(hours))
globalcloud, diffusecloud, directcloud = cloud_adjust_radiation(CLDs / 100., Diffuse, Direct, Zenith, doy)
#Global = Global .* (0.36 .+ 0.64 * (1.0 .- (CLDs / 100.0))) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992

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
deleteat!(globalcloud, hours_remove)
deleteat!(Direct, hours_remove)
deleteat!(Scattered, hours_remove)
rows_to_remove = 25:25:300
λGlobal = λGlobal[setdiff(1:end, rows_to_remove), :]
λScattered = λScattered[setdiff(1:end, rows_to_remove), :]
λDirect = λDirect[setdiff(1:end, rows_to_remove), :]
λRayleigh = λRayleigh[setdiff(1:end, rows_to_remove), :]
plot(Zenith, ylabel="Zenith angle", legend=false)
plot!(metout_NMR.ZEN, linestyle = :dash)
plot(globalcloud, ylabel="Radiation", label="solrad.jl")
plot!(metout_NMR.SOLR, linestyle = :dash, label="NMR")

month2do = 6
hour2do = 12
i = (month2do - 1) * 24 + hour2do

plot(λ, [λDirect[i, :] λScattered[i, :] λRayleigh[i, :]], xlabel="Wavelength", ylabel="Spectral Irradiance", label=["Direct" "Scattered" "Rayleigh"])
plot!(λ, [λDirect_NMR_units[i, :] λScattered_NMR_units[i, :] λRayleigh_NMR_units[i, :]], xlabel="Wavelength", ylabel="Spectral Irradiance", label=["Direct" "Scattered" "Rayleigh"], linestyle=[:dash :dash :dash])

@testset "solar radiation comparisons" begin
    @test τA ≈ τA_NMR atol=1e-7
    @test ustrip.(u"°", Zenith) ≈ metout_NMR.ZEN atol=1e-4
    @test all(isapprox.(ustrip.(u"W/m^2", Global), metout_NMR.SOLR; atol=0.5))
    @test λDirect ≈ λDirect_NMR_units atol=1e-4u"W/nm/m^2"
    @test λScattered ≈ λScattered_NMR_units atol=1e-6u"W/nm/m^2"
    @test λRayleigh ≈ λRayleigh_NMR_units atol=1e-4u"W/nm/m^2"
end  
