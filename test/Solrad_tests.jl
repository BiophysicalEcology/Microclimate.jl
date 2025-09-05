using Microclimate
using Unitful
using Unitful: °, rad, R, kg, m
using Plots
using CSV, DataFrames, Polynomials, Statistics
using Test

# TODO - compare with McCullough and Porter plots and Insolation.jl and CloudScat.jl? 
# TODO check how CliMA do atmospheric response - looks like it's in RRTMGP.jl

testdir = realpath(joinpath(dirname(pathof(Microclimate)), "../test"))

# NicheMapR simulation results
direct_spectra_nmr = Matrix(DataFrame(CSV.File("$testdir/data/drlam_monthly.csv"))[:, 4:114])
rayleigh_spectra_nmr = Matrix(DataFrame(CSV.File("$testdir/data/drrlam_monthly.csv"))[:, 4:114])
diffuse_spectra_nmr = Matrix(DataFrame(CSV.File("$testdir/data/srlam_monthly.csv"))[:, 4:114])
metout_nmr = DataFrame(CSV.File("$testdir/data/metout_monthly.csv"))
direct_spectra_nmr_units = direct_spectra_nmr*u"W/m^2/nm"
rayleigh_spectra_nmr_units = rayleigh_spectra_nmr*u"W/m^2/nm"
diffuse_spectra_nmr_units = diffuse_spectra_nmr*u"W/m^2/nm"

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

longlat = (DataFrame(CSV.File("$testdir/data/init_monthly/longlat.csv"))[:, 2] * 1.0)
latitude = longlat[2]*1.0u"°" # latitude
longitude =  longlat[1]*1.0u"°" # longitude
slope = (microinput[:slope])*1.0u"°" # slope
aspect = (microinput[:azmuth])*1.0u"°" # aspect
elevation = (microinput[:ALTT])*1.0u"m" # elevation
horizon_angles = (DataFrame(CSV.File("$testdir/data/init_monthly/HORI.csv"))[:, 2])*1.0u"°"#fill(0.0u"°", 24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
albedos = (DataFrame(CSV.File("$testdir/data/init_monthly/REFLS.csv"))[:, 2]*1.0) # set up vector of soil albedoectances for each day (decimal %)
τA_nmr = (DataFrame(CSV.File("$testdir/data/init_monthly/TAI.csv"))[:, 2]*1.0)
minima_times = [microinput[:TIMINS1], microinput[:TIMINS2], microinput[:TIMINS3], microinput[:TIMINS4]] # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
maxima_times = [microinput[:TIMAXS1], microinput[:TIMAXS2], microinput[:TIMAXS3], microinput[:TIMAXS4]] # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
cloud_min = (DataFrame(CSV.File("$testdir/data/init_monthly/CCMINN.csv"))[:, 2] * 1.0) # min cloud cover (%)
cloud_max = (DataFrame(CSV.File("$testdir/data/init_monthly/CCMAXX.csv"))[:, 2] * 1.0) # max cloud cover (c%)
iuv = Bool(Int(microinput[:IUV])) # this makes it take ages if true!

hours = collect(0.:1:24.)
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]*1.0

# test gads function
# run gads and get mean interpolated output over seasons
optdep_summer = gads(ustrip(latitude), ustrip(longitude), 1, 0)
optdep_winter = gads(ustrip(latitude), ustrip(longitude), 1, 1)
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
plot!(λ, τA_nmr, linecolor="grey")

solrad_out = @inferred solrad(;
    days,               # days of year
    hours,              # hours of day
    latitude,           # latitude (degrees)
    elevation,          # elevation (m)
    horizon_angles,     # horizon angles 0 degrees azimuth (north) clockwise in 15 degree intervals
    slope,              # slope (degrees, range 0-90)
    aspect,             # aspect (degrees, 0 = North, range 0-360)
    albedos,            # substrate solar albedoectivity (decimal %)
    iuv,                # use Dave_Furkawa theory for UV radiation (290-360 nm)?
    τA                  # aerosol profile from gads (global aerosol data set)
    )

# using ProfileView
# ProfileView.@profview solrad_out = @inferred solrad(;
#     days,               # days of year
#     hours,              # hours of day
#     latitude,           # latitude (degrees)
#     elevation,          # elevation (m)
#     horizon_angles,     # horizon angles 0 degrees azimuth (north) clockwise in 15 degree intervals
#     slope,              # slope (degrees, range 0-90)
#     aspect,             # aspect (degrees, 0 = North, range 0-360)
#     albedos,            # substrate solar albedoectivity (decimal %)
#     iuv = true,                # use Dave_Furkawa theory for UV radiation (290-360 nm)?
#     τA                  # aerosol profile from gads (global aerosol data set)
#     )

cloud_covers = hourly_vars(
    cloud_min,
    cloud_max,
    solrad_out,
    minima_times,
    maxima_times,
    false
)
    
# extract output
zenith_angle = solrad_out.zenith_angle
zenith_angle[zenith_angle.>90u"°"] .= 90u"°"
azimuth_angle = solrad_out.azimuth_angle
plot(azimuth_angle)
plot(zenith_angle)
hour_angle_sunrise = solrad_out.hour_angle_sunrise
hour_solar_noon = solrad_out.hour_solar_noon
global_total = solrad_out.global_total
direct_total = solrad_out.direct_total
diffuse_total = solrad_out.diffuse_total
rayleigh_total = solrad_out.rayleigh_total
# Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
day_of_year = repeat(days, inner=length(hours))
global_cloud, diffuse_cloud, direct_cloud = cloud_adjust_radiation(cloud_covers / 100., diffuse_total, direct_total, zenith_angle, day_of_year)
#Global = Global .* (0.36 .+ 0.64 * (1.0 .- (CLDs / 100.0))) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992

wavelength = solrad_out.wavelength
global_spectra = solrad_out.global_spectra
diffuse_spectra = solrad_out.diffuse_spectra
direct_spectra = solrad_out.direct_spectra
rayleigh_spectra = solrad_out.rayleigh_spectra

hours25 = repeat(hours, outer = length(days))
hours24 = repeat(0:1:23, outer = length(days))
hours_remove = findall(x->x==24.0, hours25)
deleteat!(zenith_angle, hours_remove)
deleteat!(global_total, hours_remove)
deleteat!(global_cloud, hours_remove)
deleteat!(direct_total, hours_remove)
deleteat!(diffuse_total, hours_remove)
rows_to_remove = 25:25:300
global_spectra = global_spectra[setdiff(1:end, rows_to_remove), :]
diffuse_spectra = diffuse_spectra[setdiff(1:end, rows_to_remove), :]
direct_spectra = direct_spectra[setdiff(1:end, rows_to_remove), :]
rayleigh_spectra = rayleigh_spectra[setdiff(1:end, rows_to_remove), :]
plot(zenith_angle, ylabel="Zenith angle", legend=false)
plot!(metout_nmr.ZEN, linestyle = :dash)
plot(global_cloud, ylabel="Radiation", label="solrad.jl")
plot!(metout_nmr.SOLR, linestyle = :dash, label="NMR")

month2do = 6
hour2do = 12
i = (month2do - 1) * 24 + hour2do

plot(λ, [direct_spectra[i, :] diffuse_spectra[i, :] rayleigh_spectra[i, :]], xlabel="Wavelength", ylabel="Spectral Irradiance", label=["Direct" "Diffuse" "Rayleigh"])
plot!(λ, [direct_spectra_nmr_units[i, :] diffuse_spectra_nmr_units[i, :] rayleigh_spectra_nmr_units[i, :]], xlabel="Wavelength", ylabel="Spectral Irradiance", label=["Direct" "Diffuse" "Rayleigh"], linestyle=[:dash :dash :dash])

@testset "solar radiation comparisons" begin
    @test τA ≈ τA_nmr atol=1e-7
    @test ustrip.(u"°", Zenith) ≈ metout_nmr.ZEN atol=1e-4
    @test all(isapprox.(ustrip.(u"W/m^2", solrad_out.global_total), metout_nmr.SOLR; atol=0.5))
    @test direct_spectra ≈ direct_spectra_nmr_units atol=1e-4u"W/nm/m^2"
    @test diffuse_spectra ≈ diffuse_spectra_nmr_units atol=1e-6u"W/nm/m^2"
    @test rayleigh_spectra ≈ rayleigh_spectra_nmr_units atol=1e-4u"W/nm/m^2"
end  
