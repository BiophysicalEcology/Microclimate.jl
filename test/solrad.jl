using Microclimate
using Microclimate: solar_geometry, McCulloughPorterSolarGeometry, hour_angle, SolarRadiation, solrad
using Unitful
using Unitful: °, rad, R, kg, m
using CSV, DataFrames, Statistics
using Test
sm = McCulloughPorterSolarGeometry()
h = hour_angle(12.0)[1]
latitude = [-29.0u"°", 29.0u"°"]
results = solar_geometry.(Ref(sm), latitude; d=1, h=h)
solar_model = SolarRadiation()
latitude = [-29.0u"°", 29.0u"°"]
elevation = [1000.0u"m", 1.0u"m"]
horizon_angles = fill(1.0u"°", 24)
slope = [20.0u"°", 0.0u"°"]
aspect = [0.0u"°", 0.0u"°"]
P_atmos = [101325.0u"Pa", 101325.0u"Pa"]
albedo = [0.2, 0.45]
solrad_out = solrad.(Ref(solar_model), latitude, elevation, 
    slope, aspect, P_atmos, albedo;
    days=1,               # days of year
    hours=12,              # hours of day
    horizon_angles, 
)

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
albedo = (DataFrame(CSV.File("$testdir/data/init_monthly/REFLS.csv"))[:, 2]*1.0) # set up vector of soil albedoectances for each day (decimal %)
τA_nmr = (DataFrame(CSV.File("$testdir/data/init_monthly/TAI.csv"))[:, 2]*1.0)
minima_times = [microinput[:TIMINS1], microinput[:TIMINS2], microinput[:TIMINS3], microinput[:TIMINS4]] # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
maxima_times = [microinput[:TIMAXS1], microinput[:TIMAXS2], microinput[:TIMAXS3], microinput[:TIMAXS4]] # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
cloud_min = (DataFrame(CSV.File("$testdir/data/init_monthly/CCMINN.csv"))[:, 2] * 1.0) # min cloud cover (%)
cloud_max = (DataFrame(CSV.File("$testdir/data/init_monthly/CCMAXX.csv"))[:, 2] * 1.0) # max cloud cover (c%)

hours = collect(0.0:1:23.0)
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]*1.0

solar_model = SolarRadiation(; iuv = Bool(Int(microinput[:IUV])))

terrain = Terrain(;
    slope = (microinput[:slope])*1.0u"°",
    aspect = (microinput[:azmuth])*1.0u"°",
    elevation = (microinput[:ALTT])*1.0u"m",
    horizon_angles = (DataFrame(CSV.File("$testdir/data/init_monthly/hori.csv"))[:, 2])*1.0u"°", #fill(0.0u"°", 24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
)

@time solrad_out = solrad(solar_model;
    days,               # days of year
    hours,              # hours of day
    latitude,           # latitude (degrees)
    terrain,
    albedo,            # substrate solar albedoectivity (decimal %)
);

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

hour_angle_sunrise = solrad_out.hour_angle_sunrise
hour_solar_noon = solrad_out.hour_solar_noon
global_total = solrad_out.global_total
direct_total = solrad_out.direct_total
diffuse_total = solrad_out.diffuse_total
rayleigh_total = solrad_out.rayleigh_total
# Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
day_of_year = repeat(days, inner=length(hours))
global_cloud, diffuse_cloud, direct_cloud = cloud_adjust_radiation(cloud_covers / 100., diffuse_total, direct_total, zenith_angle, day_of_year)
#global_cloud = global_total .* (0.36 .+ 0.64 * (1.0 .- (cloud_covers ./ 100.0))) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992

# diffuse spectra test needs to be 1e-2 to pass with iuv=true
# global_cloud out by a tiny amount, < 0.5 W/nm/m2
@testset "solar radiation comparisons" begin
    @test ustrip.(u"°", zenith_angle) ≈ metout_nmr.ZEN atol=1e-4
    @test all(isapprox.(ustrip.(u"W/m^2", global_cloud), metout_nmr.SOLR; atol=0.5))
    @test solrad_out.direct_spectra ≈ direct_spectra_nmr_units atol=1e-4u"W/nm/m^2"
    @test solrad_out.diffuse_spectra ≈ diffuse_spectra_nmr_units atol=1e-2u"W/nm/m^2"
    @test solrad_out.rayleigh_spectra ≈ rayleigh_spectra_nmr_units atol=1e-4u"W/nm/m^2"
end  
