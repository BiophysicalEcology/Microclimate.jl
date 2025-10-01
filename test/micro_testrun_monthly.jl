using Microclimate
using Unitful
using CSV, DataFrames
using Test

testdir = realpath(joinpath(dirname(pathof(Microclimate)), "../test"))

# read in output from NicheMapR and input variables
soiltemps_nmr = (DataFrame(CSV.File("$testdir/data/soil_monthly.csv"))[:, 4:13]) .* u"°C"
metout_nmr = DataFrame(CSV.File("$testdir/data/metout_monthly.csv"))
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

# Zip into a NamedTuple
microinput = (; zip(names, microinput_vec)...)

longlat = (DataFrame(CSV.File("$testdir/data/init_monthly/longlat.csv"))[:, 2] * 1.0)
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
LAIs = fill(0.1, length(days))
depths = ((DataFrame(CSV.File("$testdir/data/init_monthly/DEP.csv"))[:, 2]) / 100.0)u"m"
heights = [microinput[:Usrhyt], microinput[:Refhyt]]u"m" # air nodes for temperature, wind speed and humidity profile
days2do = 1:12
keywords = (;
    # locations, times, depths and heights
    latitude = longlat[2]*1.0u"°",
    days = days[days2do], # days of year for solrad
    hours = collect(0.:1:24.), # hour of day for solrad
    depths,
    heights, # air nodes for temperature, wind speed and humidity profile
    # terrain
    elevation = microinput[:ALTT] * 1.0u"m",
    horizon_angles = horizon_angles = (DataFrame(CSV.File("$testdir/data/init_monthly/hori.csv"))[:, 2]) * 1.0u"°",
    slope = microinput[:slope] * 1.0u"°",
    aspect = microinput[:azmuth] * 1.0u"°",
    roughness_height = microinput[:RUF] * 1.0u"m", # roughness height for standard mode TODO dispatch based on roughness pars
    zh = microinput[:ZH] * 1.0u"m", # heat transfer roughness height for Campbell and Norman mode
    d0 = microinput[:D0] * 1.0u"m", # zero plane displacement correction factor
    # soil thermal parameters 
    soil_mineral_conductivity = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][4]) * 1.0u"W/m/K", # soil minerals thermal conductivity (W/mC)
    soil_mineral_density = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][6]) * 1.0u"Mg/m^3", # soil minerals density (Mg/m3)
    soil_mineral_heat_capacity = c_p_m = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][5]) * 1.0u"J/kg/K", # soil minerals specific heat (J/kg-K)
    soil_bulk_density = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][2]) * 1.0u"Mg/m^3", # dry soil bulk density (Mg/m3)
    soil_saturation_moisture = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][3]) * 1.0u"m^3/m^3", # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    # daily environmental vectors
    albedos = (DataFrame(CSV.File("$testdir/data/init_monthly/REFLS.csv"))[days2do, 2] * 1.0), # substrate albedo (decimal %)
    shades = (DataFrame(CSV.File("$testdir/data/init_monthly/Minshades.csv"))[days2do, 2] * 1.0), # daily shade from vegetation (%)
    pctwets = (DataFrame(CSV.File("$testdir/data/init_monthly/PCTWET.csv"))[days2do, 2] * 1.0),
    sles = (DataFrame(CSV.File("$testdir/data/init_monthly/SLES.csv"))[days2do, 2] * 1.0), # - surface emissivity
    daily_rainfall = ((DataFrame(CSV.File("$testdir/data/init_monthly/rain.csv"))[days2do, 2] * 1.0) / 1000)u"kg/m^2", # monthly total rainfall
    air_temperature_min = (DataFrame(CSV.File("$testdir/data/init_monthly/TMINN.csv"))[days2do, 2] * 1.0)u"°C", # minimum air temperatures
    air_temperature_max = (DataFrame(CSV.File("$testdir/data/init_monthly/TMAXX.csv"))[days2do, 2] * 1.0)u"°C", # maximum air temperatures
    wind_min = (DataFrame(CSV.File("$testdir/data/init_monthly/WNMINN.csv"))[days2do, 2] * 1.0)u"m/s", # min wind speed (m/s)
    wind_max = (DataFrame(CSV.File("$testdir/data/init_monthly/WNMAXX.csv"))[days2do, 2] * 1.0)u"m/s", # max wind speed (m/s)
    humidity_min = (DataFrame(CSV.File("$testdir/data/init_monthly/RHMINN.csv"))[days2do, 2] * 1.0), # min relative humidity (%)
    humidity_max = (DataFrame(CSV.File("$testdir/data/init_monthly/RHMAXX.csv"))[days2do, 2] * 1.0), # max relative humidity (%)
    cloud_min = (DataFrame(CSV.File("$testdir/data/init_monthly/CCMINN.csv"))[days2do, 2] * 1.0), # min cloud cover (%)
    cloud_max = (DataFrame(CSV.File("$testdir/data/init_monthly/CCMAXX.csv"))[days2do, 2] * 1.0), # max cloud cover (%)
    minima_times = [microinput[:TIMINS1], microinput[:TIMINS2], microinput[:TIMINS3], microinput[:TIMINS4]], # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
    maxima_times = [microinput[:TIMAXS1], microinput[:TIMAXS2], microinput[:TIMAXS3], microinput[:TIMAXS4]], # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
    deep_soil_temperatures = (DataFrame(CSV.File("$testdir/data/init_monthly/tannulrun.csv"))[days2do, 2] * 1.0)u"°C", # daily deep soil temperatures
    # intial conditions
    initial_soil_temperature = u"K".((DataFrame(CSV.File("$testdir/data/init_monthly/soilinit.csv"))[1:length(depths), 2] * 1.0)u"°C"), # initial soil temperature
    initial_soil_moisture = (Array(DataFrame(CSV.File("$testdir/data/init_monthly/moists.csv"))[1, 2:13]) .* 1.0), # initial soil moisture
    leaf_area_index = fill(0.1, length(days)),
    iterate_day = microinput[:ndmax], # number of iterations per day
    daily = Bool(Int(microinput[:microdaily])), # doing consecutive days?
    runmoist = Bool(Int(microinput[:runmoist])), # run soil moisture algorithm?
    spinup = Bool(Int(microinput[:spinup])), # spin-up the first day by iterate_day iterations?
    iuv = Bool(Int(microinput[:IUV])), # this makes it take ages if true!
    #maximum_surface_temperature = u"K"(microinput[:maxsurf]u"°C")
);

# TODO check why deep soil temp not being outputted
@time micro_out = runmicro(; keywords...);

# subset NicheMapR predictions
vel1cm_nmr = collect(metout_nmr[:, 8]) .* 1u"m/s"
vel2m_nmr = collect(metout_nmr[:, 9]) .* 1u"m/s"
ta1cm_nmr = collect(metout_nmr[:, 4] .+ 273.15) .* 1u"K"
ta2m_nmr = collect(metout_nmr[:, 5] .+ 273.15) .* 1u"K"
rh1cm_nmr = collect(metout_nmr[:, 6])
rh2m_nmr = collect(metout_nmr[:, 7])
tskyC_nmr = collect(metout_nmr[:, 15]) .* u"°C"

# note tests seem to need to be in K rather than °C to work properly
# not all tests passing, some commented out, possibly due to different
# solvers and possibly to do with floating point error and issue with 
# the way the phase transition is being calculated
@testset "runmicro comparisons" begin
    @test_broken micro_out.relative_humidity[:, 1] ≈ rh1cm_nmr atol=0.2 # TODO make this work
    @test micro_out.relative_humidity[:, 2] ≈ rh2m_nmr atol=1e-5
    @test micro_out.wind_speed[:, 1] ≈ vel1cm_nmr atol=1e-6u"m/s"
    @test_broken micro_out.wind_speed[:, 2] ≈ vel2m_nmr atol=1e-6u"m/s" # now failing because of first day due to soil temps not being the same
    @test all(isapprox.(micro_out.wind_speed[:, 2], vel2m_nmr; atol=1e-6u"m/s"))
    @test u"K".(micro_out.air_temperature[:, 2]) ≈ ta2m_nmr atol=1e-5u"K"
    @test micro_out.sky_temperature ≈ u"K".(tskyC_nmr) atol=1u"K" # TODO make this better
    # The last column is different for these
    soiltemps_mat = reinterpret(reshape, typeof(1.0u"K"), micro_out.soil_temperature)'
    @test_broken all(isapprox.(soiltemps_mat, u"K".(Matrix(soiltemps_nmr)); atol=0.2u"K"))
    @test_broken all(isapprox.(soiltemps_mat[:, 2:10], u"K".(Matrix(soiltemps_nmr)[:, 2:10]); atol=0.5u"K")) # TODO make better!
end  
