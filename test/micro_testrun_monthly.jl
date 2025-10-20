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

terrain = Terrain(;
    elevation = microinput[:ALTT] * 1.0u"m", # elevation (m)
    horizon_angles = horizon_angles = (DataFrame(CSV.File("$testdir/data/init_monthly/hori.csv"))[:, 2]) * 1.0u"°",
    slope = microinput[:slope] * 1.0u"°",
    aspect = microinput[:azmuth] * 1.0u"°",
    roughness_height = microinput[:RUF] * 1.0u"m", # roughness height for standard mode TODO dispatch based on roughness pars
    karman_constant = 0.4, # Kármán constant
    dyer_constant = 16.0, # coefficient from Dyer and Hicks for Φ_m (momentum), γ
)

mineral_density = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][6]) * 1.0u"Mg/m^3" # soil minerals density (Mg/m3)
bulk_density = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][2]) * 1.0u"Mg/m^3" # dry soil bulk density (Mg/m3)

soil_thermal_model = CampbelldeVriesSoilThermal(;
    bulk_density, 
    mineral_density,
    deVries_shape_factor = 0.1, # de Vries shape factor, 0.33 for organic soils, 0.1 for mineral
    mineral_conductivity = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][4]) * 1.0u"W/m/K", # soil minerals thermal conductivity (W/mC)
    mineral_heat_capacity = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][5]) * 1.0u"J/kg/K", # soil minerals specific heat (J/kg-K)
    saturation_moisture = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][3]) * 1.0u"m^3/m^3", # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    recirculation_power = 4.0, # power for recirculation function
    return_flow_threshold = 0.162, # return-flow cutoff soil moisture, m^3/m^3
)

environment_daily = DailyTimeseries(;
    # daily environmental vectors
    albedo = (DataFrame(CSV.File("$testdir/data/init_monthly/REFLS.csv"))[days2do, 2] * 1.0), # substrate albedo (decimal %)
    shade = (DataFrame(CSV.File("$testdir/data/init_monthly/Minshades.csv"))[days2do, 2] * 1.0), # daily shade from vegetation (%)
    soil_wetness = (DataFrame(CSV.File("$testdir/data/init_monthly/PCTWET.csv"))[days2do, 2] * 1.0),
    surface_emissivity = (DataFrame(CSV.File("$testdir/data/init_monthly/SLES.csv"))[days2do, 2] * 1.0), # - surface emissivity
    cloud_emissivity = (DataFrame(CSV.File("$testdir/data/init_monthly/SLES.csv"))[days2do, 2] * 1.0), # - cloud emissivity
    rainfall = ((DataFrame(CSV.File("$testdir/data/init_monthly/rain.csv"))[days2do, 2] * 1.0) / 1000)u"kg/m^2", # monthly total rainfall
    deep_soil_temperature = (DataFrame(CSV.File("$testdir/data/init_monthly/tannulrun.csv"))[days2do, 2] * 1.0)u"°C", # daily deep soil temperatures
    leaf_area_index = fill(0.1, length(days)),
)

environment_minmax = MonthlyMinMaxEnvironment(;
    reference_temperature_min = (DataFrame(CSV.File("$testdir/data/init_monthly/TMINN.csv"))[days2do, 2] * 1.0)u"°C", # minimum air temperatures
    reference_temperature_max = (DataFrame(CSV.File("$testdir/data/init_monthly/TMAXX.csv"))[days2do, 2] * 1.0)u"°C", # maximum air temperatures
    reference_wind_min = (DataFrame(CSV.File("$testdir/data/init_monthly/WNMINN.csv"))[days2do, 2] * 1.0)u"m/s", # min wind speed (m/s)
    reference_wind_max = (DataFrame(CSV.File("$testdir/data/init_monthly/WNMAXX.csv"))[days2do, 2] * 1.0)u"m/s", # max wind speed (m/s)
    reference_humidity_min = (DataFrame(CSV.File("$testdir/data/init_monthly/RHMINN.csv"))[days2do, 2] * 1.0), # min relative humidity (%)
    reference_humidity_max = (DataFrame(CSV.File("$testdir/data/init_monthly/RHMAXX.csv"))[days2do, 2] * 1.0), # max relative humidity (%)
    cloud_min = (DataFrame(CSV.File("$testdir/data/init_monthly/CCMINN.csv"))[days2do, 2] * 1.0), # min cloud cover (%)
    cloud_max = (DataFrame(CSV.File("$testdir/data/init_monthly/CCMAXX.csv"))[days2do, 2] * 1.0), # max cloud cover (%)
    minima_times = [microinput[:TIMINS1], microinput[:TIMINS2], microinput[:TIMINS3], microinput[:TIMINS4]], # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
    maxima_times = [microinput[:TIMAXS1], microinput[:TIMAXS2], microinput[:TIMAXS3], microinput[:TIMAXS4]], # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
)

soil_moisture_model = example_soil_moisture_model(depths; bulk_density, mineral_density)
solar_model = SolarRadiation(; iuv = Bool(Int(microinput[:IUV])))

# now try the simulation function
problem = MicroProblem(;
    # locations, times, depths and heights 
    latitude = longlat[2]*1.0u"°",
    days = days[days2do], # days of year for solrad
    hours = collect(0.0:1:23.0), # hour of day for solrad
    depths,
    heights, # air nodes for temperature, wind speed and humidity profile
    # Objects defined above
    terrain,
    solar_model,
    soil_moisture_model,
    soil_thermal_model,
    environment_minmax,
    environment_daily,
    iterate_day = (microinput[:ndmax]), # number of iterations per day
    daily = Bool(Int(microinput[:microdaily])), # doing consecutive days?
    runmoist = Bool(Int(microinput[:runmoist])), # run soil moisture algorithm?
    hourly_rainfall = Bool(Int(microinput[:rainhourly])), # use hourly rainfall?
    spinup = Bool(Int(microinput[:spinup])), # spin-up the first day by iterate_day iterations?
    # intial conditions
    initial_soil_temperature = nothing, # initial soil temperature, # initial soil temperature
    initial_soil_moisture = (Array(DataFrame(CSV.File("$testdir/data/init_monthly/moists.csv"))[1:10, 2]) .* 1.0), # initial soil moisture
    #maximum_surface_temperature = u"K"(microinput[:maxsurf]u"°C")
)

# now try the simulation function
@time micro_out = Microclimate.solve(problem);

# subset NicheMapR predictions
vel1cm_nmr = collect(metout_nmr[:, 8]) .* 1u"m/s"
vel2m_nmr = collect(metout_nmr[:, 9]) .* 1u"m/s"
ta1cm_nmr = collect(metout_nmr[:, 4] .+ 273.15) .* 1u"K"
ta2m_nmr = collect(metout_nmr[:, 5] .+ 273.15) .* 1u"K"
rh1cm_nmr = collect(metout_nmr[:, 6])
rh2m_nmr = collect(metout_nmr[:, 7])
tskyC_nmr = collect(metout_nmr[:, 15]) .* u"°C"

air_temperature_matrix = hcat([p.air_temperature for p in micro_out.profile]...)'
humidity_matrix = hcat([p.relative_humidity for p in micro_out.profile]...)'
wind_matrix = hcat([p.wind_speed for p in micro_out.profile]...)'

@testset "runmicro comparisons" begin
    @test humidity_matrix[:, 1] ≈ rh1cm_nmr rtol=1e-1
    @test humidity_matrix[:, 2] ≈ rh2m_nmr rtol=1e-8
    @test wind_matrix[:, 1] ≈ vel1cm_nmr rtol=1e-2
    @test wind_matrix[:, 2] ≈ vel2m_nmr rtol=1e-8 
    @test u"K".(air_temperature_matrix[:, 1]) ≈ ta1cm_nmr rtol=1e-3
    @test u"K".(air_temperature_matrix[:, 2]) ≈ ta2m_nmr rtol=1e-8
    @test micro_out.sky_temperature ≈ u"K".(tskyC_nmr) rtol=1e-7
    @test all(isapprox.(micro_out.soil_temperature, u"K".(Matrix(soiltemps_nmr)); rtol=1e-2))
end  
