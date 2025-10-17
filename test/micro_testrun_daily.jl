using Microclimate
using Unitful
using CSV, DataFrames
using Test

testdir = realpath(joinpath(dirname(pathof(Microclimate)), "../test"))

# read in output from NicheMapR
soil_temperature_nmr = (DataFrame(CSV.File("$testdir/data/soil_FordDryLake.csv"))[:, 5:14]) .* u"°C"
soil_moisture_nmr = (DataFrame(CSV.File("$testdir/data/soilmoist_FordDryLake.csv"))[:, 5:14])
soil_conductivity_nmr = (DataFrame(CSV.File("$testdir/data/tcond_FordDryLake.csv"))[:, 5:14])
metout_nmr = (DataFrame(CSV.File("$testdir/data/metout_FordDryLake.csv"))[:, 2:21])

microinput_vec = DataFrame(CSV.File("$testdir/data/init_daily/microinput.csv"))[:, 2]

names = [
    :doynum, :RUF, :ERR, :Usrhyt, :Refhyt, :Numtyps, :Z01, :kZ02, :ZH1, :ZH2,
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

days = collect(1:Int(length(soil_temperature_nmr[:, 1]) / 24)) # days of year to run (for solrad)
depths = ((DataFrame(CSV.File("$testdir/data/init_daily/DEP.csv"))[:, 2]) / 100.0)u"m" # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
heights = [microinput[:Usrhyt], microinput[:Refhyt]]u"m" # air nodes for temperature, wind speed and humidity profile
days2do = 30
hours2do = days2do * 24

terrain = Terrain(;
    elevation = microinput[:ALTT] * 1.0u"m", # elevation (m)
    horizon_angles = (DataFrame(CSV.File("$testdir/data/init_daily/hori.csv"))[:, 2]) * 1.0u"°", # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
    slope = microinput[:slope] * 1.0u"°",
    aspect = microinput[:azmuth] * 1.0u"°",
    roughness_height = microinput[:RUF] * 1.0u"m", # roughness height for standard mode TODO dispatch based on roughness pars
    zh = microinput[:ZH] * 1.0u"m", # heat transfer roughness height for Campbell and Norman mode
    d0 = microinput[:D0] * 1.0u"m", # zero plane displacement correction factor
    karman_constant = 0.4, # Kármán constant
)

soil_thermal_model = CampbelldeVriesSoilThermal(;
    deVries_shape_factor = 0.1, # de Vries shape factor, 0.33 for organic soils, 0.1 for mineral
    mineral_conductivity = (CSV.File("$testdir/data/init_daily/soilprop.csv")[1, 1][4]) * 1.0u"W/m/K", # soil minerals thermal conductivity (W/mC)
    mineral_density = (CSV.File("$testdir/data/init_daily/soilprop.csv")[1, 1][6]) * 1.0u"Mg/m^3", # soil minerals density (Mg/m3)
    mineral_heat_capacity = (CSV.File("$testdir/data/init_daily/soilprop.csv")[1, 1][5]) * 1.0u"J/kg/K", # soil minerals specific heat (J/kg-K)
    bulk_density = (CSV.File("$testdir/data/init_daily/soilprop.csv")[1, 1][2]) * 1.0u"Mg/m^3", # dry soil bulk density (Mg/m3)
    saturation_moisture = (CSV.File("$testdir/data/init_daily/soilprop.csv")[1, 1][3]) * 1.0u"m^3/m^3", # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    recirculation_power = 4.0, # power for recirculation function
    return_flow_threshold = 0.162, # return-flow cutoff soil moisture, m^3/m^3
)

environment_minmax = nothing
# MonthlyMinMaxEnvironment(;
#     air_temperature_min = (DataFrame(CSV.File("$testdir/data/init_daily/TMINN.csv"))[1:days2do, 2] * 1.0)u"°C", # minimum air temperatures (°C)
#     air_temperature_max = (DataFrame(CSV.File("$testdir/data/init_daily/TMAXX.csv"))[1:days2do, 2] * 1.0)u"°C", # maximum air temperatures (°C)
#     wind_min = nothing,
#     wind_max = nothing,
#     humidity_min = nothing,
#     humidity_max = nothing,
#     cloud_min = nothing,
#     cloud_max = nothing,
#     minima_times = [0, 0, 1, 1], # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
#     maxima_times = [1, 1, 0, 0], # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
# )

soil_moisture_model = SoilMoistureModel(; 
    # soil moisture model soil parameters
    air_entry_water_potential = (DataFrame(CSV.File("$testdir/data/init_daily/PE.csv"))[:, 2] * 1.0u"J/kg"), # set up vector of ground emissivities for each day
    saturated_hydraulic_conductivity = (DataFrame(CSV.File("$testdir/data/init_daily/KS.csv"))[:, 2] * 1.0u"kg*s/m^3"), # set up vector of ground emissivities for each day
    Campbells_b_parameter = (DataFrame(CSV.File("$testdir/data/init_daily/BB.csv"))[:, 2] * 1.0), # set up vector of ground emissivities for each day
    soil_bulk_density2 = (DataFrame(CSV.File("$testdir/data/init_daily/BD.csv"))[:, 2] * 1.0u"Mg/m^3"), # set up vector of ground emissivities for each day
    soil_mineral_density2 = (DataFrame(CSV.File("$testdir/data/init_daily/DD.csv"))[:, 2] * 1.0u"Mg/m^3"), # set up vector of ground emissivities for each day
    # soil moisture plant parameters
    root_density = DataFrame(CSV.File("$testdir/data/init_daily/L.csv"))[:, 2] * u"m/m^3", # root density at each node, mm/m3 (from Campell 1985 Soil Physics with Basic, p. 131) # max depth for water pooling on the surface, mm (to account for runoff)
    root_resistance = microinput[:RW] * u"m^3/kg/s", # resistance per unit length of root, m3 kg-1 s-1
    stomatal_closure_potential = -microinput[:PC] * u"J/kg", # critical leaf water potential for stomatal closure, J kg-1
    leaf_resistance = microinput[:RL] * u"m^4/kg/s", # resistance per unit length of leaf, m3 kg-1 s-1
    stomatal_stability_parameter = microinput[:SP], # stability parameter, -
    root_radius = microinput[:R1]u"m", # root radius, m    
    # soil moisture simulation controls
    moist_error = microinput[:IM]u"kg/m^2/s", # maximum overall mass balance error allowed, kg
    moist_count = microinput[:MAXCOUNT],
    moist_step = microinput[:moiststep]u"s",    
    maxpool = microinput[:maxpool] * 1000.0u"kg/m^2", # max depth for water pooling on the surface, mm (to account for runoff)
)

environment_daily = DailyTimeseries(;
    # daily environmental vectors
    albedo = (DataFrame(CSV.File("$testdir/data/init_daily/REFLS.csv"))[1:days2do, 2] * 1.0), # substrate albedo (decimal %)
    shade = (DataFrame(CSV.File("$testdir/data/init_daily/Minshades.csv"))[1:days2do, 2] * 1.0), # daily shade (%)
    soil_wetness = (DataFrame(CSV.File("$testdir/data/init_daily/PCTWET.csv"))[1:days2do, 2] * 1.0),
    surface_emissivity = (DataFrame(CSV.File("$testdir/data/init_daily/SLES.csv"))[1:days2do, 2] * 1.0), # - surface emissivity
    cloud_emissivity = (DataFrame(CSV.File("$testdir/data/init_daily/SLES.csv"))[1:days2do, 2] * 1.0), # - cloud emissivity
    rainfall = ((DataFrame(CSV.File("$testdir/data/init_daily/rain.csv"))[1:days2do, 2] * 1.0))u"kg/m^2",
    deep_soil_temperature = (DataFrame(CSV.File("$testdir/data/init_daily/tannulrun.csv"))[1:days2do, 2] * 1.0)u"°C", # daily deep soil temperatures
    leaf_area_index = (DataFrame(CSV.File("$testdir/data/init_daily/LAI.csv"))[:, 2] * 1.0u"Mg/m^3"), # leaf area indices per day
)

environment_hourly = HourlyTimeseries(;
    reference_temperature = Float64.(CSV.File("$testdir/data/init_daily/TAIRhr.csv").x[1:hours2do])u"°C",
    reference_humidity = clamp.(Float64.(CSV.File("$testdir/data/init_daily/RHhr.csv").x[1:hours2do]), 0, 100),
    reference_wind_speed = clamp.(Float64.(CSV.File("$testdir/data/init_daily/WNhr.csv").x[1:hours2do])u"m/s", 0.1u"m/s", (Inf)u"m/s"),
    solar_radiation = Float64.(CSV.File("$testdir/data/init_daily/SOLRhr.csv").x[1:hours2do])u"W/m^2",
    cloud_cover = clamp.(Float64.(CSV.File("$testdir/data/init_daily/CLDhr.csv").x[1:hours2do]), 0, 100),
    rainfall = clamp.(Float64.(CSV.File("$testdir/data/init_daily/RAINhr.csv").x[1:hours2do])u"kg/m^2", 0u"kg/m^2", Inf * u"kg/m^2"),
    zenith_angle=nothing,
    longwave_radiation=nothing,
)

solar_model = SolarRadiation(; iuv = Bool(Int(microinput[:IUV])))

# now try the simulation function
problem = MicroProblem(;
    # locations, times, depths and heights
    latitude = (microinput[:ALAT] + microinput[:AMINUT] / 60) * 1.0u"°", # latitude
    days = days[1:days2do], # days of year to simulate - TODO leap years
    hours = 0:1:23, # hour of day for solrad # TODO how and in what context would users change this
    depths, # soil nodes - keep spacing close near the surface
    heights, # air nodes for temperature, wind speed and humidity profile
    # Objects defined above
    terrain,
    solar_model,
    soil_moisture_model,
    soil_thermal_model,
    environment_minmax,
    environment_daily,
    environment_hourly,
    iterate_day = (microinput[:ndmax]), # number of iterations per day
    daily = Bool(Int(microinput[:microdaily])), # doing consecutive days?
    runmoist = Bool(Int(microinput[:runmoist])), # run soil moisture algorithm?
    hourly_rainfall = Bool(Int(microinput[:rainhourly])), # use hourly rainfall?
    spinup = Bool(Int(microinput[:spinup])), # spin-up the first day by iterate_day iterations?
    # intial conditions
    initial_soil_temperature = u"K".((DataFrame(CSV.File("$testdir/data/init_daily/soilinit.csv"))[1:length(depths), 2] * 1.0)u"°C"), # initial soil temperature
    initial_soil_moisture = (Array(DataFrame(CSV.File("$testdir/data/init_daily/moists.csv"))[1:10, 2]) .* 1.0), # initial soil moisture
    #maximum_surface_temperature = u"K"(microinput[:maxsurf]u"°C")
)

# now try the simulation function
@time micro_out = Microclimate.solve(problem);

# TODO test plotting again at some stage, but it slows down CI a lot
# plot(micro_out)

@testset "runmicro comparisons" begin
    @test all(isapprox.(micro_out.soil_temperature[:, 1:10], u"K".(Matrix(soil_temperature_nmr[1:hours2do, 1:10])); rtol=1e-2)) # TODO make better!
    @test all(isapprox.(micro_out.soil_moisture[:, 1:10], Matrix(soil_moisture_nmr[1:hours2do, 1:10]); rtol=1e1)) # TODO make better!
    @test all(isapprox.(micro_out.soil_thermal_conductivity[:, 1:10], Matrix(soil_conductivity_nmr[1:hours2do, 1:10])u"W * m^-1 * K^-1"; rtol=1e1)) # TODO make better!
end 
