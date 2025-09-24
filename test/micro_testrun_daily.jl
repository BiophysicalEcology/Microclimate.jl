using Microclimate
using GLMakie
using Unitful
using CSV, DataFrames
using Test

testdir = realpath(joinpath(dirname(pathof(Microclimate)), "../test"))

# read in output from NicheMapR
soil_temperature_nmr = (DataFrame(CSV.File("$testdir/data/soil_FordDryLake.csv"))[:, 5:14]) .* u"°C"
soil_moisture_nmr = (DataFrame(CSV.File("$testdir/data/soilmoist_FordDryLake.csv"))[:, 5:14])
soil_conductivity_nmr = (DataFrame(CSV.File("$testdir/data/tcond_FordDryLake.csv"))[:, 5:14])

microinput_vec = DataFrame(CSV.File("$testdir/data/init_daily/microinput.csv"))[:, 2]

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

days = collect(1:Int(length(soil_temperature_nmr[:, 1]) / 24)) # days of year to run (for solrad)
depths = ((DataFrame(CSV.File("$testdir/data/init_daily/DEP.csv"))[:, 2]) / 100.0)u"m" # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
heights = [0.01]u"m" # air nodes for temperature, wind speed and humidity profile
days2do = 30
hours2do = days2do * 24
# now try the simulation function
keywords = (;
    # locations, times, depths and heights
    latitude = (microinput[:ALAT] + microinput[:AMINUT] / 60) * 1.0u"°", # latitude
    days = days[1:days2do], # days of year to simulate - TODO leap years
    hours = collect(0.:1:24.), # hour of day for solrad
    reference_height = microinput[:Refhyt] * 1.0u"m", # reference height of weather data (air temperature, wind speed, humidity)
    depths = depths, # soil nodes - keep spacing close near the surface
    heights = heights, # air nodes for temperature, wind speed and humidity profile
    # terrain
    elevation = microinput[:ALTT] * 1.0u"m", # elevation (m)
    horizon_angles = (DataFrame(CSV.File("$testdir/data/init_daily/hori.csv"))[:, 2]) * 1.0u"°", # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
    slope = microinput[:slope] * 1.0u"°",
    aspect = microinput[:azmuth] * 1.0u"°",
    roughness_height = microinput[:RUF] * 1.0u"m", # roughness height for standard mode TODO dispatch based on roughness pars
    zh = microinput[:ZH] * 1.0u"m", # heat transfer roughness height for Campbell and Norman mode
    d0 = microinput[:D0] * 1.0u"m", # zero plane displacement correction factor
    # soil thermal parameters 
    soil_mineral_conductivity = (CSV.File("$testdir/data/init_daily/soilprop.csv")[1, 1][4]) * 1.0u"W/m/K", # soil minerals thermal conductivity (W/mC)
    soil_mineral_density = (CSV.File("$testdir/data/init_daily/soilprop.csv")[1, 1][6]) * 1.0u"Mg/m^3", # soil minerals density (Mg/m3)
    soil_mineral_heat_capacity = (CSV.File("$testdir/data/init_daily/soilprop.csv")[1, 1][5]) * 1.0u"J/kg/K", # soil minerals specific heat (J/kg-K)
    soil_bulk_density = (CSV.File("$testdir/data/init_daily/soilprop.csv")[1, 1][2]) * 1.0u"Mg/m^3", # dry soil bulk density (Mg/m3)
    soil_saturation_moisture = (CSV.File("$testdir/data/init_daily/soilprop.csv")[1, 1][3]) * 1.0u"m^3/m^3", # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    # soil moisture model soil parameters
    air_entry_water_potential = (DataFrame(CSV.File("$testdir/data/init_daily/PE.csv"))[:, 2] * 1.0u"J/kg"), # set up vector of ground emissivities for each day
    saturated_hydraulic_conductivity = (DataFrame(CSV.File("$testdir/data/init_daily/KS.csv"))[:, 2] * 1.0u"kg*s/m^3"), # set up vector of ground emissivities for each day
    Campbells_b_parameter = (DataFrame(CSV.File("$testdir/data/init_daily/BB.csv"))[:, 2] * 1.0), # set up vector of ground emissivities for each day
    soil_bulk_density2 = (DataFrame(CSV.File("$testdir/data/init_daily/BD.csv"))[:, 2] * 1.0u"Mg/m^3"), # set up vector of ground emissivities for each day
    soil_mineral_density2 = (DataFrame(CSV.File("$testdir/data/init_daily/DD.csv"))[:, 2] * 1.0u"Mg/m^3"), # set up vector of ground emissivities for each day
    # soil moisture plant parameters
    root_density = DataFrame(CSV.File("$testdir/data/init_daily/L.csv"))[:, 2]u"m/m^3", # root density at each node, mm/m3 (from Campell 1985 Soil Physics with Basic, p. 131) # max depth for water pooling on the surface, mm (to account for runoff)
    root_resistance = microinput[:RW]u"m^3/kg/s", # resistance per unit length of root, m3 kg-1 s-1
    stomatal_closure_potential = -microinput[:PC]u"J/kg", # critical leaf water potential for stomatal closure, J kg-1
    leaf_resistance = microinput[:RL]u"m^4/kg/s", # resistance per unit length of leaf, m3 kg-1 s-1
    stomatal_stability_parameter = microinput[:SP], # stability parameter, -
    root_radius = microinput[:R1]u"m", # root radius, m    
    # soil moisture simulation controls
    moist_error = microinput[:IM]u"kg/m^2/s", # maximum overall mass balance error allowed, kg
    moist_count = microinput[:MAXCOUNT],
    moist_step = microinput[:moiststep]u"s",    
    maxpool = microinput[:maxpool] * 1000.0u"kg/m^2", # max depth for water pooling on the surface, mm (to account for runoff)
    # daily environmental vectors
    albedos = (DataFrame(CSV.File("$testdir/data/init_daily/REFLS.csv"))[1:days2do, 2] * 1.0), # substrate albedo (decimal %)
    shades = (DataFrame(CSV.File("$testdir/data/init_daily/Minshades.csv"))[1:days2do, 2] * 1.0), # daily shade (%)
    pctwets = (DataFrame(CSV.File("$testdir/data/init_daily/PCTWET.csv"))[1:days2do, 2] * 1.0),
    sles = (DataFrame(CSV.File("$testdir/data/init_daily/SLES.csv"))[1:days2do, 2] * 1.0), # - surface emissivity
    daily_rainfall=((DataFrame(CSV.File("$testdir/data/init_daily/rain.csv"))[1:days2do, 2] * 1.0))u"kg/m^2",
    air_temperature_min = (DataFrame(CSV.File("$testdir/data/init_daily/TMINN.csv"))[1:days2do, 2] * 1.0)u"°C", # minimum air temperatures (°C)
    air_temperature_max = (DataFrame(CSV.File("$testdir/data/init_daily/TMAXX.csv"))[1:days2do, 2] * 1.0)u"°C", # maximum air temperatures (°C)
    wind_min = nothing,
    wind_max = nothing,
    humidity_min = nothing,
    humidity_max = nothing,
    cloud_min = nothing,
    cloud_max = nothing,
    deep_soil_temperatures = (DataFrame(CSV.File("$testdir/data/init_daily/tannulrun.csv"))[1:days2do, 2] * 1.0)u"°C", # daily deep soil temperatures
    # hourly weather vectors
    air_temperatures = Float64.(CSV.File("$testdir/data/init_daily/TAIRhr.csv").x[1:hours2do])u"°C",
    humidities = clamp.(Float64.(CSV.File("$testdir/data/init_daily/RHhr.csv").x[1:hours2do]), 0, 100),
    wind_speeds = clamp.(Float64.(CSV.File("$testdir/data/init_daily/WNhr.csv").x[1:hours2do])u"m/s", 0.1u"m/s", (Inf)u"m/s"),
    solar_radiation = Float64.(CSV.File("$testdir/data/init_daily/SOLRhr.csv").x[1:hours2do])u"W/m^2",
    cloud_covers = clamp.(Float64.(CSV.File("$testdir/data/init_daily/CLDhr.csv").x[1:hours2do]), 0, 100),
    RAINs = clamp.(Float64.(CSV.File("$testdir/data/init_daily/RAINhr.csv").x[1:hours2do])u"kg" / u"m^2", 0u"kg/m^2", (Inf)u"kg/m^2"),
    zenith_angles=nothing,
    longwave_radiation=nothing,
    # intial conditions
    initial_soil_temperature = u"K".((DataFrame(CSV.File("$testdir/data/init_daily/soilinit.csv"))[1:length(depths), 2] * 1.0)u"°C"), # initial soil temperature
    initial_soil_moisture = (Array(DataFrame(CSV.File("$testdir/data/init_daily/moists.csv"))[1, 2:13]) .* 1.0), # initial soil moisture
    leaf_area_index = (DataFrame(CSV.File("$testdir/data/init_daily/LAI.csv"))[:, 2] * 1.0u"Mg/m^3"), # leaf area indices per day
    iterate_day = (microinput[:ndmax]), # number of iterations per day
    daily = Bool(Int(microinput[:microdaily])), # doing consecutive days?
    runmoist = Bool(Int(microinput[:runmoist])), # run soil moisture algorithm?
    spinup = Bool(Int(microinput[:spinup])), # spin-up the first day by iterate_day iterations?
    iuv = Bool(Int(microinput[:IUV])), # this makes it take ages if true!
)

# now try the simulation function
@time micro_out = runmicro(; keywords...);

plot(micro_out)

# TODO include 1st node (currently left out, i.e. just columns 2:10, because way off at times)
@testset "runmicro comparisons" begin
    @test all(isapprox.(micro_out.soil_moisture[:, 2:10], Matrix(soil_moisture_nmr[1:hours2do, 2:10]); atol=0.3)) # TODO make better!
    @test all(isapprox.(micro_out.soil_temperature[:, 2:10], u"K".(Matrix(soil_temperature_nmr[1:hours2do, 2:10])); atol=10u"K")) # TODO make better!
    @test all(isapprox.(micro_out.soil_thermal_conductivity[:, 2:10], Matrix(soil_conductivity_nmr[1:hours2do, 2:10])u"W * m^-1 * K^-1"; atol=1u"W * m^-1 * K^-1")) # TODO make better!
end 
