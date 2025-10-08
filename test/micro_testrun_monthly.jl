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
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349] # days of year to simulate - TODO leap years - why not use real dates?
LAIs = fill(0.1, length(days))
depths = ((DataFrame(CSV.File("$testdir/data/init_monthly/DEP.csv"))[:, 2]) / 100.0)u"m"
heights = [microinput[:Usrhyt], microinput[:Refhyt]]u"m" # air nodes for temperature, wind speed and humidity profile
soil_saturation_moisture = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][3]) * 1.0u"m^3/m^3" # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)

terrain = Terrain(;
    elevation = microinput[:ALTT] * 1.0u"m", # elevation (m)
    horizon_angles = (DataFrame(CSV.File("$testdir/data/init_monthly/hori.csv"))[:, 2]) * 1.0u"°", # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
    slope = microinput[:slope] * 1.0u"°",
    aspect = microinput[:azmuth] * 1.0u"°",
    roughness_height = microinput[:RUF] * 1.0u"m",
    karman_constant = 0.4, # Kármán constant
    dyer_constant = 16, # coefficient from Dyer and Hicks for Φ_m (momentum), γ
)

soil_thermal_model = CampbelldeVriesSoilThermal(;
    deVries_shape_factor = 0.1, # de Vries shape factor, 0.33 for organic soils, 0.1 for mineral
    mineral_conductivity = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][4]) * 1.0u"W/m/K", # soil minerals thermal conductivity (W/mC)
    mineral_density = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][6]) * 1.0u"Mg/m^3", # soil minerals density (Mg/m3)
    mineral_heat_capacity = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][5]) * 1.0u"J/kg/K", # soil minerals specific heat (J/kg-K)
    bulk_density = (CSV.File("$testdir/data/init_monthly/soilprop.csv")[1, 1][2]) * 1.0u"Mg/m^3", # dry soil bulk density (Mg/m3)
    recirculation_power = 4.0, # power for recirculation function
    return_flow_threshold = 0.162, # return-flow cutoff soil moisture, m^3/m^3
)

environment_minmax = MonthlyMinMaxEnvironment(;
    reference_temperature_min = (DataFrame(CSV.File("$testdir/data/init_monthly/TMINN.csv"))[:, 2] * 1.0)u"°C", # minimum air temperatures (°C)
    reference_temperature_max = (DataFrame(CSV.File("$testdir/data/init_monthly/TMAXX.csv"))[:, 2] * 1.0)u"°C", # maximum air temperatures (°C)
    reference_wind_min = (DataFrame(CSV.File("$testdir/data/init_monthly/WNMINN.csv"))[:, 2] * 1.0)u"m/s", # min wind speed (m/s)
    reference_wind_max = (DataFrame(CSV.File("$testdir/data/init_monthly/WNMAXX.csv"))[:, 2] * 1.0)u"m/s", # max wind speed (m/s)
    reference_humidity_min = (DataFrame(CSV.File("$testdir/data/init_monthly/RHMINN.csv"))[:, 2] * 1.0), # min relative humidity (%)
    reference_humidity_max = (DataFrame(CSV.File("$testdir/data/init_monthly/RHMAXX.csv"))[:, 2] * 1.0), # max relative humidity (%)
    cloud_min = (DataFrame(CSV.File("$testdir/data/init_monthly/CCMINN.csv"))[:, 2] * 1.0), # min cloud cover (%)
    cloud_max = (DataFrame(CSV.File("$testdir/data/init_monthly/CCMAXX.csv"))[:, 2] * 1.0), # max cloud cover (%)
    minima_times = [0, 0, 1, 1], # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
    maxima_times = [1, 1, 0, 0], # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
)

soil_moisture_model = SoilMoistureModel(; 
    # soil moisture model soil parameters
    air_entry_water_potential = (DataFrame(CSV.File("$testdir/data/init_monthly/PE.csv"))[:, 2] * 1.0u"J/kg"), # set up vector of ground emissivities for each day
    saturated_hydraulic_conductivity = (DataFrame(CSV.File("$testdir/data/init_monthly/KS.csv"))[:, 2] * 1.0u"kg*s/m^3"), # set up vector of ground emissivities for each day
    Campbells_b_parameter = (DataFrame(CSV.File("$testdir/data/init_monthly/BB.csv"))[:, 2] * 1.0), # set up vector of ground emissivities for each day
    soil_bulk_density2 = (DataFrame(CSV.File("$testdir/data/init_monthly/BD.csv"))[:, 2] * 1.0u"Mg/m^3"), # set up vector of ground emissivities for each day
    soil_mineral_density2 = (DataFrame(CSV.File("$testdir/data/init_monthly/DD.csv"))[:, 2] * 1.0u"Mg/m^3"), # set up vector of ground emissivities for each day
    # soil moisture plant parameters
    root_density = DataFrame(CSV.File("$testdir/data/init_monthly/L.csv"))[:, 2] * u"m/m^3", # root density at each node, mm/m3 (from Campell 1985 Soil Physics with Basic, p. 131) # max depth for water pooling on the surface, mm (to account for runoff)
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
    albedo = (DataFrame(CSV.File("$testdir/data/init_monthly/REFLS.csv"))[:, 2] * 1.0), # substrate albedo (decimal %)
    shade = (DataFrame(CSV.File("$testdir/data/init_monthly/Minshades.csv"))[:, 2] * 1.0), # daily shade (%)
    surface_emissivity = (DataFrame(CSV.File("$testdir/data/init_monthly/SLES.csv"))[:, 2] * 1.0), # - surface emissivity
    cloud_emissivity = fill(0.96, length(days)), # - surface emissivity
    rainfall = ((DataFrame(CSV.File("$testdir/data/init_monthly/rain.csv"))[:, 2] * 1.0))u"kg/m^2",
    deep_soil_temperature = (DataFrame(CSV.File("$testdir/data/init_monthly/tannulrun.csv"))[:, 2] * 1.0)u"°C", # daily deep soil temperatures
    leaf_area_index = (DataFrame(CSV.File("$testdir/data/init_monthly/LAI.csv"))[:, 2] * 1.0u"Mg/m^3"), # leaf area indices per day
)

environment_hourly = HourlyTimeseries(;
    reference_temperature = Float64.(CSV.File("$testdir/data/init_monthly/TAIRhr.csv").x)u"°C",
    reference_humidity = clamp.(Float64.(CSV.File("$testdir/data/init_monthly/RHhr.csv").x), 0, 100),
    reference_wind_speed = clamp.(Float64.(CSV.File("$testdir/data/init_monthly/WNhr.csv").x)u"m/s", 0.1u"m/s", (Inf)u"m/s"),
    solar_radiation = Float64.(CSV.File("$testdir/data/init_monthly/SOLRhr.csv").x)u"W/m^2",
    cloud_cover = clamp.(Float64.(CSV.File("$testdir/data/init_monthly/CLDhr.csv").x), 0, 100),
    rainfall = clamp.(Float64.(CSV.File("$testdir/data/init_monthly/RAINhr.csv").x)u"kg" / u"m^2", 0u"kg/m^2", Inf * u"kg/m^2"),
    zenith_angle=nothing,
    longwave_radiation=nothing,
)

solar_model = SolarRadiation(; iuv = Bool(Int(microinput[:IUV])))

# now try the simulation function
problem = MicroProblem(;
    # locations, times, depths and heights
    latitude = (microinput[:ALAT] + microinput[:AMINUT] / 60) * 1.0u"°", # latitude
    days, # days of year to simulate - TODO leap years
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
    hourly_rainfall = Bool(Int(microinput[:rainhourly])), # is hourly rainfall to be used?
    spinup = Bool(Int(microinput[:spinup])), # spin-up the first day by iterate_day iterations?
    # intial conditions
    initial_soil_temperature = u"K".((DataFrame(CSV.File("$testdir/data/init_monthly/soilinit.csv"))[1:length(depths), 2] * 1.0)u"°C"), # initial soil temperature
    initial_soil_moisture = (Array(DataFrame(CSV.File("$testdir/data/init_monthly/moists.csv"))[1:10, 2]) .* 1.0), # initial soil moisture
    #maximum_surface_temperature = u"K"(microinput[:maxsurf]u"°C")
)

solrad_out = solve_solar(problem)
ndays = length(days)
hours = 0:1:23
nhours = length(hours)
nsteps = ndays * nhours
numnodes_a = length(depths) # number of soil nodes for temperature calcs and final output
output = MicroResult(nsteps, numnodes_a)
    reference_temperature, reference_wind_speed, reference_humidity, cloud_cover = 
        hourly_vars(environment_minmax, solrad_out, environment_daily)
interpolate_minmax!(output, environment_minmax, environment_daily, environment_hourly, solrad_out)

# TODO allow vector of pre-calculated soil moisture to be provided as input
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

# note tests seem to need to be in K rather than °C to work properly
# not all tests passing, some commented out, possibly due to different
# solvers and possibly to do with floating point error and issue with 
# the way the phase transition is being calculated
@testset "runmicro comparisons" begin
    @test_broken micro_out.relative_humidity[:, 1] ≈ rh1cm_nmr atol=0.2 # TODO make this work
    @test micro_out.relative_humidity[:, 2] ≈ rh2m_nmr atol=1e-5
    @test micro_out.wind_speed[:, 1] ≈ vel1cm_nmr atol=2e-1u"m/s" # now failing because of first day due to soil temps not being the same
    @test micro_out.wind_speed[:, 2] ≈ vel2m_nmr atol=1e-6u"m/s" 
    @test u"K".(micro_out.air_temperature[:, 2]) ≈ ta2m_nmr atol=1e-5u"K"
    @test micro_out.sky_temperature ≈ u"K".(tskyC_nmr) atol=1u"K" # TODO make this better
    # The last column is different for these
    soiltemps_mat = reinterpret(reshape, typeof(1.0u"K"), micro_out.soil_temperature)'
    @test_broken all(isapprox.(soiltemps_mat, u"K".(Matrix(soiltemps_nmr)); atol=0.2u"K"))
    @test_broken all(isapprox.(soiltemps_mat[:, 2:10], u"K".(Matrix(soiltemps_nmr)[:, 2:10]); atol=0.5u"K")) # TODO make better!
end  
