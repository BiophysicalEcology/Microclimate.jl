using Microclimate
using SolarRadiation
using FluidProperties
using Unitful
using CSV, DataFrames
using Test

testdir = realpath(joinpath(dirname(pathof(Microclimate)), "../test"))

# read in output from NicheMapR and input variables
soiltemps_nmr = (DataFrame(CSV.File("$testdir/data/soil_monthly_snow.csv"))[:, 4:13]) .* u"°C"
metout_nmr = DataFrame(CSV.File("$testdir/data/metout_monthly_snow.csv"; missingstring="NA"))
microinput_vec = DataFrame(CSV.File("$testdir/data/init_monthly_snow/microinput.csv"))[:, 2]

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

longlat = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/longlat.csv"))[:, 2] * 1.0)
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
LAIs = fill(0.1, length(days))
depths = ((DataFrame(CSV.File("$testdir/data/init_monthly_snow/DEP.csv"))[:, 2]) / 100.0)u"m"
heights = [microinput[:Usrhyt], microinput[:Refhyt]]u"m" # air nodes for temperature, wind speed and humidity profile
days2do = 1:12

precomputed_soil_moisture = (Array(DataFrame(CSV.File("$testdir/data/init_monthly/moists.csv"))[:, 2:13]) .* 1.0)

site = Site(;
    latitude = longlat[2] * 1.0u"°",
    longitude = longlat[1] * 1.0u"°",
    elevation = microinput[:ALTT] * 1.0u"m",
    slope = microinput[:slope] * 1.0u"°",
    aspect = microinput[:azmuth] * 1.0u"°",
    horizon_angles = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/hori.csv"))[:, 2]) * 1.0u"°",
    sky_view_fraction = 1.0,
    albedo = DataFrame(CSV.File("$testdir/data/init_monthly_snow/REFLS.csv"))[1, 2] * 1.0,
    roughness_height = microinput[:RUF] * 1.0u"m",
    atmospheric_pressure = atmospheric_pressure(microinput[:ALTT] * 1.0u"m"),
)
boundary_layer_model = MoninObukhov(; karman_constant=0.4, dyer_constant=16.0)

mineral_density = (CSV.File("$testdir/data/init_monthly_snow/soilprop.csv")[1, 1][6]) * 1.0u"Mg/m^3" # soil minerals density (Mg/m3)
bulk_density = (CSV.File("$testdir/data/init_monthly_snow/soilprop.csv")[1, 1][2]) * 1.0u"Mg/m^3" # dry soil bulk density (Mg/m3)

soil_properties_model = CampbelldeVriesSoilProperties(;
    de_vries_shape_factor = 0.1, # de Vries shape factor, 0.33 for organic soils, 0.1 for mineral
    mineral_conductivity = (CSV.File("$testdir/data/init_monthly_snow/soilprop.csv")[1, 1][4]) * 1.0u"W/m/K", # soil minerals thermal conductivity (W/mC)
    mineral_heat_capacity = (CSV.File("$testdir/data/init_monthly_snow/soilprop.csv")[1, 1][5]) * 1.0u"J/kg/K", # soil minerals specific heat (J/kg-K)
    recirculation_power = 4.0, # power for recirculation function
    return_flow_threshold = 0.162, # return-flow cutoff soil moisture, m^3/m^3
)

rainfall = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/rain.csv"))[days2do, 2] * 1.0)u"kg/m^2"

environment_hourly = HourlyTimeseries(;
    pressure = fill(atmospheric_pressure((microinput[:ALTT])*1.0u"m"), length(days2do)*24),
    reference_temperature = nothing,
    reference_humidity = nothing,
    reference_wind_speed = nothing,
    global_radiation = nothing,
    cloud_cover = nothing,
    rainfall = repeat(rainfall ./ 24, inner=24),
    zenith_angle = nothing,
    longwave_radiation = nothing,
)

environment_daily = DailyTimeseries(;
    # daily environmental vectors
    shade = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/Minshades.csv"))[days2do, 2] * 1.0) ./ 100.0, # daily shade from vegetation (fractional)
    soil_wetness = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/PCTWET.csv"))[days2do, 2] * 1.0) ./ 100.0, # daily soil wetness (fractional)
    surface_emissivity = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/SLES.csv"))[days2do, 2] * 1.0), # - surface emissivity
    cloud_emissivity = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/SLES.csv"))[days2do, 2] * 1.0), # - cloud emissivity
    rainfall, # monthly total rainfall (mm = kg/m²)
    deep_soil_temperature = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/tannulrun.csv"))[days2do, 2] * 1.0)u"°C", # daily deep soil temperatures
    leaf_area_index = fill(0.1, length(days)),
)

environment_minmax = MonthlyMinMaxEnvironment(;
    reference_temperature_min = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/TMINN.csv"))[days2do, 2] * 1.0)u"°C", # minimum air temperatures
    reference_temperature_max = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/TMAXX.csv"))[days2do, 2] * 1.0)u"°C", # maximum air temperatures
    reference_wind_min = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/WNMINN.csv"))[days2do, 2] * 1.0)u"m/s", # min wind speed (m/s)
    reference_wind_max = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/WNMAXX.csv"))[days2do, 2] * 1.0)u"m/s", # max wind speed (m/s)
    reference_humidity_min = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/RHMINN.csv"))[days2do, 2] * 1.0) ./ 100.0, # min relative humidity (fractional)
    reference_humidity_max = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/RHMAXX.csv"))[days2do, 2] * 1.0) ./ 100.0, # max relative humidity (fractional)
    cloud_min = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/CCMINN.csv"))[days2do, 2] * 1.0) ./ 100.0, # min cloud cover (fractional)
    cloud_max = (DataFrame(CSV.File("$testdir/data/init_monthly_snow/CCMAXX.csv"))[days2do, 2] * 1.0) ./ 100.0, # max cloud cover (fractional)
    minima_times = [microinput[:TIMINS1], microinput[:TIMINS2], microinput[:TIMINS3], microinput[:TIMINS4]], # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
    maxima_times = [microinput[:TIMAXS1], microinput[:TIMAXS2], microinput[:TIMAXS3], microinput[:TIMAXS4]], # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
)

_runmoist = Bool(Int(microinput[:runmoist]))
soil_profile = SoilProfile(;
    bulk_density = fill(bulk_density, length(depths)),
    mineral_density = fill(mineral_density, length(depths)),
)
soil_hydraulic_model = example_soil_hydraulic_model(depths;
    root_density = fill(0.0, length(depths))u"m/m^3")
radiation = RadiationModel(;
    solar_radiation_model = SolarProblem(;
        diffuse_model = Bool(Int(microinput[:IUV])) ? SolarRadiation.ChandrasekharScattering() : SolarRadiation.NoScattering(),
    ),
)

# Set up time mode from the daily/spinup flags
_daily = Bool(Int(microinput[:microdaily]))
_spinup = Bool(Int(microinput[:spinup]))
time_mode = _daily ? ConsecutiveDayMode(; spinup_first_day=_spinup) :
    NonConsecutiveDayMode(; iterations_per_day=Int(microinput[:ndmax]))

snow_model = SnowModel(;
    snow_temperature_threshold = microinput[:snowtemp] * u"°C",
    snow_density = microinput[:snowdens] * u"g/cm^3",
    snow_melt_factor = microinput[:snowmelt],
    undercatch = microinput[:undercatch],
    rain_multiplier = microinput[:rainmult],
    rain_melt_factor = microinput[:rainmelt],
    density_function = (microinput[:densfun1], microinput[:densfun2], microinput[:densfun3], microinput[:densfun4]),
    snow_conductivity = microinput[:snowcond] * u"W/m/K",
    canopy_interception = microinput[:intercept],
)

config = MicroConfig(;
    time_mode,
    convergence = FixedSoilTemperatureIterations(Int(microinput[:ndmax])),
    rainfall_schedule = Bool(Int(microinput[:rainhourly])) ? HourlyRainfall() : DailyRainfall(),
    soil_moisture_strategy = _runmoist ? DynamicSoilMoisture() :
        PrescribedSoilMoisture(; precomputed_soil_moisture),
)

model = MicroModel(;
    days = days[days2do],
    hours = collect(0.0:1:23.0),
    depths,
    heights,
    radiation,
    soil_profile,
    soil_properties_model,
    soil_hydraulic_model,
    snow_model,
    boundary_layer_model,
    config,
)
inputs = MicroInputs(;
    site,
    environment_minmax,
    environment_daily,
    environment_hourly,
    initial_soil_temperature = nothing,
    initial_soil_moisture = precomputed_soil_moisture[1:10, 1],
)
problem = MicroProblem(model, inputs)

@time micro_out = Microclimate.solve(problem);

# subset NicheMapR predictions
vel2m_nmr = collect(metout_nmr[:, 9]) .* 1u"m/s"
rh2m_nmr = collect(metout_nmr[:, 7]) ./ 100.0
tskyC_nmr = collect(metout_nmr[:, 15]) .* u"°C"
solr_nmr = collect(metout_nmr[:, 14]) .* u"W/m^2"
# Snow columns may contain missing values (NA in R output)
snowfall_nmr = metout_nmr[:, 18] .* 1u"cm/hr"
snowdepth_nmr = metout_nmr[:, 19] .* 1u"cm"
snowdensity_nmr = metout_nmr[:, 20] .* 1u"g/cm^3"

air_temperature_matrix = micro_out.profile.air_temperature
humidity_matrix = micro_out.profile.relative_humidity
wind_matrix = micro_out.profile.wind_speed
snow_depth_matrix = micro_out.snow_depth
snow_density_matrix = micro_out.snow_density

# Find rows where NicheMapR has non-missing snow data (first hour of each day)
snow_valid = .!ismissing.(snowdepth_nmr)

@testset "runmicro comparisons" begin
    @test humidity_matrix[:, 2] ≈ rh2m_nmr rtol=1e-8
    @test wind_matrix[:, 2] ≈ vel2m_nmr rtol=1e-8
    @test micro_out.sky_temperature ≈ u"K".(tskyC_nmr) rtol=1e-7
    @test micro_out.global_radiation ≈ solr_nmr rtol=1e-4
    # Snow outputs: tolerances follow PR #102. Snow physics across non-consecutive
    # days is approximate (precipitation distribution, density carry-over, albedo
    # — see PR description). Depth/density compared only over the first 3 days
    # where Julia ↔ Fortran tracks essentially exactly.
    @test ustrip.(u"cm/hr", micro_out.snow_fall[snow_valid]) ≈ ustrip.(u"cm/hr", snowfall_nmr[snow_valid]) rtol=1e-3
    @test ustrip.(u"cm", micro_out.snow_depth[1:72]) ≈ ustrip.(u"cm", snowdepth_nmr[1:72]) rtol=1e-2
    @test ustrip.(u"g/cm^3", micro_out.snow_density[1:72]) ≈ ustrip.(u"g/cm^3", snowdensity_nmr[1:72]) rtol=1e-2
end

# Visual comparisons — run manually (not in CI)
#=
using Plots
let
    t = 1:length(days2do)*24
    depth_labels = ["$(round(ustrip(u"cm", depths[i]); digits=1)) cm" for i in 1:length(depths)]

    # Soil temperature (°C)
    p_st = plot(layout=(2, 5), size=(1400, 600), title=reshape(depth_labels, 1, :))
    for col in 1:length(depths)
        plot!(p_st, t, ustrip.(u"°C", micro_out.soil_temperature[t, col]); sp=col, label="Julia",     color=:red,   ylabel="°C")
        plot!(p_st, t, collect(soiltemps_nmr[t, col]);                     sp=col, label="NicheMapR", color=:black)
    end
    display(p_st)

    # Atmospheric profiles
    p_atm = plot(layout=(4, 2), size=(900, 700))
    plot!(p_atm, t, humidity_matrix[t, 2];                               sp=1, label="Julia",     color=:red,   title="RH 2m")
    plot!(p_atm, t, rh2m_nmr[t];                                         sp=1, label="NicheMapR", color=:black)
    plot!(p_atm, t, wind_matrix[t, 2];                                   sp=2, label="Julia",     color=:red,   title="Wind 2m")
    plot!(p_atm, t, vel2m_nmr[t];                                        sp=2, label="NicheMapR", color=:black)
    plot!(p_atm, t, snow_depth_matrix[t];                                sp=3, label="Julia",     color=:red,   title="Snow depth")
    plot!(p_atm, t, snowdepth_nmr[t];                                    sp=3, label="NicheMapR", color=:black)
    plot!(p_atm, t, snow_density_matrix[t];                              sp=4, label="Julia",     color=:red,   title="Snow density")
    plot!(p_atm, t, snowdensity_nmr[t];                                  sp=4, label="NicheMapR", color=:black)
    display(p_atm)
end
=#
