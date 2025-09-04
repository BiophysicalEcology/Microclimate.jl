using Microclimate
using Plots
using Unitful
using CSV, DataFrames

# read in output from NicheMapR and input variables
soiltemps_NMR = (DataFrame(CSV.File("test/data/soil_monthly.csv"))[:, 4:13]) .* u"°C"
metout_NMR = DataFrame(CSV.File("test/data/metout_monthly.csv"))
microinput_vec = DataFrame(CSV.File("test/data/init_monthly/microinput.csv"))[:, 2]

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

longlat = (DataFrame(CSV.File("test/data/init_monthly/longlat.csv"))[:, 2] * 1.0)
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
LAIs = fill(0.1, length(days))
depths = ((DataFrame(CSV.File("test/data/init_monthly/DEP.csv"))[:, 2]) / 100.0)u"m"
heights = [1.0,]u"cm" # air nodes for temperature, wind speed and humidity profile

micro_out = runmicro(;
    # locations, times, depths and heights
    latitude = longlat[2]*1.0u"°",
    days, # days of year for solrad
    hours = hours = collect(0.:1:24.), # hour of day for solrad
    reference_height = microinput[:Refhyt] * 1.0u"m",
    depths,
    heights, # air nodes for temperature, wind speed and humidity profile
    # terrain
    elevation = microinput[:ALTT] * 1.0u"m",
    horizon_angles = horizon_angles = (DataFrame(CSV.File("test/data/init_monthly/hori.csv"))[:, 2]) * 1.0u"°",
    slope = microinput[:slope] * 1.0u"°",
    aspect = microinput[:azmuth] * 1.0u"°",
    roughness_height = microinput[:RUF] * 1.0u"m", # roughness height for standard mode TODO dispatch based on roughness pars
    zh = microinput[:ZH] * 1.0u"m", # heat transfer roughness height for Campbell and Norman mode
    d0 = microinput[:D0] * 1.0u"m", # zero plane displacement correction factor
    # soil thermal parameters 
    soil_mineral_conductivity = (CSV.File("test/data/init_monthly/soilprop.csv")[1, 1][4]) * 1.0u"W/m/K", # soil minerals thermal conductivity (W/mC)
    soil_mineral_density = (CSV.File("test/data/init_monthly/soilprop.csv")[1, 1][6]) * 1.0u"Mg/m^3", # soil minerals density (Mg/m3)
    soil_mineral_heat_capacity = c_p_m = (CSV.File("test/data/init_monthly/soilprop.csv")[1, 1][5]) * 1.0u"J/kg/K", # soil minerals specific heat (J/kg-K)
    soil_bulk_density = (CSV.File("test/data/init_monthly/soilprop.csv")[1, 1][2]) * 1.0u"Mg/m^3", # dry soil bulk density (Mg/m3)
    soil_saturation_moisture = (CSV.File("test/data/init_monthly/soilprop.csv")[1, 1][3]) * 1.0u"m^3/m^3", # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    # daily environmental vectors
    albedos = (DataFrame(CSV.File("test/data/init_monthly/REFLS.csv"))[:, 2] * 1.0), # substrate albedo (decimal %)
    shades = (DataFrame(CSV.File("test/data/init_monthly/MINSHADES.csv"))[:, 2] * 1.0), # daily shade (%)
    pctwets = (DataFrame(CSV.File("test/data/init_monthly/PCTWET.csv"))[:, 2] * 1.0),
    sles = (DataFrame(CSV.File("test/data/init_monthly/SLES.csv"))[:, 2] * 1.0), # - surface emissivity
    daily_rainfall = ((DataFrame(CSV.File("test/data/init_monthly/rain.csv"))[:, 2] * 1.0) / 1000)u"kg/m^2", # monthly total rainfall
    air_temperature_min = (DataFrame(CSV.File("test/data/init_monthly/TMINN.csv"))[:, 2] * 1.0)u"°C", # minimum air temperatures
    air_temperature_max = (DataFrame(CSV.File("test/data/init_monthly/TMAXX.csv"))[:, 2] * 1.0)u"°C", # maximum air temperatures
    wind_min = (DataFrame(CSV.File("test/data/init_monthly/WNMINN.csv"))[:, 2] * 1.0)u"m/s", # min wind speed (m/s)
    wind_max = (DataFrame(CSV.File("test/data/init_monthly/WNMAXX.csv"))[:, 2] * 1.0)u"m/s", # max wind speed (m/s)
    humidity_min = (DataFrame(CSV.File("test/data/init_monthly/RHMINN.csv"))[:, 2] * 1.0), # min relative humidity (%)
    humidity_max = (DataFrame(CSV.File("test/data/init_monthly/RHMAXX.csv"))[:, 2] * 1.0), # max relative humidity (%)
    cloud_min = (DataFrame(CSV.File("test/data/init_monthly/CCMINN.csv"))[:, 2] * 1.0), # min cloud cover (%)
    cloud_max = (DataFrame(CSV.File("test/data/init_monthly/CCMAXX.csv"))[:, 2] * 1.0), # max cloud cover (%)
    minima_times = [microinput[:TIMINS1], microinput[:TIMINS2], microinput[:TIMINS3], microinput[:TIMINS4]], # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
    maxima_times = [microinput[:TIMAXS1], microinput[:TIMAXS2], microinput[:TIMAXS3], microinput[:TIMAXS4]], # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
    # intial conditions
    initial_soil_temperature = u"K".((DataFrame(CSV.File("test/data/init_monthly/soilinit.csv"))[1:length(depths), 2] * 1.0)u"°C"), # initial soil temperature
    initial_soil_moisture = (Array(DataFrame(CSV.File("test/data/init_monthly/moists.csv"))[1, 2:13]) .* 1.0), # initial soil moisture
    leaf_area_index = fill(0.1, length(days)),
    iterate_day = microinput[:ndmax], # number of iterations per day
    daily = Bool(Int(microinput[:microdaily])), # doing consecutive days?
    runmoist = Bool(Int(microinput[:runmoist])), # run soil moisture algorithm?
    spinup = Bool(Int(microinput[:spinup])), # spin-up the first day by iterate_day iterations?
    iuv = Bool(Int(microinput[:IUV])), # this makes it take ages if true!
)

using ProfileView
@profview micro_out = runmicro(;)

plot(micro_out.soil_temperature, legend=false)
plot!(Matrix(soiltemps_NMR);
        xlabel="time", ylabel="soil temperature", lw=2,
        linestyle=:dash, linecolor="grey"
    )
dayplot=2
sub=((dayplot-1)*24+1):(dayplot*24)
plt = plot(u"°C".(micro_out.soil_temperature[sub,:]), xlabel="Time", ylabel="Soil Temperature", lw=2, legend=false, ylims=[-20, 50])
plot!(plt, Matrix(soiltemps_NMR[sub,:]);
        xlabel="time", ylabel="soil temperature", lw=2,
        linestyle=:dash, linecolor="grey"
    )

# subset NicheMapR predictions
vel1cm_NMR = collect(metout_NMR[:, 8]) .* 1u"m/s"
vel2m_NMR = collect(metout_NMR[:, 9]) .* 1u"m/s"
ta1cm_NMR = collect(metout_NMR[:, 4] .+ 273.15) .* 1u"K"
ta2m_NMR = collect(metout_NMR[:, 5] .+ 273.15) .* 1u"K"
rh1cm_NMR = collect(metout_NMR[:, 6])
rh2m_NMR = collect(metout_NMR[:, 7])
tskyC_NMR = collect(metout_NMR[:, 15]) .* u"°C"

plot(u"°C".(micro_out.sky_temperature), xlabel="time", ylabel="sky temperature", lw=2)
plot!(tskyC_NMR, xlabel="time", ylabel="sky temperature", lw=2, linestyle=:dash, linecolor="grey")

plot(micro_out.wind_speed[:, 2:3], xlabel="time", ylabel="wind speed", lw=2)
plot!(vel1cm_NMR, xlabel="time", ylabel="wind speed", lw=2, label="1cm NMR", linestyle=:dash, linecolor="grey")
plot!(vel2m_NMR, xlabel="time", ylabel="wind speed", lw=2, label="200cm NMR", linestyle=:dash, linecolor="grey")

plot(micro_out.air_temperature[:, 2:3], xlabel="time", ylabel="air temperature", lw=2)
plot!(ta1cm_NMR, xlabel="time", ylabel="air temperature", lw=2, label="1cm NMR", linestyle=:dash, linecolor="grey")
plot!(ta2m_NMR, xlabel="time", ylabel="air temperature", lw=2, label="200cm NMR", linestyle=:dash, linecolor="grey")

plot(micro_out.relative_humidity[:, 2:3], xlabel="time", ylabel="humidity (%)", lw=2)
plot!(rh1cm_NMR, xlabel="time", ylabel="humidity (%)", lw=2, label="1cm NMR", linestyle=:dash, linecolor="grey")
plot!(rh2m_NMR, xlabel="time", ylabel="humidity (%)", lw=2, label="200cm NMR", linestyle=:dash, linecolor="grey")
