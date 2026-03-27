using Microclimate
using FluidProperties
using SolarRadiation
using Unitful
using CSV, DataFrames
using Test

testdir = realpath(joinpath(dirname(pathof(Microclimate)), "../test"))

# read in output from NicheMapR
soil_temperature_nmr = DataFrame(CSV.File("$testdir/data/soil_FordDryLake.csv"))[:, 5:14]
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

days = collect(1:Int(length(soil_temperature_nmr[:, 1]) / 24)) # days of year to run (for solar_radiation)
coarse_depths = ((DataFrame(CSV.File("$testdir/data/init_daily/DEP.csv"))[:, 2]) / 100.0)u"m" # coarse soil nodes from R implementation
depths = let n = length(coarse_depths) # expand to fine grid by inserting midpoints between each pair
    result = Vector{eltype(coarse_depths)}(undef, 2n - 1)
    for i in 1:n; result[2i-1] = coarse_depths[i]; end
    for i in 1:n-1; result[2i] = (coarse_depths[i] + coarse_depths[i+1]) / 2; end
    result
end
heights = [microinput[:Usrhyt], microinput[:Refhyt]]u"m" # air nodes for temperature, wind speed and humidity profile
days2do = 30
hours2do = days2do * 24

#TODO make one terrain object via BiophysicalEcologyBase or BiophysicalGrids
micro_terrain = MicroTerrain(;
    elevation = microinput[:ALTT] * 1.0u"m", # elevation (m)
    roughness_height = microinput[:RUF] * 1.0u"m", # roughness height for standard mode TODO dispatch based on roughness pars
    karman_constant = 0.4, # Kármán constant
    dyer_constant = 16.0, # coefficient from Dyer and Hicks for Φ_m (momentum), γ
    viewfactor = 1.0, # view factor to sky
)

solar_terrain = SolarTerrain(;
    slope = (microinput[:slope])*1.0u"°",
    aspect = (microinput[:azmuth])*1.0u"°",
    elevation = (microinput[:ALTT])*1.0u"m",
    horizon_angles = (DataFrame(CSV.File("$testdir/data/init_daily/hori.csv"))[:, 2])*1.0u"°",
    albedo = (DataFrame(CSV.File("$testdir/data/init_daily/REFLS.csv"))[1, 2] * 1.0),
    atmospheric_pressure = atmospheric_pressure((microinput[:ALTT])*1.0u"m"),
    latitude = (microinput[:ALAT] + microinput[:AMINUT] / 60) * 1.0u"°",
    longitude = (microinput[:ALONG] + microinput[:ALMINT] / 60) * 1.0u"°",
)

soil_thermal_model = CampbelldeVriesSoilThermal(;
    de_vries_shape_factor = 0.1, # de Vries shape factor, 0.33 for organic soils, 0.1 for mineral
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
    campbell_b_parameter = (DataFrame(CSV.File("$testdir/data/init_daily/BB.csv"))[:, 2] * 1.0), # set up vector of ground emissivities for each day
    soil_bulk_density2 = (DataFrame(CSV.File("$testdir/data/init_daily/BD.csv"))[:, 2] * 1.0u"Mg/m^3"), # set up vector of ground emissivities for each day
    soil_mineral_density2 = (DataFrame(CSV.File("$testdir/data/init_daily/DD.csv"))[:, 2] * 1.0u"Mg/m^3"), # set up vector of ground emissivities for each day
    # soil moisture plant parameters
    root_density = DataFrame(CSV.File("$testdir/data/init_daily/L.csv"))[:, 2] * u"m/m^3", # root density at each node, mm/m3 (from Campell 1985 Soil Physics with Basic, p. 131) # max depth for water pooling on the surface, mm (to account for runoff)
    root_resistance = microinput[:RW] * u"m^3/kg/s", # resistance per unit length of root, m3 kg-1 s-1
    stomatal_closure_potential = -microinput[:PC] * u"J/kg", # critical leaf water potential for stomatal closure, J kg-1
    leaf_resistance = microinput[:RL] * u"m^4/kg/s", # resistance per unit length of leaf, m4 kg-1 s-1
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
    shade = (DataFrame(CSV.File("$testdir/data/init_daily/Minshades.csv"))[1:days2do, 2] * 1.0) ./ 100.0, # daily shade (fractional)
    soil_wetness = (DataFrame(CSV.File("$testdir/data/init_daily/PCTWET.csv"))[1:days2do, 2] * 1.0) ./ 100.0, # daily soil wetness (fractional)
    surface_emissivity = (DataFrame(CSV.File("$testdir/data/init_daily/SLES.csv"))[1:days2do, 2] * 1.0), # - surface emissivity
    cloud_emissivity = (DataFrame(CSV.File("$testdir/data/init_daily/SLES.csv"))[1:days2do, 2] * 1.0), # - cloud emissivity
    rainfall = ((DataFrame(CSV.File("$testdir/data/init_daily/rain.csv"))[1:days2do, 2] * 1.0))u"kg/m^2",
    deep_soil_temperature = (DataFrame(CSV.File("$testdir/data/init_daily/tannulrun.csv"))[1:days2do, 2] * 1.0)u"°C", # daily deep soil temperatures
    leaf_area_index = (DataFrame(CSV.File("$testdir/data/init_daily/LAI.csv"))[:, 2] * 1.0u"Mg/m^3"), # leaf area indices per day
)

environment_hourly = HourlyTimeseries(;
    pressure = fill(atmospheric_pressure((microinput[:ALTT])*1.0u"m"), length(hours2do)),
    reference_temperature = Float64.(CSV.File("$testdir/data/init_daily/TAIRhr.csv").x[1:hours2do])u"°C",
    reference_humidity = clamp.(Float64.(CSV.File("$testdir/data/init_daily/RHhr.csv").x[1:hours2do]), 0, 100) ./ 100.0,
    reference_wind_speed = clamp.(Float64.(CSV.File("$testdir/data/init_daily/WNhr.csv").x[1:hours2do])u"m/s", 0.1u"m/s", (Inf)u"m/s"),
    global_radiation = Float64.(CSV.File("$testdir/data/init_daily/SOLRhr.csv").x[1:hours2do])u"W/m^2",
    cloud_cover = clamp.(Float64.(CSV.File("$testdir/data/init_daily/CLDhr.csv").x[1:hours2do]), 0, 100) ./ 100.0,
    rainfall = clamp.(Float64.(CSV.File("$testdir/data/init_daily/RAINhr.csv").x[1:hours2do])u"kg/m^2", 0u"kg/m^2", Inf * u"kg/m^2"),
    zenith_angle=nothing,
    longwave_radiation=nothing,
)

solar_model = SolarProblem(; scattered_uv = Bool(Int(microinput[:IUV])))

# now try the simulation function
problem = MicroProblem(;
    # locations, times, depths and heights
    latitude = (microinput[:ALAT] + microinput[:AMINUT] / 60) * 1.0u"°", # latitude
    days = days[1:days2do], # days of year to simulate - TODO leap years
    hours = 0:1:23, # hour of day for solar_radiation # TODO how and in what context would users change this
    depths, # soil nodes - keep spacing close near the surface
    heights, # air nodes for temperature, wind speed and humidity profile
    # Objects defined above
    solar_model,
    solar_terrain,
    micro_terrain, #TODO combine terrains via a generic terrain in BiophysicalEcologyBase
    soil_moisture_model,
    soil_thermal_model,
    environment_minmax,
    environment_daily,
    environment_hourly,
    iterate_day = Int(microinput[:ndmax]), # number of iterations per day
    daily = Bool(Int(microinput[:microdaily])), # doing consecutive days?
    runmoist = Bool(Int(microinput[:runmoist])), # run soil moisture algorithm?
    hourly_rainfall = Bool(Int(microinput[:rainhourly])), # use hourly rainfall?
    spinup = Bool(Int(microinput[:spinup])), # spin-up the first day by iterate_day iterations?
    # intial conditions
    initial_soil_temperature = let coarse = u"K".((DataFrame(CSV.File("$testdir/data/init_daily/soilinit.csv"))[:, 2] * 1.0)u"°C"), n = length(coarse)
        result = Vector{eltype(coarse)}(undef, 2n - 1)
        for i in 1:n; result[2i-1] = coarse[i]; end
        for i in 1:n-1; result[2i] = (coarse[i] + coarse[i+1]) / 2; end
        result
    end,
    initial_soil_moisture = let coarse = Array(DataFrame(CSV.File("$testdir/data/init_daily/moists.csv"))[:, 2]) .* 1.0, n = length(coarse)
        result = Vector{Float64}(undef, 2n - 1)
        for i in 1:n; result[2i-1] = coarse[i]; end
        for i in 1:n-1; result[2i] = (coarse[i] + coarse[i+1]) / 2; end
        result
    end,
    #maximum_surface_temperature = u"K"(microinput[:maxsurf]u"°C")
)

# now try the simulation function
@time micro_out = Microclimate.solve(problem);
# using Profile, ProfileView
# @profile Microclimate.solve(problem)
# ProfileView.view()

# TODO test plotting again at some stage, but it slows down CI a lot
# plot(micro_out)

# subset NicheMapR predictions
vel1cm_nmr = collect(metout_nmr[:, 8]) .* 1u"m/s"
vel2m_nmr = collect(metout_nmr[:, 9]) .* 1u"m/s"
ta1cm_nmr = collect(metout_nmr[:, 4] .+ 273.15) .* 1u"K"
ta2m_nmr = collect(metout_nmr[:, 5] .+ 273.15) .* 1u"K"
rh1cm_nmr = collect(metout_nmr[:, 6]) ./ 100.0
rh2m_nmr = collect(metout_nmr[:, 7]) ./ 100.0
tskyC_nmr = collect(metout_nmr[:, 15]) .* u"°C"

air_temperature_matrix = micro_out.profile.air_temperature
humidity_matrix = micro_out.profile.relative_humidity
wind_matrix = micro_out.profile.wind_speed

coarse_indices = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19] # indices of original 10 coarse nodes within 19-node fine grid
@testset "runmicro comparisons" begin
    @test all(isapprox.(micro_out.soil_temperature[:, coarse_indices[1:3]], u"K".(Matrix(soil_temperature_nmr[1:hours2do, 1:3]) .* u"°C"); rtol=1e-1))
    @test all(isapprox.(micro_out.soil_moisture[:, coarse_indices[1:8]], Matrix(soil_moisture_nmr[1:hours2do, 1:8]); rtol=1e1)) # TODO make better!
    @test all(isapprox.(micro_out.soil_thermal_conductivity[:, coarse_indices[1:3]], Matrix(soil_conductivity_nmr[1:hours2do, 1:3])u"W * m^-1 * K^-1"; rtol=1e1)) # TODO make better!
    @test humidity_matrix[:, 1] ≈ rh1cm_nmr[1:hours2do] rtol=1e-1
    @test humidity_matrix[:, 2] ≈ rh2m_nmr[1:hours2do] rtol=1e-5
    @test wind_matrix[:, 1] ≈ vel1cm_nmr[1:hours2do] rtol=1e-2
    @test wind_matrix[:, 2] ≈ vel2m_nmr[1:hours2do] rtol=1e-5 
    @test u"K".(air_temperature_matrix[:, 1]) ≈ ta1cm_nmr[1:hours2do] rtol=1e-3
    @test u"K".(air_temperature_matrix[:, 2]) ≈ ta2m_nmr[1:hours2do] rtol=1e-5
end

# Visual comparisons — run manually (not in CI)
# using Plots
# let
#     t = 1:hours2do
#     depth_labels = ["$(round(ustrip(u"cm", coarse_depths[i]); digits=1)) cm" for i in 1:length(coarse_depths)]

#     # Soil temperature (°C)
#     p_st = plot(layout=(2, 5), size=(1400, 600), title=reshape(depth_labels, 1, :))
#     for (col, i) in enumerate(coarse_indices)
#         plot!(p_st, t, ustrip.(u"°C", micro_out.soil_temperature[t, i]); sp=col, label="Julia", color=:red, ylabel="°C")
#         plot!(p_st, t, soil_temperature_nmr[t, col];                      sp=col, label="NicheMapR", color=:black)
#     end
#     display(p_st)

#     # Soil moisture (m³/m³)
#     p_sm = plot(layout=(2, 5), size=(1400, 600), title=reshape(depth_labels, 1, :))
#     for (col, i) in enumerate(coarse_indices)
#         plot!(p_sm, t, micro_out.soil_moisture[t, i];    sp=col, label="Julia", color=:red, ylabel="m³/m³")
#         plot!(p_sm, t, soil_moisture_nmr[t, col];        sp=col, label="NicheMapR", color=:black)
#     end
#     display(p_sm)

#     # Soil thermal conductivity (W/m/K)
#     p_tc = plot(layout=(2, 5), size=(1400, 600), title=reshape(depth_labels, 1, :))
#     for (col, i) in enumerate(coarse_indices)
#         plot!(p_tc, t, ustrip.(u"W/m/K", micro_out.soil_thermal_conductivity[t, i]); sp=col, label="Julia", color=:red, ylabel="W/m/K")
#         plot!(p_tc, t, soil_conductivity_nmr[t, col];                                sp=col, label="NicheMapR", color=:black)
#     end
#     display(p_tc)

#     # Atmospheric profiles
#     p_atm = plot(layout=(3, 2), size=(900, 700))
#     plot!(p_atm, t, humidity_matrix[t, 1];                          sp=1, label="Julia",     color=:red,   title="RH 1cm",        ylabel="–")
#     plot!(p_atm, t, rh1cm_nmr[t];                                   sp=1, label="NicheMapR", color=:black)
#     plot!(p_atm, t, humidity_matrix[t, 2];                          sp=2, label="Julia",     color=:red,   title="RH 2m")
#     plot!(p_atm, t, rh2m_nmr[t];                                    sp=2, label="NicheMapR", color=:black)
#     plot!(p_atm, t, ustrip.(u"m/s", wind_matrix[t, 1]);             sp=3, label="Julia",     color=:red,   title="Wind 1cm",      ylabel="m/s")
#     plot!(p_atm, t, ustrip.(u"m/s", vel1cm_nmr[t]);                 sp=3, label="NicheMapR", color=:black)
#     plot!(p_atm, t, ustrip.(u"m/s", wind_matrix[t, 2]);             sp=4, label="Julia",     color=:red,   title="Wind 2m")
#     plot!(p_atm, t, ustrip.(u"m/s", vel2m_nmr[t]);                  sp=4, label="NicheMapR", color=:black)
#     plot!(p_atm, t, ustrip.(u"°C", u"K".(air_temperature_matrix[t, 1])); sp=5, label="Julia", color=:red, title="Air temp 1cm",  ylabel="°C")
#     plot!(p_atm, t, ustrip.(u"°C", ta1cm_nmr[t]);                   sp=5, label="NicheMapR", color=:black)
#     plot!(p_atm, t, ustrip.(u"°C", u"K".(air_temperature_matrix[t, 2])); sp=6, label="Julia", color=:red, title="Air temp 2m")
#     plot!(p_atm, t, ustrip.(u"°C", ta2m_nmr[t]);                    sp=6, label="NicheMapR", color=:black)
#     display(p_atm)
# end
