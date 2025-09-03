using Microclimate
using Unitful
using Plots
using Statistics
using Interpolations
using OrdinaryDiffEq
using CSV, DataFrames, Dates
using Test

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

reference_height = microinput[:Refhyt] * 1.0u"m"
depths = ((DataFrame(CSV.File("test/data/init_monthly/DEP.csv"))[:, 2]) / 100.0)u"m"#[0.0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0, 2.0]u"m" # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
hours = collect(0.:1:24.) # hour of day for solrad
longlat = (DataFrame(CSV.File("test/data/init_monthly/longlat.csv"))[:, 2] * 1.0)
latitude = longlat[2]*1.0u"°" # latitude
longitude = longlat[1]*1.0u"°" # longitude
iuv = Bool(Int(microinput[:IUV])) # this makes it take ages if true!
ndmax = (microinput[:ndmax])
elevation = microinput[:ALTT] * 1.0u"m" # elevation (m)
horizon_angles = (DataFrame(CSV.File("test/data/init_monthly/hori.csv"))[:, 2]) * 1.0u"°"#fill(0.0u"°", 24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
slope = microinput[:slope] * 1.0u"°" # slope (degrees, range 0-90)
aspect = microinput[:azmuth] * 1.0u"°" # aspect (degrees, 0 = North, range 0-360)
ruf = microinput[:RUF] * 1.0u"m" # m roughness height
zh = microinput[:ZH] * 1.0u"m" # m heat transfer roughness height
d0 = microinput[:D0] * 1.0u"m" # zero plane displacement correction factor
# soil properties# soil thermal parameters 
ρ_b_dry = (CSV.File("test/data/init_monthly/soilprop.csv")[1, 1][2]) * 1.0u"Mg/m^3" # dry soil bulk density (Mg/m3)
θ_sat = (CSV.File("test/data/init_monthly/soilprop.csv")[1, 1][3]) * 1.0u"m^3/m^3" # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
λ_m = (CSV.File("test/data/init_monthly/soilprop.csv")[1, 1][4]) * 1.0u"W/m/K" # soil minerals thermal conductivity (W/mC)
c_p_m = (CSV.File("test/data/init_monthly/soilprop.csv")[1, 1][5]) * 1.0u"J/kg/K" # soil minerals specific heat (J/kg-K)
ρ_m = (CSV.File("test/data/init_monthly/soilprop.csv")[1, 1][6]) * 1.0u"Mg/m^3" # soil minerals density (Mg/m3)
# Time varying environmental data
TIMINS = [microinput[:TIMINS1], microinput[:TIMINS2], microinput[:TIMINS3], microinput[:TIMINS4]] # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
TIMAXS = [microinput[:TIMAXS1], microinput[:TIMAXS2], microinput[:TIMAXS3], microinput[:TIMAXS4]] # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
TMINN = (DataFrame(CSV.File("test/data/init_monthly/TMINN.csv"))[:, 2] * 1.0)u"°C" # minimum air temperatures (°C)
TMAXX = (DataFrame(CSV.File("test/data/init_monthly/TMAXX.csv"))[:, 2] * 1.0)u"°C" # maximum air temperatures (°C)
RHMINN = (DataFrame(CSV.File("test/data/init_monthly/RHMINN.csv"))[:, 2] * 1.0) # min relative humidity (%)
RHMAXX = (DataFrame(CSV.File("test/data/init_monthly/RHMAXX.csv"))[:, 2] * 1.0) # max relative humidity (%)
WNMINN = (DataFrame(CSV.File("test/data/init_monthly/WNMINN.csv"))[:, 2] * 1.0)u"m/s" # min wind speed (m/s)
WNMAXX = (DataFrame(CSV.File("test/data/init_monthly/WNMAXX.csv"))[:, 2] * 1.0)u"m/s" # max wind speed (m/s)
CCMINN = (DataFrame(CSV.File("test/data/init_monthly/CCMINN.csv"))[:, 2] * 1.0) # min cloud cover (%)
CCMAXX = (DataFrame(CSV.File("test/data/init_monthly/CCMAXX.csv"))[:, 2] * 1.0) # max cloud cover (%)
RAINFALL = ((DataFrame(CSV.File("test/data/init_monthly/rain.csv"))[:, 2] * 1.0) / 1000)u"m" # monthly mean rainfall (mm)
SoilMoist = (DataFrame(CSV.File("test/data/init_monthly/moists.csv"))[:, 2:13] .* 1.0)#fill(0.0, ndays)
daily = Bool(Int(microinput[:microdaily]))
runmoist = Bool(Int(microinput[:runmoist]))

# creating the arrays of environmental variables that are assumed not to change with month for this simulation
ndays = length(days)
shades = (DataFrame(CSV.File("test/data/init_monthly/MINSHADES.csv"))[:, 2] * 1.0) # daily shade (%)
sles = (DataFrame(CSV.File("test/data/init_monthly/SLES.csv"))[:, 2] * 1.0) # set up vector of ground emissivities for each day
albedos = (DataFrame(CSV.File("test/data/init_monthly/REFLS.csv"))[:, 2] * 1.0) # set up vector of soil albedos for each day (decimal %)
pctwets = (DataFrame(CSV.File("test/data/init_monthly/PCTWET.csv"))[:, 2] * 1.0) # set up vector of soil wetness for each day (%)
tannul = mean(Unitful.ustrip.(vcat(TMAXX, TMINN)))u"°C" # annual mean temperature for getting monthly deep soil temperature (°C)
tannulrun = fill(tannul, ndays) # monthly deep soil temperature (2m) (°C)

# defining view factor based on horizon angles
viewfactor = 1 - sum(sin.(horizon_angles)) / length(horizon_angles) # convert horizon angles to radians and calc view factor(s)

# Soil properties
# set up a profile of soil properites with depth for each day to be run
numtyps = 1 # number of soil types
numnodes = length(depths) # number of soil nodes
nodes_day = zeros(numnodes, ndays) # array of all possible soil nodes
nodes_day[1, 1:ndays] .= 10 # deepest node for first substrate type
# Create an empty 10×5 matrix that can store any type (including different units)
soilprops = Matrix{Any}(undef, numnodes, 5)
# Fill row 1 (top layer) with the defined values
soilprops[1, 1] = ρ_b_dry
soilprops[1, 2] = θ_sat
soilprops[1, 3] = λ_m
soilprops[1, 4] = c_p_m
soilprops[1, 5] = ρ_m
# Copy the same properties to all other layers (if desired)
for i in 2:numnodes
    soilprops[i, :] .= soilprops[1, :]
end
soillayers = init_soillayers(numnodes)  # only once
soilinit = (DataFrame(CSV.File("test/data/init_monthly/soilinit.csv"))[1:numnodes, 2] * 1.0)u"°C" # set up vector of soil wetness for each day (%)
∑phase = zeros(Float64, numnodes)u"J"

# compute solar radiation
solrad_out = solrad(;
    days,
    hours,
    latitude,
    elevation,
    horizon_angles,
    slope,
    aspect,
    albedos,
    iuv,
)
solrad_out.Zenith[solrad_out.Zenith.>90u"°"] .= 90u"°"
solrad_out.ZenithSlope[solrad_out.ZenithSlope.>90u"°"] .= 90u"°"

# interpolate air temperature to hourly
TAIRs, VELs, RHs, CLDs = hourly_vars(
    TMINN,
    TMAXX,
    WNMINN,
    WNMAXX,
    RHMINN,
    RHMAXX,
    CCMINN,
    CCMAXX,
    solrad_out,
    TIMINS,
    TIMAXS,
    daily
)
RHs[RHs.>100] .= 100
CLDs[CLDs.>100] .= 100

TAIRs25 = TAIRs # keeping a copy to get mean monthly over 25 hrs as in fortran
skip25 = setdiff(1:length(solrad_out.Zenith), 25:25:length(solrad_out.Zenith))
ZENRs = solrad_out.Zenith[skip25] # remove every 25th output
ZSLs = solrad_out.ZenithSlope[skip25] # remove every 25th output
SOLRs = solrad_out.Global[skip25] # remove every 25th output
TAIRs = TAIRs[skip25]
VELs = VELs[skip25]
RHs = RHs[skip25]
CLDs = CLDs[skip25]
# Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
SOLRs = SOLRs .* (0.36 .+ 0.64 * (1.0 .- (CLDs / 100.0))) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992

# simulate a day
iday = 1
nsteps = length(hours)
T_soils = Array{Float64}(undef, nsteps, numnodes)u"K"
soiltemps = nothing
T0 = u"K".(soilinit)
∑phase = zeros(Float64, numnodes)u"J"

for iday in 1:ndays
    #iday=2
    sub = (iday*24-24+1):(iday*24)#(iday*24-23):(iday*24)
    sub2 = (iday*25-25+1):(iday*25) # for getting mean monthly over the 25 hrs as in fortran version
    shade = shades[iday] # daily shade (%)
    sle = sles[iday] # set up vector of ground emissivities for each day
    albedo = albedos[iday]
    slep = sle # - cloud emissivity
    pctwet = pctwets[iday] # set up vector of soil wetness for each day
    tdeep = u"K"(tannulrun[iday]) # annual mean temperature for getting daily deep soil temperature (°C)
    nodes = nodes_day[:, iday]

    # get today's weather
    SOLR = SOLRs[sub]
    ZENR = ZENRs[sub]
    ZSL = ZSLs[sub]
    TAIR = TAIRs[sub]
    VEL = VELs[sub]
    RH = RHs[sub]
    CLD = CLDs[sub]

    # Initial conditions
    if !daily
        soilinit = u"K"(mean(ustrip(TAIRs25[sub2]))u"°C") # make initial soil temps equal to mean annual temperature
        T0 = fill(soilinit, numnodes)
        ∑phase = zeros(Float64, numnodes)u"J"
    end
    θ_soil = SoilMoist[:, iday]#collect(fill(SoilMoist, numnodes)) # initial soil moisture, constant for now
    #T0[numnodes] = tdeep
    # create forcing weather variable splines
    #tspan = 0.:60:1440
    tspan = 0.0:60:(60*24)
    tmin = tspan .* u"minute"
    interpSOLR = interpolate([SOLR; SOLR[end]], BSpline(Linear()))
    interpZENR = interpolate([ZENR; ZENR[end]], BSpline(Linear()))
    interpZSL = interpolate([ZSL; ZSL[end]], BSpline(Linear()))
    interpTAIR = interpolate(u"K".([TAIR; TAIR[end]]), BSpline(Linear()))
    interpVEL = interpolate([VEL; VEL[end]], BSpline(Linear()))
    interpRH = interpolate([RH; RH[end]], BSpline(Linear()))
    interpCLD = interpolate([CLD; CLD[end]], BSpline(Linear()))
    SOLRt = scale(interpSOLR, tspan)
    ZENRt = scale(interpZENR, tspan)
    ZSLt = scale(interpZSL, tspan)
    TAIRt = scale(interpTAIR, tspan)
    VELt = scale(interpVEL, tspan)
    RHt = scale(interpRH, tspan)
    CLDt = scale(interpCLD, tspan)

    # Parameters
    params = MicroParams(
        soilprops,
        depths,
        reference_height,
        ruf,
        d0,
        zh,
        slope,
        shade,
        viewfactor,
        elevation,
        albedo,
        sle,
        slep, # TODO check if this is what it should be - sle vs. slep (set as 1 in PAR in Fortran but then changed to user SLE later)
        pctwet,
        nodes,
        tdeep,
        θ_soil,
        runmoist,
    )

    forcing = MicroForcing(;
        SOLRt,
        ZENRt,
        ZSLt,
        TAIRt,
        VELt,
        RHt,
        CLDt,
    )

    input = MicroInputs(
        params,
        forcing,
        soillayers
    )

    # loop through hours of day
    niter = ndmax # number of interations for steady periodic
    iter = 0
    for iter in 1:niter
        #iter = 1
        step = 2
        T_soils[1, :] = T0
        for i in 1:length(hours)-1
            #i=1
            tspan = ((0.0 + (step - 2) * 60)u"minute", (60.0 + (step - 2) * 60)u"minute")  # 1 hour
            #tspan = 0.0:60:(60*24)
            prob = ODEProblem(soil_energy_balance!, T0, tspan, input)
            sol = solve(prob, Tsit5(); saveat=60.0u"minute", reltol=1e-6u"K", abstol=1e-8u"K")
            soiltemps = hcat(sol.u...)
            if iter == niter #|| iday == ndays
                ∑phase, qphase, T0 = phase_transition(soiltemps[:, 2], soiltemps[:, 1], ∑phase, θ_soil, depths)
            else
               T0 = soiltemps[:, 2]
            end
            T_soils[step, :] = T0
            step += 1
        end
    end

    t = hours[1:24]u"hr"
    plt = plot(t, u"°C".(T_soils)[1:24, :], xlabel="Time", ylabel="Soil Temperature", lw=2, legend=false, ylims=[-20, 50])
    plot!(plt, t, Matrix(soiltemps_NMR[sub, :]);
        xlabel="time", ylabel="soil temperature", lw=2,
        label=string.(depths'), linestyle=:dash, linecolor="grey"
    )
    #plot!(plt, t, u"°C".(micro_out.T_soils[sub,:]), xlabel="Time", ylabel="Soil Temperature", lw=2, legend=false, ylims=[-20, 50])


#     # Display the plot
    display(plt)


    #     # now get wind air temperature and humidity profiles
    # heights = [0.01] .* u"m"
    # profiles = Vector{NamedTuple}(undef, nsteps)  # or whatever you have from your loop
    # for i in 1:nsteps
    #     profiles[i] = get_profile(
    #         reference_height = reference_height,
    #         ruf = ruf,
    #         zh = zh,
    #         d0 = d0,
    #         TAREF = TAIR[i],
    #         VREF = VEL[i],
    #         rh = RH[i],
    #         D0cm = u"°C"(T_soils'[1, i]),  # top layer temp at time i
    #         ZEN = ZENR[i],
    #         heights = heights,
    #         elevation = elevation
    #     )
    # end
    # Tskys = zeros(24)u"K"  # or whatever you have from your loop
    # for i in 1:nsteps
    #     Tskys[i] = Microclimate.get_longwave(
    #         elevation = elevation, 
    #         rh = RH[i], 
    #         tair = TAIR[i], 
    #         tsurf = u"°C"(T_soils'[1, i]), 
    #         slep = slep, 
    #         sle = sle, 
    #         cloud = CLD[i], 
    #         viewfactor = viewfactor,
    #         shade = shade
    #         ).Tsky
    # end
    # # subset NicheMapR predictions
    # vel1cm_NMR = collect(metout_NMR[sub, 8]).*1u"m/s"
    # vel2m_NMR = collect(metout_NMR[sub, 9]).*1u"m/s"
    # ta1cm_NMR = collect(metout_NMR[sub, 4] .+ 273.15).*1u"K"
    # ta2m_NMR = collect(metout_NMR[sub, 5] .+ 273.15).*1u"K"
    # rh1cm_NMR = collect(metout_NMR[sub, 6])
    # rh2m_NMR = collect(metout_NMR[sub, 7])
    # tskyC_NMR = collect(metout_NMR[sub, 15]).*u"°C"

    # plot(t, u"°C".(Tskys), xlabel="time", ylabel="sky temperature", lw=2)
    # plot!(t, tskyC_NMR, xlabel="time", ylabel="sky temperature", lw=2, linestyle = :dash, linecolor="grey")

    # VEL_matrix = hcat([getfield.(profiles, :VELs)[i] for i in 1:length(profiles)]...)    # velocities at all time steps and heights
    # TA_matrix   = hcat([getfield.(profiles, :TAs)[i]   for i in 1:length(profiles)]...)  # air temperatures at all time steps and heights
    # RH_matrix   = hcat([getfield.(profiles, :RHs)[i]   for i in 1:length(profiles)]...)       # relative humidity at all time steps and heights
    # Qconv_vec   = [getfield(p, :QCONV) for p in profiles]    # convective fluxes at all time steps
    # ustar_vec   = [getfield(p, :USTAR) for p in profiles]    # friction velocities
    # heights_vec = getfield(profiles[1], :heights)                     # constant across time
    # plot(t, VEL_matrix', xlabel="time", ylabel="wind speed", lw=2, label = string.(first(heights_vec, 11)'))
    # plot!(t, vel1cm_NMR, xlabel="time", ylabel="wind speed", lw=2, label = "1cm NMR", linestyle = :dash, linecolor="grey")
    # plot!(t, vel2m_NMR, xlabel="time", ylabel="wind speed", lw=2, label = "200cm NMR", linestyle = :dash, linecolor="grey")

    # plot(t, TA_matrix', xlabel="time", ylabel="air temperature", lw=2, label = string.(first(heights_vec, 11)'))
    # plot!(t, ta1cm_NMR, xlabel="time", ylabel="air temperature", lw=2, label = "1cm NMR", linestyle = :dash, linecolor="grey")
    # plot!(t, ta2m_NMR, xlabel="time", ylabel="air temperature", lw=2, label = "200cm NMR", linestyle = :dash, linecolor="grey")

    # plot(t, RH_matrix', xlabel="time", ylabel="humidity (%)", lw=2, label = string.(first(heights_vec, 11)'))
    # plot!(t, rh1cm_NMR, xlabel="time", ylabel="humidity (%)", lw=2, label = "1cm NMR", linestyle = :dash, linecolor="grey")
    # plot!(t, rh2m_NMR, xlabel="time", ylabel="humidity (%)", lw=2, label = "200cm NMR", linestyle = :dash, linecolor="grey")

end

# now try the simulation function 'runmicro' which has the above inputs as default
micro_out = runmicro(;)
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