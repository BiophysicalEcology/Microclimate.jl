#using Microclimate
using Unitful
using Unitful: °, rad, R, kg, m
using Plots
using Statistics
using Interpolations
using DifferentialEquations
using CSV, DataFrames, Dates

# read in output from Norman
soiltemps_NMR = (DataFrame(CSV.File("data/soil.csv"))[:, 4:13]).*u"°C"
metout_NMR = DataFrame(CSV.File("data/metout.csv"))

DEP = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm" # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
refhyt = 2u"m"
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
hours = collect(0.:1:24.) # hour of day for solrad
lat = 43.1379° # latitude
iuv = false # this makes it take ages if true!
elev = 226.0u"m" # elevation (m)
hori = fill(0.0°, 24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
slope = 0.0° # slope (degrees, range 0-90)
aspect = 0.0° # aspect (degrees, 0 = North, range 0-360)
refl = 0.10 # substrate solar reflectivity (decimal %)
shade = 0.0 # % shade cast by vegetation
pctwet = 0.0 # % surface wetness
sle = 0.96 # - surface emissivity
ruf = 0.004u"m" # m roughness height
zh = 0u"m" # m heat transfer roughness height
d0 = 0u"m" # zero plane displacement correction factor

# soil properties

# soil thermal parameters 
λ_m = 1.25u"W/m/K" # soil minerals thermal conductivity (W/mC)
ρ_m = 2.560u"Mg/m^3" # soil minerals density (Mg/m3)
cp_m = 870.0u"J/kg/K" # soil minerals specific heat (J/kg-K)
ρ_b_dry = 1.3u"Mg/m^3" # dry soil bulk density (Mg/m3)
θ_sat = 0.26u"m^3/m^3" # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)

# soil moisture model parameters
runmoist = true
PE = fill(0.7, 19)u"J/kg" # air entry potential (J/kg) (19 values descending through soil for specified soil nodes in parameter DEP and points half way between)
KS = fill(0.0058, 19)u"kg*s/m^3" # saturated conductivity, (kg s/m3) (19 values descending through soil for specified soil nodes in parameter DEP and points half way between)
BB = fill(1.7, 19) # Campbell's soil 'b' parameter (-) (19 values descending through soil for specified soil nodes in parameter DEP and points half way between)
BD = fill(ρ_b_dry, 19)u"Mg/m^3" # soil bulk density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter DEP and points half way between)
DD = fill(ρ_m, 19)u"Mg/m^3" # soil density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter DEP and points half way between)
L = [0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0]*10000u"mm/m^3" # root density at each node, mm/m3 (from Campell 1985 Soil Physics with Basic, p. 131)
rw =  2.5E+10u"m^3/kg/s" # resistance per unit length of root, m3 kg-1 s-1
pc = -1500.0u"J/kg" # critical leaf water potential for stomatal closure, J kg-1
rl = 2.0e6u"m^4/kg/s" # resistance per unit length of leaf, m3 kg-1 s-1
sp = 10.0 # stability parameter, -
r1 = 0.001u"m" # root radius, m
im = 1e-6u"kg/m^2/s" # maximum overall mass balance error allowed, kg
maxcount = 500
timestep = 360.0u"s"

# Time varying environmental data
TIMINS = [0, 0, 1, 1] # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
TIMAXS = [1, 1, 0, 0] # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
TMINN = [-14.3, -12.1, -5.1, 1.2, 6.9, 12.3, 15.2, 13.6, 8.9, 3, -3.2, -10.6]u"°C" # minimum air temperatures (°C)
TMAXX = [-3.2, 0.1, 6.8, 14.6, 21.3, 26.4, 29, 27.7, 23.3, 16.6, 7.8, -0.4]u"°C" # maximum air temperatures (°C)
RHMINN = [50.2, 48.4, 48.7, 40.8, 40, 42.1, 45.5, 47.3, 47.6, 45, 51.3, 52.8] # min relative humidity (%)
RHMAXX = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0] # max relative humidity (%)
WNMINN = 0.1 .* [4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8]u"m/s" # min wind speed (m/s)
WNMAXX = [4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8]u"m/s" # max wind speed (m/s)
CCMINN = 0.0 .* [50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1] # min cloud cover (%)
CCMAXX = 0.0 .* [50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1] # max cloud cover (%)
RAINFALL = ([28, 28.2, 54.6, 79.7, 81.3, 100.1, 101.3, 102.5, 89.7, 62.4, 54.9, 41.2]) / 1000u"m" # monthly mean rainfall (mm)
SoilMoist = fill(0.2, length(TMAXX))
LAIs = fill(0.1, length(TMAXX))

daily = false

# creating the arrays of environmental variables that are assumed not to change with month for this simulation
ndays = length(days)
SHADES = fill(shade, ndays) # daily shade (%)
SLES = fill(sle, ndays) # set up vector of ground emissivities for each day
REFLS = fill(refl, ndays) # set up vector of soil reflectances for each day
PCTWETS = fill(pctwet, ndays) # set up vector of soil wetness for each day
tannul = mean(Unitful.ustrip.(vcat(TMAXX, TMINN)))u"°C" # annual mean temperature for getting monthly deep soil temperature (°C)
tannulrun = fill(tannul, ndays) # monthly deep soil temperature (2m) (°C)

# defining view factor based on horizon angles
viewf = 1 - sum(sin.(hori)) / length(hori) # convert horizon angles to radians and calc view factor(s)

# Soil properties
# set up a profile of soil properites with depth for each day to be run
numtyps = 1 # number of soil types
numnodes = length(DEP) # number of soil nodes
nodes_day = zeros(numnodes, ndays) # array of all possible soil nodes
nodes_day[1, 1:ndays] .= 10 # deepest node for first substrate type
# Create an empty 10×5 matrix that can store any type (including different units)
soilprops = Matrix{Any}(undef, numnodes, 5)
# Fill row 1 (top layer) with the defined values
soilprops[1, 1] = ρ_b_dry
soilprops[1, 2] = θ_sat
soilprops[1, 3] = λ_m
soilprops[1, 4] = cp_m
soilprops[1, 5] = ρ_m
# Copy the same properties to all other layers (if desired)
for i in 2:numnodes
    soilprops[i, :] .= soilprops[1, :]
end

# compute solar radiation (need to make refl time varying)
solrad_out = solrad(days = days, hours = hours, lat = lat, elev = elev, hori = hori, slope = slope, aspect = aspect, refl = refl, iuv = iuv)

# interpolate air temperature to hourly
TAIRs, WNs, RHs, CLDs = hourly_vars(
    TMINN=TMINN, 
    TMAXX=TMAXX, 
    WNMINN=WNMINN, 
    WNMAXX=WNMAXX, 
    RHMINN=RHMINN, 
    RHMAXX=RHMAXX, 
    CCMINN=CCMINN, 
    CCMAXX=CCMAXX,
    solrad_out=solrad_out, 
    TIMINS=TIMINS, 
    TIMAXS=TIMAXS, 
    daily=daily
    )
    
# Initial conditions
soilinit = u"K"(mean(ustrip(TAIR))u"°C") # make initial soil temps equal to mean daily temperature
T0 = fill(soilinit, numnodes)
θ_soil = collect(fill(SoilMoist[iday], numnodes)) # initial soil moisture, constant for now
# intitial soil moisture
θ_soil18 = similar(θ_soil, 18)  # preallocate vector of length 18
j = 1
for i in 1:18
    if isodd(i)
        θ_soil18[i] = θ_soil[j]
        j += 1
    else
        θ_soil18[i] = θ_soil18[i-1]
    end
end
θ_soil0 = θ_soil18
θ_soil0[1] = 1 - BD[1] / DD[1]

# output arrays
nsteps = length(days) * (length(hours) - 1)
T_soils = Array{Float64}(undef, nsteps, numnodes)u"K"
θ_soils = Array{Float64}(undef, nsteps, numnodes)
ψ_soils = Array{Float64}(undef, nsteps, numnodes)u"J/kg"
rh_soils = Array{Float64}(undef, nsteps, numnodes)

# simulate all days
step = 0
for j in 1:length(days)
    iday = j
    LAI = LAIs[iday]
    # vel1cm_NMR = collect(metout_NMR[(iday*24-23):(iday*24), 8]) .* 1u"m/s"
    # vel2m_NMR = collect(metout_NMR[(iday*24-23):(iday*24), 9]) .* 1u"m/s"
    # ta1cm_NMR = collect(metout_NMR[(iday*24-23):(iday*24), 4] .+ 273.15) .* 1u"K"
    # ta2m_NMR = collect(metout_NMR[(iday*24-23):(iday*24), 5] .+ 273.15) .* 1u"K"
    # rh1cm_NMR = collect(metout_NMR[(iday*24-23):(iday*24), 6])
    # rh2m_NMR = collect(metout_NMR[(iday*24-23):(iday*24), 7])

    sub = (iday*25-24):(iday*25)
    REFL = REFLS[iday]
    SHADE = SHADES[iday] # daily shade (%)
    SLE = SLES[iday] # set up vector of ground emissivities for each day
    PCTWET = PCTWETS[iday] # set up vector of soil wetness for each day
    tdeep = u"K"(tannulrun[iday]) # annual mean temperature for getting daily deep soil temperature (°C)
    nodes = nodes_day[:, iday]

    # get today's weather
    SOLR = solrad_out.Global[sub]
    ZENR = solrad_out.Zenith[sub]
    ZSL = solrad_out.ZenithSlope[sub]
    TAIR = TAIRs[sub]
    VEL = WNs[sub]
    RH = RHs[sub]
    CLD = CLDs[sub]

    # create forcing weather variable splines
    tspan = 0.0:60:1440
    tmin = tspan .* u"minute"
    interpSOLR = interpolate(SOLR, BSpline(Cubic(Line(OnGrid()))))
    interpZENR = interpolate(ZENR, BSpline(Cubic(Line(OnGrid()))))
    interpZSL = interpolate(ZSL, BSpline(Cubic(Line(OnGrid()))))
    interpTAIR = interpolate(u"K".(TAIR), BSpline(Cubic(Line(OnGrid()))))
    interpVEL = interpolate(VEL, BSpline(Cubic(Line(OnGrid()))))
    interpRH = interpolate(RH, BSpline(Cubic(Line(OnGrid()))))
    interpCLD = interpolate(CLD, BSpline(Cubic(Line(OnGrid()))))
    SOLRt = scale(interpSOLR, tspan)
    ZENRt = scale(interpZENR, tspan)
    ZSLt = scale(interpZSL, tspan)
    TAIRt = scale(interpTAIR, tspan)
    VELt = scale(interpVEL, tspan)
    RHt = scale(interpRH, tspan)
    CLDt = scale(interpCLD, tspan)

    forcing = MicroForcing(
        SOLRt=SOLRt,
        ZENRt=ZENRt,
        ZSLt=ZSLt,
        TAIRt=TAIRt,
        VELt=VELt,
        RHt=RHt,
        CLDt=CLDt
    )

    # initial soil temperatures
    #T0 = fill(soilinit, numnodes)
    #T0[numnodes]=tdeep
    #T0[numnodes-2]=(u"K"(TMINN[iday])+u"K"(TMAXX[iday]))/2.0
    #T0[numnodes-1]=(T0[numnodes]+T0[numnodes-2])/2.0

    # loop through hours of day
    for i in 1:length(hours)-1
            # Parameters
    params = MicroParams(
        soilprops=soilprops,
        dep=DEP,
        refhyt=refhyt,
        ruf=ruf,
        d0=d0,
        zh=zh,
        slope=slope,
        shade=shade,
        viewf=viewf,
        elev=elev,
        refl=refl,
        sle=sle,
        slep=sle, # check if this is what it should be - sle vs. slep (set as 1 in PAR in Fortran but then changed to user SLE later)
        pctwet=pctwet,
        nodes=nodes,
        tdeep=tdeep,
        θ_soil=θ_soil,
        runmoist=false,
        runsnow=false
    )
    input = MicroInput(
        params,
        forcing
    )

        tspan = ((0.0 + (i - 1) * 60)u"minute", (60.0 + (i - 1) * 60)u"minute")  # 1 hour
        # if j == 1 && i == 1
        #     niter_soil = 3
        # else
        #     niter_soil = 1
        # end
        prob = ODEProblem(soil_energy_balance!, T0, tspan, input)
        sol = solve(prob, Tsit5(); saveat=60.0u"minute")
        soiltemps = hcat(sol.u...)
        T0 = soiltemps[:, 2] # new initial soil temps

        # compute scalar profiles
        heights = [0.01] .* u"m"
        profile_out = get_profile(
            refhyt=refhyt,
            ruf=ruf,
            zh=zh,
            d0=d0,
            TAREF=TAIR[i],
            VREF=VEL[i],
            rh=RH[i],
            D0cm=u"°C"(T0[1]),  # top layer temp at time i
            ZEN=ZENR[i],
            heights=heights,
            elev=elev,
            warn=true
        )

        # convection
        qconv = profile_out.QCONV

        # evaporation
        P_atmos = get_pressure(elev)
        rh_loc = profile_out.RHs[2]
        hc = max(abs(qconv / (T0[1] - u"K"(TAIR[i]))), 0.5u"W/m^2/K")
        wet_air_out = wet_air(u"K"(TAIR[i]); rh=RH[i], P_atmos=P_atmos)
        cp_air = wet_air_out.cp
        ρ_air = wet_air_out.ρ_air
        hd = (hc / (cp_air * ρ_air)) * (0.71 / 0.60)^0.666
        qevap, gwsurf = evap(tsurf=u"K"(T0[1]), tair=u"K"(TAIR[i]), rh=RH[i], rhsurf=100.0, hd=hd, elev=elev, pctwet=100.0, sat=false)
        λ_evap = get_λ_evap(T0[1])
        EP = qevap / λ_evap # evaporation potential, mm/s (kg/s)

        # run infiltration algorithm
        #T0 = u"K".([-14.947471488535426, -13.401168200078198, -12.482770331140387, -11.211433151563928, -10.799624779589632, -10.81365342874745, -11.104267771849331, -11.010673704296568, -10.134577542510105, 7.741666666666674]u"°C")
        #θ_soil0[1] = 1 - BD[1] / DD[1]
        # θ_soil = θ_soil0
        T10 = T0
        # EP = 3.544e-6u"kg/m^2/s"
        # ET = EP

        niter = ustrip(3600 / timestep) # number of interations for soil moisture calc
        for iter in 1:1
            infil_out = soil_water_balance(
                PE=PE,
                KS=KS,
                BB=BB,
                BD=BD,
                DD=DD,
                rh_loc=rh_loc,
                θ_soil=θ_soil0,
                ET=EP,
                T10=T0,
                depth=DEP,
                dt=timestep,
                elev=elev,
                L=L,
                rw=rw,
                pc=pc,
                rl=rl,
                sp=sp,
                r1=r1,
                lai=LAI,
                im=im,
                maxcount=maxcount
            )
            θ_soil0 = infil_out.θ_soil
        end
        step += 1
        T_soils[step, :] = T10
        sub = vcat(findall(isodd, 1:18), 18)
        θ_soils[step, :] = infil_out.θ_soil[sub]
        ψ_soils[step, :] = infil_out.ψ_soil[sub]
        rh_soils[step, :] = infil_out.rh_soil[sub]
    end
end
plot(1:nsteps, u"°C".(T_soils), xlabel="time", ylabel="soil temperature", lw=2, label=string.(DEP'))
plot!(1:nsteps, Matrix(soiltemps_NMR), xlabel="time", ylabel="soil temperature", lw=2, label = string.(DEP'), linestyle = :dash, linecolor="grey")
plot(1:nsteps, θ_soils, xlabel="time", ylabel="soil moisture (m^3/m^3)", lw=2, label = string.(DEP'))
plot(1:nsteps, ψ_soils, xlabel="time", ylabel="soil water potential", lw=2, label = string.(DEP'), ylim = (-1500, 0))
plot(1:nsteps, rh_soils, xlabel="time", ylabel="soil humidity", lw=2, label = string.(DEP'), ylim = (0, 1))
