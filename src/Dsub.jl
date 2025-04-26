using DifferentialEquations
using Interpolations
using Statistics
using Unitful
using Unitful: °
using UnitfulMoles
using Dates
using Plots
@compound H2O
@compound O2
@compound CO2
@compound N2

DEP = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm" # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
refhyt = 2u"m"
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
#days = [196]
hours = collect(0.:1:24.) # hour of day for solrad
lat = 43.1379° # latitude
elev = 0u"m" # elevation (m)
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

iuv = false
# soil properties# soil thermal parameters 
λ_m = 1.25u"W/m/K" # soil minerals thermal conductivity (W/mC)
ρ_m = 2.560u"Mg/m^3" # soil minerals density (Mg/m3)
cp_m = 870.0u"J/kg/K" # soil minerals specific heat (J/kg-K)
ρ_b_dry = 2.56u"Mg/m^3" # dry soil bulk density (Mg/m3)
θ_sat = 0.26u"m^3/m^3" # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
# Time varying environmental data
TIMINS = [0, 0, 1, 1] # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
TIMAXS = [1, 1, 0, 0] # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
TMINN = [-14.3, -12.1, -5.1, 1.2, 6.9, 12.3, 15.2, 13.6, 8.9, 3, -3.2, -10.6]u"°C" # minimum air temperatures (°C)
TMAXX = [-3.2, 0.1, 6.8, 14.6, 21.3, 26.4, 29, 27.7, 23.3, 16.6, 7.8, -0.4]u"°C" # maximum air temperatures (°C)
RHMINN = [50.2, 48.4, 48.7, 40.8, 40, 42.1, 45.5, 47.3, 47.6, 45, 51.3, 52.8] # min relative humidity (%)
RHMAXX = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0] # max relative humidity (%)
WNMINN = [4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8]u"m/s" # min wind speed (m/s)
WNMAXX = [4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8]u"m/s" # max wind speed (m/s)
CCMINN = [50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1] # min cloud cover (%)
CCMAXX = [50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1] # max cloud cover (%)
RAINFALL = ([28, 28.2, 54.6, 79.7, 81.3, 100.1, 101.3, 102.5, 89.7, 62.4, 54.9, 41.2]) / 1000u"m" # monthly mean rainfall (mm)
SoilMoist = fill(0.2, 10)
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
nodes = zeros(numnodes, ndays) # array of all possible soil nodes
nodes[1, 1:ndays] .= 10 # deepest node for first substrate type
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

# example application of soil_properties getting bulk properties from soil temperature, moisture and mineral properties
T_soil = collect(fill(u"K"(tannul), numnodes)) # for test of soil_properties
θ_soil = collect(fill(SoilMoist[1], numnodes)) # for test of soil_properties
node = nodes[:, 1]
λ_b, cp_b, ρ_b = soil_properties(T_soil, θ_soil, node, soilprops, numtyps, elev)

# compute solar radiation (need to make refl time varying)
solrad_out = solrad(days = days, hours = hours, lat = lat, elev = elev, hori = hori, slope = slope, aspect = aspect, refl = refl, iuv = iuv)
plot(solrad_out.Global)

# interpolate air temperature to hourly
TAIRs, WNs, RHs, CLDs = hourly_vars(TMINN=TMINN, TMAXX=TMAXX, WNMINN=WNMINN, WNMAXX=WNMAXX, RHMINN=RHMINN, RHMAXX=RHMAXX, CCMINN=CCMINN, CCMAXX=CCMAXX,solrad_out=solrad_out, TIMINS=TIMINS, TIMAXS=TIMAXS, daily=daily)
plot(TAIRs)
plot(WNs)
plot(RHs)
plot(CLDs)

# simulate a day
iday = 6
sub = (iday*25-24):(iday*25)
REFL = REFLS[iday]
SHADE = SHADES[iday] # daily shade (%)
SLE = SLES[iday] # set up vector of ground emissivities for each day
PCTWET = PCTWETS[iday] # set up vector of soil wetness for each day
tdeep = u"K"(tannulrun[iday]) # annual mean temperature for getting daily deep soil temperature (°C)
node = nodes[:, iday]

# get today's weather
SOLR = solrad_out.Global[sub]
ZENR = solrad_out.Zenith[sub]
ZSL = solrad_out.ZenithSlope[sub]
TAIR = TAIRs[sub]
VEL = WNs[sub]
RH = RHs[sub]
CLD = CLDs[sub]

# Initial conditions
soilinit = u"K"(mean(ustrip(TAIR))u"°C") # make initial soil temps equal to mean daily temperature
θ_soil = collect(fill(SoilMoist[iday], numnodes)) # initial soil moisture, constant for now

# create forcing weather variable splines
tspan = 0.:60:1440
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

plot(tspan, SOLRt(tspan))
plot(tspan, ZENRt(tspan))
plot(tspan, TAIRt(tspan))
plot(tspan, VELt(tspan))
plot(tspan, RHt(tspan))
plot(tspan, CLDt(tspan))

# Parameters
params = MicroParams(
    soilprops = soilprops,
    dep = DEP,
    refhyt = refhyt,
    ruf = ruf,
    d0 = d0,
    zh = zh,
    slope = slope,
    shade = shade,
    viewf = viewf,
    elev = elev,
    refl = refl,
    sle = sle,
    slep = sle, # check if this is what it should be - sle vs. slep (set as 1 in PAR in Fortran but then changed to user SLE later)
    pctwet = pctwet,
    tdeep = tdeep,
    θ_soil = θ_soil,
    runmoist = false,
    runsnow = false
)

# initial soil temperatures
T0 = fill(soilinit, numnodes)
T0[numnodes]=tdeep
T0[numnodes-2]=(u"K"(TMINN[iday])+u"K"(TMAXX[iday]))/2.0
T0[numnodes-1]=(T0[numnodes]+T0[numnodes-2])/2.0

tspan = (0.0u"minute", 1440.0u"minute")  # 1 hour
prob = ODEProblem(soil_energy_balance!, T0, tspan, params)
sol = solve(prob, Tsit5(); saveat=60.0u"minute")
soiltemps = hcat(sol.u...)
#plot(u"hr".(sol.t), u"°C".(soiltemps'), xlabel="Time", ylabel="Soil Temperature", lw=2)
T0 = soiltemps[:, 25] # new initial soil temps

# iterate through to get final soil temperature profile
niter = 2 # number of interations for steady periodic
for iter in 1:niter
    prob = ODEProblem(soil_energy_balance!, T0, tspan, params)
    sol = solve(prob, Tsit5(); saveat=60.0u"minute")
    soiltemps = hcat(sol.u...)
    #plot(u"hr".(sol.t), u"°C".(soiltemps'), xlabel="Time", ylabel="Soil Temperature", lw=2)
    T0 = soiltemps[:, 25] # new initial soil temps
end
labels = ["$(d)" for d in DEP]
plot(u"hr".(sol.t), u"°C".(soiltemps'), xlabel="Time", ylabel="Soil Temperature", lw=2, label = string.(DEP'))

# now get wind air temperature and humidity profiles
nsteps = 24
heights = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0] .* u"m"
profiles = Vector{NamedTuple}(undef, nsteps)  # or whatever you have from your loop
for i in 1:nsteps
    profiles[i] = get_profile(
        refhyt = refhyt,
        ruf = ruf,
        zh = zh,
        d0 = d0,
        TAREF = TAIRs[i],
        VREF = VEL[i],
        rh = RH[i],
        D0cm = u"°C"(soiltemps[1, i]),  # top layer temp at time i
        ZEN = ZENR[i],
        heights = heights,
        elev = elev
    )
end

VEL_matrix = hcat([getfield.(profiles, :VELs)[i] for i in 1:length(profiles)]...)    # velocities at all time steps and heights
TA_matrix   = hcat([getfield.(profiles, :TAs)[i]   for i in 1:length(profiles)]...)  # air temperatures at all time steps and heights
RH_matrix   = hcat([getfield.(profiles, :RHs)[i]   for i in 1:length(profiles)]...)       # relative humidity at all time steps and heights
Qconv_vec   = [getfield(p, :QCONV) for p in profiles]    # convective fluxes at all time steps
ustar_vec   = [getfield(p, :USTAR) for p in profiles]    # friction velocities
heights_vec = getfield(profiles[1], :heights)                     # constant across time
plot(u"hr".(sol.t[1:24]), VEL_matrix', xlabel="Time", ylabel="wind speed", lw=2, label = string.(last(heights_vec, 11)'))
plot(u"hr".(sol.t[1:24]), TA_matrix', xlabel="Time", ylabel="air temperature", lw=2, label = string.(last(heights_vec, 11)'))
plot(u"hr".(sol.t[1:24]), RH_matrix', xlabel="Time", ylabel="humidity (%)", lw=2, label = string.(last(heights_vec, 11)'))

numnodes
dT = zeros(Float64, numnodes)*u"K/minute"
N = numnodes
T = T0
t = 0.0
T_soil=T0