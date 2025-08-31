using Microclimate
using Unitful, UnitfulMoles
using Plots
using Statistics
using Interpolations
using OrdinaryDiffEq
using CSV, DataFrames, Dates
using Test

# read in output from NicheMapR
soiltemps_NMR = (DataFrame(CSV.File("test/data/soil_FordDryLake.csv"))[:, 5:14]) .* u"°C"
soilmoists_NMR = (DataFrame(CSV.File("test/data/soilmoist_FordDryLake.csv"))[:, 5:14])
soilconds_NMR = (DataFrame(CSV.File("test/data/tcond_FordDryLake.csv"))[:, 5:14])
date = DateTime(2015, 1, 1):Hour(1):DateTime(2015, 12, 31, 23)

microinput_vec = DataFrame(CSV.File("test/data/init_daily/microinput.csv"))[:, 2]

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

# Time varying environmental data
TAIRs = Float64.(CSV.File("test/data/init_daily/TAIRhr.csv").x)u"°C"
RHs = Float64.(CSV.File("test/data/init_daily/RHhr.csv").x)
VELs = Float64.(CSV.File("test/data/init_daily/WNhr.csv").x)u"m/s"
SOLRs = Float64.(CSV.File("test/data/init_daily/SOLRhr.csv").x)u"W/m^2"
CLDs = Float64.(CSV.File("test/data/init_daily/CLDhr.csv").x)
RAINs = Float64.(CSV.File("test/data/init_daily/RAINhr.csv").x)u"kg" / u"m^2"
RHs .= clamp.(RHs, 0, 100)
VELs .= clamp.(VELs, 0.1u"m/s", (Inf)u"m/s")
CLDs .= clamp.(CLDs, 0, 100)
RAINs .= clamp.(RAINs, 0u"kg/m^2", (Inf)u"kg/m^2")

refhyt = microinput[:Refhyt] * 1.0u"m"
depths = ((DataFrame(CSV.File("test/data/init_daily/DEP.csv"))[:, 2]) / 100.0)u"m" # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
days = collect(1:Int(length(TAIRs) / 24)) # days of year to run (for solrad)
hours = collect(0.:1:24.) # hour of day for solrad
lat = (microinput[:ALAT] + microinput[:AMINUT] / 60) * 1.0u"°" # latitude
iuv = Bool(Int(microinput[:IUV])) # this makes it take ages if true!
ndmax = (microinput[:ndmax]) # number of iterations per day
elev = microinput[:ALTT] * 1.0u"m" # elevation (m)
hori = (DataFrame(CSV.File("test/data/init_daily/hori.csv"))[:, 2]) * 1.0u"°" # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
slope = microinput[:slope] * 1.0u"°" # slope (degrees, range 0-90)
aspect = microinput[:azmuth] * 1.0u"°" # aspect (degrees, 0 = North, range 0-360)
ruf = microinput[:RUF] * 1.0u"m" # m roughness height
zh = microinput[:ZH] * 1.0u"m" # m heat transfer roughness height
d0 = microinput[:D0] * 1.0u"m" # zero plane displacement correction factor
# soil properties# soil thermal parameters 
ρ_b_dry = (CSV.File("test/data/init_daily/soilprop.csv")[1, 1][2]) * 1.0u"Mg/m^3" # dry soil bulk density (Mg/m3)
θ_sat = (CSV.File("test/data/init_daily/soilprop.csv")[1, 1][3]) * 1.0u"m^3/m^3" # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
λ_m = (CSV.File("test/data/init_daily/soilprop.csv")[1, 1][4]) * 1.0u"W/m/K" # soil minerals thermal conductivity (W/mC)
c_p_m = (CSV.File("test/data/init_daily/soilprop.csv")[1, 1][5]) * 1.0u"J/kg/K" # soil minerals specific heat (J/kg-K)
ρ_m = (CSV.File("test/data/init_daily/soilprop.csv")[1, 1][6]) * 1.0u"Mg/m^3" # soil minerals density (Mg/m3)
# Time varying environmental data
TMINN = (DataFrame(CSV.File("test/data/init_daily/TMINN.csv"))[:, 2] * 1.0)u"°C" # minimum air temperatures (°C)
TMAXX = (DataFrame(CSV.File("test/data/init_daily/TMAXX.csv"))[:, 2] * 1.0)u"°C" # maximum air temperatures (°C)
RAINFALL = ((DataFrame(CSV.File("test/data/init_daily/rain.csv"))[:, 2] * 1.0) / 1000)u"m" # monthly mean rainfall (mm)
SoilMoist = (DataFrame(CSV.File("test/data/init_daily/moists.csv"))[:, 2:13] .* 1.0)#fill(0.0, ndays)
daily = Bool(Int(microinput[:microdaily]))
runmoist = Bool(Int(microinput[:runmoist]))
spinup = Bool(Int(microinput[:spinup]))

# creating the arrays of environmental variables that are assumed not to change with month for this simulation
ndays = length(days)
SHADES = (DataFrame(CSV.File("test/data/init_daily/MINSHADES.csv"))[:, 2] * 1.0) # daily shade (%)
SLES = (DataFrame(CSV.File("test/data/init_daily/SLES.csv"))[:, 2] * 1.0) # set up vector of ground emissivities for each day
refls = (DataFrame(CSV.File("test/data/init_daily/REFLS.csv"))[:, 2] * 1.0) # set up vector of soil reflectances for each day (decimal %)
PCTWETS = (DataFrame(CSV.File("test/data/init_daily/PCTWET.csv"))[:, 2] * 1.0) # set up vector of soil wetness for each day (%)
tannul = mean(Unitful.ustrip.(TAIRs))u"°C" # annual mean temperature for getting monthly deep soil temperature (°C)
tannulrun = fill(tannul, ndays) # monthly deep soil temperature (2m) (°C)

# defining view factor based on horizon angles
viewf = 1 - sum(sin.(hori)) / length(hori) # convert horizon angles to radians and calc view factor(s)

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
soilinit = (DataFrame(CSV.File("test/data/init_daily/soilinit.csv"))[1:numnodes, 2] * 1.0)u"°C" # set up vector of soil wetness for each day (%)
∑phase = zeros(Float64, numnodes)u"J"

# soil properties

# soil moisture model parameters
PE = (DataFrame(CSV.File("test/data/init_daily/PE.csv"))[:, 2] * 1.0u"J/kg") # set up vector of ground emissivities for each day
KS = (DataFrame(CSV.File("test/data/init_daily/KS.csv"))[:, 2] * 1.0u"kg*s/m^3") # set up vector of ground emissivities for each day
BB = (DataFrame(CSV.File("test/data/init_daily/BB.csv"))[:, 2] * 1.0) # set up vector of ground emissivities for each day
BD = (DataFrame(CSV.File("test/data/init_daily/BD.csv"))[:, 2] * 1.0u"Mg/m^3") # set up vector of ground emissivities for each day
DD = (DataFrame(CSV.File("test/data/init_daily/DD.csv"))[:, 2] * 1.0u"Mg/m^3") # set up vector of ground emissivities for each day

maxpool = microinput[:maxpool] * 1000.0u"kg/m^2" # max depth for water pooling on the surface, mm (to account for runoff)
L = DataFrame(CSV.File("test/data/init_daily/L.csv"))[:, 2]u"m/m^3" # root density at each node, mm/m3 (from Campell 1985 Soil Physics with Basic, p. 131) # max depth for water pooling on the surface, mm (to account for runoff)
rw = microinput[:RW]u"m^3/kg/s" # resistance per unit length of root, m3 kg-1 s-1
pc = -microinput[:PC]u"J/kg" # critical leaf water potential for stomatal closure, J kg-1
rl = microinput[:RL]u"m^4/kg/s" # resistance per unit length of leaf, m3 kg-1 s-1
sp = microinput[:SP] # stability parameter, -
r1 = microinput[:R1]u"m" # root radius, m
im = microinput[:IM]u"kg/m^2/s" # maximum overall mass balance error allowed, kg

maxcount = microinput[:MAXCOUNT]
timestep = microinput[:moiststep]u"s"

τA = CSV.File("test/data/init_daily/TAI.csv").x

raindf = DataFrame(date=date, rainfall=RAINs)
raindf.day = Date.(raindf.date)
raindfdaily = combine(groupby(raindf, :day), :rainfall => sum => :daily_rainfall)
RAINdailys = raindfdaily.daily_rainfall
SoilMoist = Matrix((DataFrame(CSV.File("test/data/init_daily/moists.csv"))[:, 2:(ndays-1)] .* 1.0))
LAIs = (DataFrame(CSV.File("test/data/init_daily/LAI.csv"))[:, 2] * 1.0u"Mg/m^3")

# Soil properties
# set up a profile of soil properites with depth for each day to be run
numtyps = 1 # number of soil types
numnodes_a = length(depths) # number of soil nodes for temperature calcs and final output
numnodes_b = numnodes_a * 2 - 2 # number of soil nodes for soil moisture calcs
nodes_day = zeros(numnodes_a, ndays) # array of all possible soil nodes
nodes_day[1, 1:ndays] .= numnodes_a # deepest node for first substrate type
# Create an empty 10×5 matrix that can store any type (including different units)
soilprops = Matrix{Any}(undef, numnodes_a, 5)
# Fill row 1 (top layer) with the defined values
soilprops[1, 1] = ρ_b_dry
soilprops[1, 2] = θ_sat
soilprops[1, 3] = λ_m
soilprops[1, 4] = c_p_m
soilprops[1, 5] = ρ_m
# Copy the same properties to all other layers (if desired)
for i in 2:numnodes_a
    soilprops[i, :] .= soilprops[1, :]
end

# compute solar radiation (need to make refl time varying)
solrad_out = solrad(;
    days,
    hours,
    lat,
    elev,
    hori,
    slope,
    aspect,
    refls,
    iuv,
    τA)
# only need zenith angles in this case
solrad_out.Zenith[solrad_out.Zenith.>90u"°"] .= 90u"°"
solrad_out.ZenithSlope[solrad_out.ZenithSlope.>90u"°"] .= 90u"°"

skip25 = setdiff(1:length(solrad_out.Zenith), 25:25:length(solrad_out.Zenith))
ZENRs = solrad_out.Zenith[skip25] # remove every 25th output
ZSLs = solrad_out.ZenithSlope[skip25] # remove every 25th output

# create forcing weather variable splines
tspan = 0.0:60:(60*24*ndays)
tmin = tspan .* u"minute"
interpSOLR = interpolate([SOLRs; SOLRs[end]], BSpline(Linear()))
interpZENR = interpolate([ZENRs; ZENRs[end]], BSpline(Linear()))
interpZSL = interpolate([ZSLs; ZSLs[end]], BSpline(Linear()))
interpTAIR = interpolate(u"K".([TAIRs; TAIRs[end]]), BSpline(Linear()))
interpVEL = interpolate([VELs; VELs[end]], BSpline(Linear()))
interpRH = interpolate([RHs; RHs[end]], BSpline(Linear()))
interpCLD = interpolate([CLDs; CLDs[end]], BSpline(Linear()))
SOLRt = scale(interpSOLR, tspan)
ZENRt = scale(interpZENR, tspan)
ZSLt = scale(interpZSL, tspan)
TAIRt = scale(interpTAIR, tspan)
VELt = scale(interpVEL, tspan)
RHt = scale(interpRH, tspan)
CLDt = scale(interpCLD, tspan)
forcing = MicroForcing(;
    SOLRt,
    ZENRt,
    ZSLt,
    TAIRt,
    VELt,
    RHt,
    CLDt,
)

soillayers = init_soillayers(numnodes_b)  # only once
moistlayers = init_moistlayers(numnodes_b)  # only once

# Initial conditions
soilinit = u"K"(mean(ustrip(TAIRs))u"°C") # make initial soil temps equal to mean daily temperature
T0 = fill(soilinit, numnodes_a)
θ_soil0_a = collect(fill(SoilMoist[1], numnodes_a)) # initial soil moisture
# intitial soil moisture
θ_soil0_b = similar(θ_soil0_a, numnodes_b)  # preallocate vector of length numnodes_b
jj = 1
for ii in 1:numnodes_b
    if isodd(ii)
        θ_soil0_b[ii] = θ_soil0_a[jj]
        jj += 1
    else
        θ_soil0_b[ii] = θ_soil0_b[ii-1]
    end
end

# output arrays
nsteps = ndays * (length(hours)-1)
T_soils = Array{Float64}(undef, nsteps, numnodes_a)u"K"
θ_soils = Array{Float64}(undef, nsteps, numnodes_a)
ψ_soils = Array{Float64}(undef, nsteps, numnodes_a)u"J/kg"
rh_soils = Array{Float64}(undef, nsteps, numnodes_a)
λ_bulk = Array{Float64}(undef, nsteps, numnodes_a)u"W/m/K"
c_p_bulk = Array{Float64}(undef, nsteps, numnodes_a)u"J/kg/K"
ρ_bulk = Array{Float64}(undef, nsteps, numnodes_a)u"kg/m^3"
pools = Array{Float64}(undef, nsteps)u"kg/m^2"
T_skys = Array{Float64}(undef, nsteps)u"K"

# initialise outputs
T_soils[1, :] = T0
θ_soils[1, :] = θ_soil0_a
λ_b, c_p_b, ρ_b = soil_properties(T0, θ_soil0_a, nodes_day[:, 1], soilprops, elev, runmoist, false)
λ_bulk[1, :] = λ_b
c_p_bulk[1, :] = c_p_b
ρ_bulk[1, :] = ρ_b
sub = vcat(findall(isodd, 1:numnodes_b), numnodes_b)
θ_sat = 1.0 .- BD ./ DD
ψ_soils[1, :] = PE[sub] .* (θ_sat[sub] ./ θ_soil0_a) .^ BB[sub]
MW = 0.01801528u"kg/mol" # molar mass of water, kg/mol
rh_soils[1, :] = clamp.(exp.(MW .* ψ_soils[1, :] ./ (Unitful.R .* T0)), 0, 1)
pools[1] = 0.0u"kg/m^2"
# sky temperature
longwave_out = get_longwave(
    elev=elev,
    rh=RHs[1],
    tair=TAIRs[1],
    tsurf=T0[1],
    slep=SLES[1],
    sle=SLES[1],
    cloud=CLDs[1],
    viewf=viewf,
    shade=SHADES[1]
)
Tsky = longwave_out.Tsky
T_skys[1] = Tsky
# simulate all days
step = 2
pool = 0.0u"kg/m^2"
heights = [0.01] .* u"m"
niter = ustrip(3600 / timestep)
∑phase = zeros(Float64, numnodes_a)u"J"
for j in 1:ndays
    iday = j
    lai = LAIs[iday]
    refl = refls[iday]
    shade = SHADES[iday] # daily shade (%)
    sle = SLES[iday] # set up vector of ground emissivities for each day
    slep = sle # - cloud emissivity
    pctwet = PCTWETS[iday] # set up vector of soil wetness for each day
    tdeep = u"K"(tannulrun[iday]) # annual mean temperature for getting daily deep soil temperature (°C)
    nodes = nodes_day[:, iday]
    rainfall = RAINdailys[iday]
    # loop through hours of day
    if step == 2
        istart = 2
    else
        istart = 1
    end
    for i in istart:length(hours)-1
        if i == 1
            pool += rainfall#RAINs[step-1]
        end
        pool = clamp(pool, 0.0u"kg/m^2", maxpool)
        # Parameters
        dep = depths
        params = MicroParams(;
            soilprops,
            dep,
            refhyt,
            ruf,
            d0,
            zh,
            slope,
            shade,
            viewf,
            elev,
            refl,
            sle,
            slep,
            pctwet,
            nodes,
            tdeep,
            θ_soil=θ_soil0_a,
            runmoist
        )
        input = MicroInputs(
            params,
            forcing,
            soillayers
        )
        tspan = ((0.0 + (step - 2) * 60)u"minute", (60.0 + (step - 2) * 60)u"minute")  # 1 hour
        prob = ODEProblem(soil_energy_balance!, T0, tspan, input)
        sol = solve(prob, Tsit5(); saveat=60.0u"minute", reltol=1e-6u"K", abstol=1e-8u"K")
        soiltemps = hcat(sol.u...)
        # account for any phase transition of water in soil
        ∑phase, qphase, T0 = phase_transition(soiltemps[:, 2], soiltemps[:, 1], ∑phase, θ_soil0_a, depths)
        if spinup && j == 1 && i == 1
            for _ in 1:2
                prob = ODEProblem(soil_energy_balance!, T0, tspan, input)
                sol = solve(prob, Tsit5(); saveat=60.0u"minute", reltol=1e-6u"K", abstol=1e-8u"K")
                soiltemps = hcat(sol.u...)
                # account for any phase transition of water in soil
                ∑phase, qphase, T0 = phase_transition(soiltemps[:, 2], soiltemps[:, 1], ∑phase, θ_soil0_a, depths)
            end
        end
        T_soils[step, :] = T0

        if runmoist
            # compute scalar profiles
            profile_out = get_profile(
                refhyt=refhyt,
                ruf=ruf,
                zh=zh,
                d0=d0,
                TAREF=TAIRs[step-1],
                VREF=VELs[step-1],
                rh=RHs[step-1],
                D0cm=u"°C"(T0[1]),  # top layer temp at time i
                ZEN=ZENRs[step-1],
                heights=heights,
                elev=elev,
                warn=true
            )

            # convection
            qconv = profile_out.QCONV

            # evaporation
            P_atmos = get_pressure(elev)
            rh_loc = min(0.99, profile_out.RHs[2] / 100)
            hc = max(abs(qconv / (T0[1] - u"K"(TAIRs[step-1]))), 0.5u"W/m^2/K")
            wet_air_out = wet_air(u"K"(TAIRs[step-1]); rh=RHs[step-1], P_atmos=P_atmos)
            c_p_air = wet_air_out.c_p
            ρ_air = wet_air_out.ρ_air
            hd = (hc / (c_p_air * ρ_air)) * (0.71 / 0.60)^0.666
            qevap, gwsurf = evap(tsurf=u"K"(T0[1]), tair=u"K"(TAIRs[step-1]), rh=RHs[step-1], rhsurf=100.0, hd=hd, elev=elev, pctwet=pctwet, sat=true)
            λ_evap = get_λ_evap(T0[1])
            EP = max(1e-7u"kg/m^2/s", qevap / λ_evap) # evaporation potential, mm/s (kg/m2/s)

            if pool > 0.0u"kg/m^2" # surface is wet - saturate it for infiltration
                θ_soil0_b[1] = 1 - BD[1] / DD[1]
            end
            # run infiltration algorithm
            infil_out = soil_water_balance(;
                PE,
                KS,
                BB,
                BD,
                DD,
                rh_loc,
                θ_soil=θ_soil0_b,
                ET=EP,
                T10=T0,
                depth=depths,
                dt=timestep,
                elev,
                L,
                rw,
                pc,
                rl,
                sp,
                r1,
                lai,
                im,
                maxcount,
                ml=moistlayers
            )
            θ_soil0_b = infil_out.θ_soil
            surf_evap = max(0.0u"kg/m^2", infil_out.evap)
            Δ_H2O = max(0.0u"kg/m^2", infil_out.Δ_H2O)
            pool = clamp(pool - Δ_H2O - surf_evap, 0.0u"kg/m^2", maxpool) # pooling surface water
            if pool > 0.0u"kg/m^2" # surface is wet - saturate it for infiltration
                θ_soil0_b[1] = 1 - BD[1] / DD[1]
            end
            for _ in 1:(niter-1)
                infil_out = soil_water_balance(;
                    PE,
                    KS,
                    BB,
                    BD,
                    DD,
                    rh_loc=rh_loc,
                    θ_soil=θ_soil0_b,
                    ET=EP,
                    T10=T0,
                    depth=depths,
                    dt=timestep,
                    elev=elev,
                    L,
                    rw,
                    pc,
                    rl,
                    sp,
                    r1,
                    lai,
                    im,
                    maxcount,
                    ml=moistlayers
                )
                θ_soil0_b = infil_out.θ_soil
                surf_evap = max(0.0u"kg/m^2", infil_out.evap)
                Δ_H2O = max(0.0u"kg/m^2", infil_out.Δ_H2O)
                pool = clamp(pool - Δ_H2O - surf_evap, 0.0u"kg/m^2", maxpool)
                if pool > 0.0u"kg/m^2"
                    θ_soil0_b[1] = 1 - BD[1] / DD[1]
                end
            end
        end
        if runmoist
            pctwet = clamp(abs(surf_evap / (EP * timestep) * 100), 0, 100)
        end
        pools[step] = pool
        longwave_out = get_longwave(;
            elev,
            rh=RHs[step-1],
            tair=TAIRs[step-1],
            tsurf=T0[1],
            slep,
            sle,
            cloud=CLDs[step-1],
            viewf,
            shade,
        )
        Tsky = longwave_out.Tsky
        T_skys[step] = Tsky
        sub = vcat(findall(isodd, 1:numnodes_b), numnodes_b)
        θ_soil0_a = θ_soil0_b[sub]
        λ_b, c_p_b, ρ_b = soil_properties(T0, θ_soil0_a, nodes, soilprops, elev, runmoist, false)
        λ_bulk[step, :] = λ_b
        c_p_bulk[step, :] = c_p_b
        ρ_bulk[step, :] = ρ_b
        if runmoist
            θ_soils[step, :] = infil_out.θ_soil[sub]
            ψ_soils[step, :] = infil_out.ψ_soil[sub]
            rh_soils[step, :] = infil_out.rh_soil[sub]
        end
        step += 1
    end
end

pstart = 1
pfinish = ndays*24
pstart = 24*180
pfinish = 24*187
plot(pstart:pfinish, u"°C".(T_soils[pstart:pfinish, :]), xlabel="time", ylabel="soil temperature", lw=2, label=string.(depths'), legend=:none, ylim=(-10u"°C", 85u"°C"))
plot!(pstart:pfinish, Matrix(soiltemps_NMR[pstart:pfinish, :]), xlabel="time", ylabel="soil temperature", lw=2, label=string.(depths'), legend=:none, ylim=(-10u"°C", 85u"°C"), linestyle=:dash, linecolor="grey")
plot(pstart:pfinish, θ_soils[pstart:pfinish, :], xlabel="time", ylabel="soil moisture (m^3/m^3)", lw=2, label = string.(depths'), legend = :none, ylim = (0, 0.5))
plot!(pstart:pfinish, Matrix(soilmoists_NMR[pstart:pfinish, :]), xlabel="time", ylabel="soil moisture", legend = :none, lw=2, label = string.(depths'), ylim = (0, 0.5), linestyle=:dash, linecolor="grey")

# now try the simulation function
micro_out = runmicro(;
    refhyt,
    depths,
    days,
    hours,
    lat,
    iuv,
    ndmax,
    elev,
    hori,
    slope,
    aspect,
    ruf,
    zh,
    d0,
    ρ_b_dry,
    θ_sat = 0.26, # overwritten to be a vector above, so manually inputting
    λ_m,
    c_p_m,
    ρ_m,
    TMINN,
    TMAXX,
    RAINdailys = ustrip(RAINFALL)u"kg/m^2",
    SoilMoist,
    daily,
    runmoist,
    spinup,
    TAIRs,
    RHs,
    VELs,
    SOLRs,
    CLDs,
    RAINs,
    SHADES,
    SLES,
    refls,
    PCTWETS,
    PE,
    KS,
    BB,
    BD,
    DD,
    maxpool,
    L,
    rw,
    pc,
    rl,
    sp,
    r1,
    im,
    maxcount,
    timestep,
    τA,
    soilinit = fill(soilinit, numnodes_a),
    LAIs,
)

plot(pstart:pfinish, u"°C".(micro_out.T_soils[pstart:pfinish, :]), xlabel="time", ylabel="soil temperature", lw=2, label=string.(depths'), legend=:none, ylim=(-10u"°C", 85u"°C"))
plot!(pstart:pfinish, Matrix(soiltemps_NMR[pstart:pfinish, :]), xlabel="time", ylabel="soil temperature", lw=2, label=string.(depths'), legend=:none, ylim=(-10u"°C", 85u"°C"), linestyle=:dash, linecolor="grey")

plot(pstart:pfinish, micro_out.θ_soils[pstart:pfinish, :], xlabel="time", ylabel="soil moisture (m^3/m^3)", lw=2, label = string.(depths'), legend = :none, ylim = (0, 0.5))
plot!(pstart:pfinish, Matrix(soilmoists_NMR[pstart:pfinish, :]), xlabel="time", ylabel="soil moisture", legend = :none, lw=2, label = string.(depths'), ylim = (0, 0.5), linestyle=:dash, linecolor="grey")
