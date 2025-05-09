#using Microclimate
using Unitful
using Unitful: °, rad, R, kg, m
using Plots
using Statistics
using Interpolations
using DifferentialEquations
using CSV, DataFrames, Dates

# read in output from NicheMapR
soiltemps_NMR = (DataFrame(CSV.File("data/soil_FordDryLake.csv"))[:, 5:14]).*u"°C"
soilmoists_NMR = (DataFrame(CSV.File("data/soilmoist_FordDryLake.csv"))[:, 5:14])
#soilpots_NMR = (DataFrame(CSV.File("data/soilpot_FordDryLake.csv"))[:, 5:14])
#metout_NMR = DataFrame(CSV.File("data/metout_FordDryLake.csv"))

DEP = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm" # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
refhyt = 2u"m"
hours = collect(0.:1:24.) # hour of day for solrad
lat = 33.6547° # latitude
iuv = false # this makes it take ages if true!
elev = 120.0912u"m" # elevation (m)
hori = fill(0.0°, 24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
slope = 0.0° # slope (degrees, range 0-90)
aspect = 0.0° # aspect (degrees, 0 = North, range 0-360)
refl = 0.20 # substrate solar reflectivity (decimal %)
shade = 0.0 # % shade cast by vegetation
pctwet = 0.0 # % surface wetness
sle = 0.96 # - surface emissivity
ruf = 0.004u"m" # m roughness height
zh = 0u"m" # m heat transfer roughness height
d0 = 0u"m" # zero plane displacement correction factor

# soil properties

# soil thermal parameters 
λ_m = 2.5u"W/m/K" # soil minerals thermal conductivity (W/mC)
ρ_m = 2.560u"Mg/m^3" # soil minerals density (Mg/m3)
cp_m = 870.0u"J/kg/K" # soil minerals specific heat (J/kg-K)
ρ_b_dry = 1.3u"Mg/m^3" # dry soil bulk density (Mg/m3)
θ_sat = 0.26u"m^3/m^3" # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)

# soil moisture model parameters
runmoist = true
CampNormTbl9_1 = DataFrame(CSV.File("data/CampNormTbl9_1.csv"))
soiltype = 3 # 3 = sandy loam
PE = fill(CampNormTbl9_1[soiltype, 4], 19)u"J/kg" #air entry potential J/kg
KS = fill(CampNormTbl9_1[soiltype, 6], 19)u"kg*s/m^3" #saturated conductivity, kg s/m3
BB = fill(CampNormTbl9_1[soiltype, 5], 19) #soil 'b' parameter
BD = fill(ρ_b_dry, 19)u"Mg/m^3" # soil bulk density, Mg/m3
DD = fill(ρ_m, 19)u"Mg/m^3" # soil mineral density, Mg/m3
soiltype = 5 # change deeper nodes to 5 = a silt loam
PE[10:19] = fill(CampNormTbl9_1[soiltype, 4]u"J/kg", 10) #air entry potential J/kg
KS[10:19] = fill(CampNormTbl9_1[soiltype, 6]u"kg*s/m^3", 10) #saturated conductivity, kg s/m3
BB[10:19] = fill(CampNormTbl9_1[soiltype, 5], 10) #soil 'b' parameter
maxpool = 1.0e4u"kg/m^2"
L = [0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0]*1e4u"m/m^3" # root density at each node, mm/m3 (from Campell 1985 Soil Physics with Basic, p. 131)
rw =  2.5e+10u"m^3/kg/s" # resistance per unit length of root, m3 kg-1 s-1
pc = -1500.0u"J/kg" # critical leaf water potential for stomatal closure, J kg-1
rl = 2.0e6u"m^4/kg/s" # resistance per unit length of leaf, m3 kg-1 s-1
sp = 10.0 # stability parameter, -
r1 = 0.001u"m" # root radius, m
im = 1e-6u"kg/m^2/s" # maximum overall mass balance error allowed, kg
maxcount = 500
timestep = 360.0u"s"

τA = CSV.File("data/TAI_FordDryLake.csv").x

# Time varying environmental data
TAIRs = (CSV.File("data/TAIRhr_FordDryLake.csv").x)u"°C"
RHs = CSV.File("data/RHhr_FordDryLake.csv").x
VELs = (CSV.File("data/WNhr_FordDryLake.csv").x)u"m/s"
SOLRs = (CSV.File("data/SOLRhr_FordDryLake.csv").x)u"W/m^2"
CLDs = (CSV.File("data/CLDhr_FordDryLake.csv").x)
RAINs = (CSV.File("data/RAINhr_FordDryLake.csv").x/1000)u"m" * 1u"kg"/u"m^3"
RHs .= clamp.(RHs, 0, 100)
VELs .= clamp.(VELs, 0.1u"m/s", (Inf)u"m/s")
CLDs .= clamp.(CLDs, 0, 100)
RAINs .= clamp.(RAINs, 0u"kg/m^2", (Inf)u"kg/m^2")

days = collect(1.0:365.0)
SoilMoist = fill(0.2, length(days))
LAIs = fill(0.1, length(days))

daily = true
runmoist = true
# creating the arrays of environmental variables that are assumed not to change with month for this simulation
ndays = length(days)
SHADES = fill(shade, ndays) # daily shade (%)
SLES = fill(sle, ndays) # set up vector of ground emissivities for each day
REFLS = fill(refl, ndays) # set up vector of soil reflectances for each day
PCTWETS = fill(pctwet, ndays) # set up vector of soil wetness for each day
tannul = mean(Unitful.ustrip.(TAIRs))u"°C" # annual mean temperature for getting monthly deep soil temperature (°C)
tannul = 24.00356u"°C" 
tannulrun = fill(tannul, ndays) # monthly deep soil temperature (2m) (°C)

# defining view factor based on horizon angles
viewf = 1 - sum(sin.(hori)) / length(hori) # convert horizon angles to radians and calc view factor(s)

# Soil properties
# set up a profile of soil properites with depth for each day to be run
numtyps = 1 # number of soil types
numnodes_a = length(DEP) # number of soil nodes for temperature calcs and final output
numnodes_b = numnodes_a * 2 - 2 # number of soil nodes for soil moisture calcs
nodes_day = zeros(numnodes_a, ndays) # array of all possible soil nodes
nodes_day[1, 1:ndays] .= numnodes_a # deepest node for first substrate type
# Create an empty 10×5 matrix that can store any type (including different units)
soilprops = Matrix{Any}(undef, numnodes_a, 5)
# Fill row 1 (top layer) with the defined values
soilprops[1, 1] = ρ_b_dry
soilprops[1, 2] = θ_sat
soilprops[1, 3] = λ_m
soilprops[1, 4] = cp_m
soilprops[1, 5] = ρ_m
# Copy the same properties to all other layers (if desired)
for i in 2:numnodes_a
    soilprops[i, :] .= soilprops[1, :]
end

# compute solar radiation (need to make refl time varying)
solrad_out = solrad(
    days = days, 
    hours = hours, 
    lat = lat, 
    elev = elev, 
    hori = hori, 
    slope = slope, 
    aspect = aspect, 
    refl = refl, 
    iuv = iuv,
    τA = τA)
ZENRs = solrad_out.Zenith[Not(25:25:end)] # remove every 25th output
ZSLs = solrad_out.ZenithSlope[Not(25:25:end)] # remove every 25th output

# create forcing weather variable splines
tspan = 0.0:60:(60*24*(365)-60)
tmin = tspan .* u"minute"
interpSOLR = interpolate(SOLRs, BSpline(Cubic(Line(OnGrid()))))
interpZENR = interpolate(ZENRs, BSpline(Cubic(Line(OnGrid()))))
interpZSL = interpolate(ZSLs, BSpline(Cubic(Line(OnGrid()))))
interpTAIR = interpolate(u"K".(TAIRs), BSpline(Cubic(Line(OnGrid()))))
interpVEL = interpolate(VELs, BSpline(Cubic(Line(OnGrid()))))
interpRH = interpolate(RHs, BSpline(Cubic(Line(OnGrid()))))
interpCLD = interpolate(CLDs, BSpline(Cubic(Line(OnGrid()))))
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

# Initial conditions
soilinit = u"K"(tannul)#u"K"(mean(ustrip(TAIRs))u"°C") # make initial soil temps equal to mean daily temperature
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
pctwet = 0.0

days = collect(1:7)
# output arrays
nsteps = length(days) * (length(hours))
T_soils = Array{Float64}(undef, nsteps, numnodes_a)u"K"
θ_soils = Array{Float64}(undef, nsteps, numnodes_a)
ψ_soils = Array{Float64}(undef, nsteps, numnodes_a)u"J/kg"
rh_soils = Array{Float64}(undef, nsteps, numnodes_a)
pools = Array{Float64}(undef, nsteps)u"kg/m^2"
T_skys = Array{Float64}(undef, nsteps)u"K"

# initialise outputs
T_soils[1, :] = T0
θ_soils[1, :] = θ_soil0_a
sub = vcat(findall(isodd, 1:numnodes_b), numnodes_b)
θ_sat = 1.0 .- BD ./ DD
ψ_soils[1, :] =  PE[sub] .* (θ_sat[sub] ./ θ_soil0_a).^BB[sub]
MW = 0.01801528u"kg/mol" # molar mass of water, kg/mol
rh_soils[1, :] = exp.(MW .* ψ_soils[1, :] ./ (Unitful.R .* T0))
pools[1] = 0.0u"kg/m^2"
longwave_out = get_longwave(
    elev = elev, 
    rh = RHs[1], 
    tair = TAIRs[1], 
    tsurf = T0[1], 
    slep = sle, 
    sle = sle, 
    cloud = CLDs[1], 
    viewf = viewf
    )
Tsky = longwave_out.Tsky
T_skys[1] = Tsky

# simulate all days
step = 2
pool = 0.0u"kg/m^2"
heights = [0.01] .* u"m"
niter = ustrip(3600 / timestep)
sumphase = zeros(Float64, numnodes_a)u"J"
for j in 1:length(days)
    iday = j
    LAI = LAIs[iday]
    REFL = REFLS[iday]
    SHADE = SHADES[iday] # daily shade (%)
    SLE = SLES[iday] # set up vector of ground emissivities for each day
    pctwet = PCTWETS[iday] # set up vector of soil wetness for each day
    tdeep = u"K"(tannulrun[iday]) # annual mean temperature for getting daily deep soil temperature (°C)
    nodes = nodes_day[:, iday]

    # loop through hours of day
    for i in 1:length(hours)-1
        pool = clamp(pool + RAINs[step], 0u"kg/m^2", maxpool)
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
            slep=sle,
            pctwet=pctwet,
            nodes=nodes,
            tdeep=tdeep,
            θ_soil=θ_soil0_a
        )
        input = MicroInput(
            params,
            forcing
        )

        tspan = ((0.0 + (step - 2) * 60)u"minute", (60.0 + (step - 2) * 60)u"minute")  # 1 hour

        T0 = solve(ODEProblem(soil_energy_balance!, T0, tspan, input), Tsit5(); saveat=60.0u"minute").u[end]

        # if j == 1 && i == 1
        #     for _ in 1:2
        #         T0 = solve(ODEProblem(soil_energy_balance!, T0, tspan, input), Tsit5(); saveat=60.0u"minute").u[end]
        #     end
        # end
        # account for any phase transition of water in soil
        phase_transition!(T0, T_soils[step, :], sumphase, θ_soil0_a, DEP)

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
        cp_air = wet_air_out.cp
        ρ_air = wet_air_out.ρ_air
        hd = (hc / (cp_air * ρ_air)) * (0.71 / 0.60)^0.666
        qevap, gwsurf = evap(tsurf=u"K"(T0[1]), tair=u"K"(TAIRs[step-1]), rh=RHs[step-1], rhsurf=100.0, hd=hd, elev=elev, pctwet=pctwet, sat=true)
        λ_evap = get_λ_evap(T0[1])
        EP = max(1e-7u"kg/m^2/s", qevap / λ_evap) # evaporation potential, mm/s (kg/m2/s)

        # run infiltration algorithm
        if runmoist
            if pool > 0.0u"kg/m^2" # surface is wet - saturate it for infiltration
                θ_soil0_b[1] = 1 - BD[1] / DD[1]
            end
            infil_out = soil_water_balance(
                PE=PE,
                KS=KS,
                BB=BB,
                BD=BD,
                DD=DD,
                rh_loc=rh_loc,
                θ_soil=θ_soil0_b,
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
            θ_soil0_b = infil_out.θ_soil
            surf_evap = max(0.0u"kg/m^2", infil_out.evap)
            Δ_H2O = max(0.0u"kg/m^2", infil_out.Δ_H2O)
            pool = clamp(pool - Δ_H2O - surf_evap, 0u"kg/m^2", maxpool) # pooling surface water
            if pool > 0.0u"kg/m^2" # surface is wet - saturate it for infiltration
                θ_soil0_b[1] = 1 - BD[1] / DD[1]
            end
            for _ in 1:niter
                infil_out = soil_water_balance(
                    PE=PE,
                    KS=KS,
                    BB=BB,
                    BD=BD,
                    DD=DD,
                    rh_loc=rh_loc,
                    θ_soil=θ_soil0_b,
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
                θ_soil0_b = infil_out.θ_soil
                surf_evap = max(0.0u"kg/m^2", infil_out.evap)
                Δ_H2O = max(0.0u"kg/m^2", infil_out.Δ_H2O)
                pool = clamp(pool - Δ_H2O - surf_evap, 0u"kg/m^2", maxpool)
                if pool > 0.0u"kg/m^2"
                    θ_soil0_b[1] = 1 - BD[1] / DD[1]
                end
            end
        end

        T_soils[step, :] = T0
        pctwet = clamp(abs(surf_evap / (EP * timestep) * 100), 0, 100)
        pools[step] = pool

        longwave_out = get_longwave(
            elev=elev,
            rh=RHs[step-1],
            tair=TAIRs[step-1],
            tsurf=T0[1],
            slep=sle,
            sle=sle,
            cloud=CLDs[step-1],
            viewf=viewf
        )
        Tsky = longwave_out.Tsky

        sub = vcat(findall(isodd, 1:numnodes_b), numnodes_b)
        θ_soil0_a = θ_soil0_b[sub]
        θ_soils[step, :] = infil_out.θ_soil[sub]
        ψ_soils[step, :] = infil_out.ψ_soil[sub]
        rh_soils[step, :] = infil_out.rh_soil[sub]
        T_skys[step] = Tsky
        step += 1
    end
end
nsteps = 24*length(days)
pstart = 1
pfinish = 7*24
plot(pstart:pfinish, u"°C".(T_soils[pstart:pfinish, :]), xlabel="time", ylabel="soil temperature", lw=2, label=string.(DEP'), legend = :none, ylim = (-10u"°C", 85u"°C"))
plot(pstart:pfinish, Matrix(soiltemps_NMR[pstart:pfinish, :]), xlabel="time", ylabel="soil temperature", lw=2, label = string.(DEP'), legend = :none, ylim = (-10u"°C", 85u"°C"))
plot(pstart:pfinish, θ_soils[pstart:pfinish, :], xlabel="time", ylabel="soil moisture (m^3/m^3)", lw=2, label = string.(DEP'), legend = :none, ylim = (0, 0.5))
plot(pstart:pfinish, Matrix(soilmoists_NMR[pstart:pfinish, :]), xlabel="time", ylabel="soil moisture", legend = :none, lw=2, label = string.(DEP'), linestyle = :dash, ylim = (0, 0.5))


plot(1:nsteps, u"°C".(T_soils[1:nsteps, :]), xlabel="time", ylabel="soil temperature", lw=2, label=string.(DEP'), legend = :none, ylim = (-10u"°C", 85u"°C"))
plot(1:nsteps, Matrix(soiltemps_NMR[1:nsteps, :]), xlabel="time", ylabel="soil temperature", lw=2, label = string.(DEP'), legend = :none, ylim = (-10u"°C", 85u"°C"))
plot(1:nsteps, θ_soils[1:nsteps, :], xlabel="time", ylabel="soil moisture (m^3/m^3)", lw=2, label = string.(DEP'), legend = :none, ylim = (0, 0.5))
plot(1:nsteps, Matrix(soilmoists_NMR[1:nsteps, :]), xlabel="time", ylabel="soil moisture", lw=2, label = string.(DEP'), legend = :none, ylim = (0, 0.5))
# plot(1:nsteps, ψ_soils[1:nsteps, :], xlabel="time", ylabel="soil water potential", lw=2, label = string.(DEP'), legend = :none, ylim = (-3e5, 0))
# plot!(1:nsteps, (Matrix(soilpots_NMR[1:nsteps, :]))u"J/kg", xlabel="time", ylabel="soil water potential", lw=2, label = string.(DEP'), linestyle = :dash, linecolor="grey")
# plot(1:nsteps, rh_soils[1:nsteps, :], xlabel="time", ylabel="soil humidity", lw=2, label = string.(DEP'), ylim = (0, 1))
#plot(1:nsteps, RAINs[1:nsteps], xlabel="time", ylabel="rainfall", lw=2, linecolor="blue")
plot(1:nsteps, pools[1:nsteps], xlabel="time", ylabel="pooling", lw=2, linecolor="blue")

# tspan = 0.0:60:(60*24*(365)-60)

# plot(1:nsteps, TAIRt(tspan)[1:nsteps])
# plot!(1:nsteps, u"K".((metout_NMR.TAREF[1:nsteps])u"°C"), col = "red")
# plot(1:nsteps, SOLRt(tspan)[1:nsteps])
# plot!(1:nsteps, (metout_NMR.SOLR[1:nsteps])u"W/m^2", col = "red")
# plot(1:nsteps, ZENRt(tspan)[1:nsteps])
# plot!(1:nsteps, (metout_NMR.ZEN[1:nsteps])u"°", col = "red")
# plot(1:nsteps, VELt(tspan)[1:nsteps])
# plot!(1:nsteps, (metout_NMR.VREF[1:nsteps])u"m/s", col = "red")
# plot(1:nsteps, VELt(tspan)[1:nsteps])
# plot(1:nsteps, u"°C".(T_skys[1:nsteps]))
# plot!(1:nsteps, ((metout_NMR.TSKYC[1:nsteps])u"°C"), col = "red")

