function runmicro(;
    depths=[0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm", # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
    refhyt=2u"m",
    hours=collect(0.:1:24.), # hour of day for solrad
    lat=33.6547u"°", # latitude
    elev=120.0912u"m", # elevation (m)
    hori=fill(0.0u"°", 24), # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
    slope=0.0u"°", # slope (degrees, range 0-90)
    aspect=0.0u"°", # aspect (degrees, 0 = North, range 0-360)
    refl=0.20, # substrate solar reflectivity (decimal %)
    shade=0.0, # % shade cast by vegetation
    pctwet=0.0, # % surface wetness
    sle=0.96, # - surface emissivity
    ruf=0.004u"m", # m roughness height
    zh=0u"m", # m heat transfer roughness height
    d0=0u"m", # zero plane displacement correction factor
    # soil properties
    # soil thermal parameters 
    λ_m=2.5u"W/m/K", # soil minerals thermal conductivity (W/mC)
    ρ_m=2.560u"Mg/m^3", # soil minerals density (Mg/m3)
    cp_m=870.0u"J/kg/K", # soil minerals specific heat (J/kg-K)
    ρ_b_dry=1.3u"Mg/m^3", # dry soil bulk density (Mg/m3)
    θ_sat=0.26u"m^3/m^3", # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    # soil moisture model parameters
    PE=fill(2.1, 19)u"J/kg", #air entry potential J/kg
    KS=fill(1.9e-4, 19)u"kg*s/m^3", #saturated conductivity, kg s/m3
    BB=fill(4.7, 19), #soil 'b' parameter
    BD=fill(ρ_b_dry, 19)u"Mg/m^3", # soil bulk density, Mg/m3
    DD=fill(ρ_m, 19)u"Mg/m^3", # soil mineral density, Mg/m3
    maxpool=1.0e4u"kg/m^2",
    L=[0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4, 0.4, 0, 0] * 1e4u"m/m^3", # root density at each node, mm/m3 (from Campell 1985 Soil Physics with Basic, p. 131)
    rw=2.5e+10u"m^3/kg/s", # resistance per unit length of root, m3 kg-1 s-1
    pc=-1500.0u"J/kg", # critical leaf water potential for stomatal closure, J kg-1
    rl=2.0e6u"m^4/kg/s", # resistance per unit length of leaf, m3 kg-1 s-1
    sp=10.0, # stability parameter, -
    r1=0.001u"m", # root radius, m
    im=1e-6u"kg/m^2/s", # maximum overall mass balance error allowed, kg
    maxcount=500,
    timestep=360.0u"s",
    τA=[0.269904738, 0.266147825, 0.262442906, 0.258789404, 0.255186744, 0.251634356, 0.248131676, 0.2412732,
        0.234606887, 0.228128378, 0.221833385, 0.215717692, 0.20977715, 0.204007681, 0.198405272, 0.187685927,
        0.177588357, 0.168082846, 0.159140695, 0.150734206, 0.142836655, 0.135422274, 0.128466227, 0.12194459,
        0.115834329, 0.110113284, 0.104760141, 0.099754417, 0.09507644, 0.090707328, 0.086628967, 0.082823998,
        0.07927579, 0.075968428, 0.072886691, 0.070016034, 0.067342571, 0.064853053, 0.062534858, 0.060375964,
        0.058364941, 0.056490925, 0.054743609, 0.053113222, 0.051590514, 0.050166738, 0.046408775, 0.045302803,
        0.044259051, 0.043271471, 0.042334415, 0.041442618, 0.040591184, 0.039775572, 0.038991583, 0.038235345,
        0.037503301, 0.036792197, 0.036099067, 0.034101935, 0.033456388, 0.032817888, 0.032184949, 0.031556287,
        0.030930816, 0.030307633, 0.029065372, 0.027825562, 0.027205981, 0.026586556, 0.025967391, 0.025348692,
        0.024114005, 0.023498886, 0.021669152, 0.021066668, 0.019292088, 0.018144698, 0.016762709, 0.015451481,
        0.014949794, 0.014224263, 0.013093462, 0.012670686, 0.012070223, 0.011164062, 0.010241734, 0.009731103,
        0.009507687, 0.009212683, 0.008965785, 0.008827751, 0.008710756, 0.008574128, 0.008462605, 0.008446967,
        0.008539475, 0.009015237, 0.009748444, 0.010586023, 0.011359647, 0.011901268, 0.012062153, 0.011735443,
        0.010882215, 0.009561062, 0.007961182, 0.006438984, 0.005558204, 0.006133532, 0.009277754
    ],
    # Time varying environmental data
    TAIRs=TAIRs,
    RHs=RHs,
    VELs=VELs,
    SOLRs=SOLRs,
    CLDs=CLDs,
    RAINs=RAINs,
    days=collect(1.0:365.0),
    SoilMoist=fill(0.2, length(days)),
    LAIs=fill(0.1, length(days)),
    daily=true,
    runmoist=true,
    spinup=false,
    iuv=false # this makes it take ages if true!
)
    # creating the arrays of environmental variables that are assumed not to change with month for this simulation
    ndays = length(days)
    SHADES = fill(shade, ndays) # daily shade (%)
    SLES = fill(sle, ndays) # set up vector of ground emissivities for each day
    REFLS = fill(refl, ndays) # set up vector of soil reflectances for each day
    PCTWETS = fill(pctwet, ndays) # set up vector of soil wetness for each day
    tannul = mean(Unitful.ustrip.(TAIRs))u"°C" # annual mean temperature for getting monthly deep soil temperature (°C)
    tannulrun = fill(tannul, ndays) # monthly deep soil temperature (2m) (°C)

    # defining view factor based on horizon angles
    viewf = 1 - sum(sin.(hori)) / length(hori) # convert horizon angles to radians and calc view factor(s)

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
    soilprops[1, 4] = cp_m
    soilprops[1, 5] = ρ_m
    # Copy the same properties to all other layers (if desired)
    for i in 2:numnodes_a
        soilprops[i, :] .= soilprops[1, :]
    end

    # compute solar radiation (need to make refl time varying)
    solrad_out = solrad(
        days=days,
        hours=hours,
        lat=lat,
        elev=elev,
        hori=hori,
        slope=slope,
        aspect=aspect,
        refl=refl,
        iuv=iuv,
        τA=τA)
        skip25 = setdiff(1:length(solrad_out.Zenith), 25:25:length(solrad_out.Zenith))
    ZENRs = solrad_out.Zenith[skip25] # remove every 25th output
    ZSLs = solrad_out.ZenithSlope[skip25] # remove every 25th output

    # create forcing weather variable splines
    tspan = 0.0:60:(60*24*(ndays))
    tmin = tspan .* u"minute"
    SOLR_ext = [SOLRs; SOLRs[end]]
    interpSOLR = interpolate([SOLRs; SOLRs[end]], BSpline(Cubic(Line(OnGrid()))))
    interpZENR = interpolate([ZENRs; ZENRs[end]], BSpline(Cubic(Line(OnGrid()))))
    interpZSL = interpolate([ZSLs; ZSLs[end]], BSpline(Cubic(Line(OnGrid()))))
    interpTAIR = interpolate(u"K".([TAIRs; TAIRs[end]]), BSpline(Cubic(Line(OnGrid()))))
    interpVEL = interpolate([VELs; VELs[end]], BSpline(Cubic(Line(OnGrid()))))
    interpRH = interpolate([RHs; RHs[end]], BSpline(Cubic(Line(OnGrid()))))
    interpCLD = interpolate([CLDs; CLDs[end]], BSpline(Cubic(Line(OnGrid()))))
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
    soillayers = init_soillayers(numnodes_b)  # only once
    moistlayers = init_moistlayers(numnodes_b)  # only once
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

    days = collect(1:ndays)
    # output arrays
    nsteps = length(days) * (length(hours))
    T_soils = Array{Float64}(undef, nsteps + 1, numnodes_a)u"K"
    θ_soils = Array{Float64}(undef, nsteps + 1, numnodes_a)
    ψ_soils = Array{Float64}(undef, nsteps + 1, numnodes_a)u"J/kg"
    rh_soils = Array{Float64}(undef, nsteps + 1, numnodes_a)
    λ_bulk = Array{Float64}(undef, nsteps + 1, numnodes_a)u"W/m/K"
    cp_bulk = Array{Float64}(undef, nsteps + 1, numnodes_a)u"J/kg/K"
    ρ_bulk = Array{Float64}(undef, nsteps + 1, numnodes_a)u"kg/m^3"
    pools = Array{Float64}(undef, nsteps + 1)u"kg/m^2"
    T_skys = Array{Float64}(undef, nsteps + 1)u"K"

    # initialise outputs
    T_soils[1, :] = T0
    θ_soils[1, :] = θ_soil0_a
    λ_b, cp_b, ρ_b = soil_properties(T0, θ_soil0_a, nodes_day[:, 1], soilprops, elev, runmoist, false)
    λ_bulk[1, :] = λ_b
    cp_bulk[1, :] = cp_b
    ρ_bulk[1, :] = ρ_b
    sub = vcat(findall(isodd, 1:numnodes_b), numnodes_b)
    θ_sat = 1.0 .- BD ./ DD
    ψ_soils[1, :] = PE[sub] .* (θ_sat[sub] ./ θ_soil0_a) .^ BB[sub]
    MW = 0.01801528u"kg/mol" # molar mass of water, kg/mol
    rh_soils[1, :] = exp.(MW .* ψ_soils[1, :] ./ (Unitful.R .* T0))
    pools[1] = 0.0u"kg/m^2"
    longwave_out = get_longwave(
        elev=elev,
        rh=RHs[1],
        tair=TAIRs[1],
        tsurf=T0[1],
        slep=sle,
        sle=sle,
        cloud=CLDs[1],
        viewf=viewf,
        shade=shade
    )
    Tsky = longwave_out.Tsky
    T_skys[1] = Tsky

    # simulate all days
    step = 2
    pool = 0.0u"kg/m^2"
    heights = [0.01] .* u"m"
    niter = ustrip(3600 / timestep)
    ∑phase = zeros(Float64, numnodes_a)u"J"
    for j in 1:length(days)
        iday = j
        lai = LAIs[iday]
        refl = REFLS[iday]
        shade = SHADES[iday] # daily shade (%)
        sle = SLES[iday] # set up vector of ground emissivities for each day
        pctwet = PCTWETS[iday] # set up vector of soil wetness for each day
        tdeep = u"K"(tannulrun[iday]) # annual mean temperature for getting daily deep soil temperature (°C)
        nodes = nodes_day[:, iday]
        #rainfall = RAINdailys[iday]
        # loop through hours of day
        for i in 1:length(hours)-1
            if i == 1
                #pool += rainfall#RAINs[step-1]
                pool += RAINs[step-1]
           end
            pool = clamp(pool, 0.0u"kg/m^2", maxpool)
            # Parameters
            params = MicroParams(
                soilprops=soilprops,
                dep=depths,
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
            input = MicroInputs(
                params,
                forcing,
                soillayers
            )

            tspan = ((0.0 + (step - 2) * 60)u"minute", (60.0 + (step - 2) * 60)u"minute")  # 1 hour
            T0 = solve(ODEProblem(soil_energy_balance!, T0, tspan, input), Tsit5(); saveat=60.0u"minute").u[end]

            if spinup && j == 1 && i == 1
                for _ in 1:2
                    T0 = solve(ODEProblem(soil_energy_balance!, T0, tspan, input), Tsit5(); saveat=60.0u"minute").u[end]
                end
            end
            # account for any phase transition of water in soil
            phase_transition!(T0, T_soils[step-1, :], ∑phase, θ_soil0_a, depths)

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
                    depth=depths,
                    dt=timestep,
                    elev=elev,
                    L=L,
                    rw=rw,
                    pc=pc,
                    rl=rl,
                    sp=sp,
                    r1=r1,
                    lai=lai,
                    im=im,
                    maxcount=maxcount,
                    ml=moistlayers
                )
                θ_soil0_b = infil_out.θ_soil
                surf_evap = max(0.0u"kg/m^2", infil_out.evap)
                Δ_H2O = max(0.0u"kg/m^2", infil_out.Δ_H2O)
                pool = clamp(pool - Δ_H2O - surf_evap, 0.0u"kg/m^2", maxpool) # pooling surface water
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
                        depth=depths,
                        dt=timestep,
                        elev=elev,
                        L=L,
                        rw=rw,
                        pc=pc,
                        rl=rl,
                        sp=sp,
                        r1=r1,
                        lai=lai,
                        im=im,
                        maxcount=maxcount,
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

            T_soils[step, :] = T0
            if runmoist
                pctwet = clamp(abs(surf_evap / (EP * timestep) * 100), 0, 100)
            end
            pools[step] = pool

            longwave_out = get_longwave(
                elev=elev,
                rh=RHs[step-1],
                tair=TAIRs[step-1],
                tsurf=T0[1],
                slep=sle,
                sle=sle,
                cloud=CLDs[step-1],
                viewf=viewf,
                shade=shade
            )
            Tsky = longwave_out.Tsky
            T_skys[step] = Tsky
            sub = vcat(findall(isodd, 1:numnodes_b), numnodes_b)
            θ_soil0_a = θ_soil0_b[sub]
            λ_b, cp_b, ρ_b = soil_properties(T0, θ_soil0_a, nodes, soilprops, elev, runmoist, false)
            λ_bulk[step, :] = λ_b
            cp_bulk[step, :] = cp_b
            ρ_bulk[step, :] = ρ_b
            if runmoist
                θ_soils[step, :] = infil_out.θ_soil[sub]
                ψ_soils[step, :] = infil_out.ψ_soil[sub]
                rh_soils[step, :] = infil_out.rh_soil[sub]
            end
            step += 1
        end
    end
    return (
        solrad_out=solrad_out,
        T_skys=T_skys,
        T_soils=T_soils,
        θ_soils=θ_soils,
        ψ_soils=ψ_soils,
        rh_soils=rh_soils,
        λ_bulk=λ_bulk,
        cp_bulk=cp_bulk,
        ρ_bulk=ρ_bulk,
        pools=pools
    )
end