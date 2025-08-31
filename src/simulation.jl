function runmicro(;
    lat=43.07305u"°", # latitude
    days=[15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349],
    ndmax=3,
    simdays=length(days), # can't be less than days or greater than days in year
    hours=collect(0.:1:24.), # hour of day for solrad
    depths=[0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm", # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
    refhyt=2u"m",
    elev=226.0u"m", # elevation (m)
    hori=fill(0.0u"°", 24), # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
    slope=0.0u"°", # slope (degrees, range 0-90)
    aspect=0.0u"°", # aspect (degrees, 0 = North, range 0-360)
    ruf=0.004u"m", # m roughness height
    zh=0u"m", # m heat transfer roughness height
    d0=0u"m", # zero plane displacement correction factor
    # soil properties
    # soil thermal parameters 
    λ_m=1.25u"W/m/K", # soil minerals thermal conductivity (W/mC)
    ρ_m=2.560u"Mg/m^3", # soil minerals density (Mg/m3)
    c_p_m=870.0u"J/kg/K", # soil minerals specific heat (J/kg-K)
    ρ_b_dry=2.56u"Mg/m^3", # dry soil bulk density (Mg/m3)
    θ_sat=0.26u"m^3/m^3", # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    # soil moisture model parameters
    PE=fill(0.7, 19)u"J/kg", #air entry potential J/kg
    KS=fill(0.0058, 19)u"kg*s/m^3", #saturated conductivity, kg s/m3
    BB=fill(1.7, 19), #soil 'b' parameter
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
    refls=fill(0.1, length(days)), # substrate solar reflectivity (decimal %)
    SHADES=fill(0.0, length(days)), # % shade cast by vegetation
    PCTWETS=fill(0.0, length(days)), # % surface wetness
    SLES=fill(0.96, length(days)), # - surface emissivity
    RAINdailys=([28, 28.2, 54.6, 79.7, 81.3, 100.1, 101.3, 102.5, 89.7, 62.4, 54.9, 41.2])u"kg/m^2",
    TMINN=([-14.3, -12.1, -5.1, 1.2, 6.9, 12.3, 15.2, 13.6, 8.9, 3, -3.2, -10.6] * 1.0)u"°C",
    TMAXX=([-3.2, 0.1, 6.8, 14.6, 21.3, 26.4, 29, 27.7, 23.3, 16.6, 7.8, -0.4] * 1.0)u"°C",
    WNMINN=([4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8] * 0.1)u"m/s",
    WNMAXX=([4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8] * 1.0)u"m/s",
    RHMINN=[50.2, 48.4, 48.7, 40.8, 40, 42.1, 45.5, 47.3, 47.6, 45, 51.3, 52.8],
    RHMAXX=[100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
    CCMINN=[50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1],
    CCMAXX=[50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1],
    TIMINS=[0, 0, 1, 1],
    TIMAXS=[1, 1, 0, 0],
    TAIRs=nothing,
    RHs=nothing,
    VELs=nothing,
    SOLRs=nothing,
    CLDs=nothing,
    RAINs=nothing,
    ZENhr=nothing,
    IRDhr=nothing,
    soilinit=u"K".((fill(7.741667, length(depths)))u"°C"),
    SoilMoist=[0.42, 0.42, 0.42, 0.43, 0.44, 0.44, 0.43, 0.42, 0.41, 0.42, 0.42, 0.43],
    LAIs=fill(0.1, length(days)),
    daily=false,
    runmoist=false,
    spinup=false,
    iuv=false, # this makes it take ages if true!
)

    ndays = length(days)
    tannul = mean(Unitful.ustrip.(vcat(TMAXX, TMINN)))u"°C" # annual mean temperature for getting monthly deep soil temperature (°C)
    tannulrun = fill(tannul, ndays) # monthly deep soil temperature (2m) (°C)
    #TODO - running mean when longer than a year

    # defining view factor based on horizon angles
    viewf = 1 - sum(sin.(hori)) / length(hori) # convert horizon angles to radians and calc view factor(s)

    # Soil properties
    # set up a profile of soil properites with depth for each day to be run
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
    ∑phase = zeros(Float64, numnodes_a)u"J" # zero phase transition for liquid water in soil

    # compute clear sky solar radiation
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
        τA,
    )
    # limit max zenith angles to 90°
    solrad_out.Zenith[solrad_out.Zenith.>90u"°"] .= 90u"°"
    solrad_out.ZenithSlope[solrad_out.ZenithSlope.>90u"°"] .= 90u"°"

    # vector for removing extra interpolated hour
    skip25 = setdiff(1:length(solrad_out.Zenith), 25:25:length(solrad_out.Zenith))

    if TAIRs === nothing
        # interpolate daily min/max forcing variables to hourly
        TAIRs25, VELs25, RHs25, CLDs25 = hourly_vars(
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
        RHs25[RHs25.>100] .= 100
        CLDs25[CLDs25.>100] .= 100
        TAIRs = TAIRs25[skip25] # remove every 25th output
        VELs = VELs25[skip25] # remove every 25th output
        RHs = RHs25[skip25] # remove every 25th output
        CLDs = CLDs25[skip25] # remove every 25th output
        SOLRs = solrad_out.Global[skip25] # remove every 25th output
        # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
        SOLRs = SOLRs .* (0.36 .+ 0.64 * (1.0 .- (CLDs / 100.0))) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
    end
    ZENRs = solrad_out.Zenith[skip25] # remove every 25th output
    ZSLs = solrad_out.ZenithSlope[skip25] # remove every 25th output

    # Initial conditions
    if !daily
        soilinit = fill(u"K"(mean(ustrip(TAIRs25[1:25]))u"°C"), numnodes_a) # mean monthly temp as in R version
    end
    soillayers = init_soillayers(numnodes_b)  # only once
    moistlayers = init_moistlayers(numnodes_b)  # only once
    θ_soil0_a = collect(fill(SoilMoist[1], numnodes_a)) # initial soil moisture
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
    nsteps = ndays * (length(hours) - 1)
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
    T0 = soilinit
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
    pool = 0.0u"kg/m^2"
    heights = [0.01] .* u"m"
    niter = ustrip(3600 / timestep)
    ∑phase = zeros(Float64, numnodes_a)u"J"
    for j in 1:ndays
        #j = 1
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
        # get today's weather
        sub1 = (iday*24-24+1):(iday*24)
        SOLR = SOLRs[sub1]
        ZENR = ZENRs[sub1]
        ZSL = ZSLs[sub1]
        TAIR = TAIRs[sub1]
        VEL = VELs[sub1]
        RH = RHs[sub1]
        CLD = CLDs[sub1]
        tspan = 0.0:60:(60*24)
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
            slep=slep, # check if this is what it should be - sle vs. slep (set as 1 in PAR in Fortran but then changed to user SLE later)
            pctwet=pctwet,
            nodes=nodes,
            tdeep=tdeep,
            θ_soil=θ_soil0_a,
            runmoist=runmoist
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
        step = 1
        # loop through hours of day
        if spinup && j == 1 && i == 1 || daily == false
            niter = ndmax # number of interations for steady periodic
        else
            niter = 1
        end
        if daily == false
            ∑phase = zeros(Float64, numnodes_a)u"J"
            sub2 = (iday*25-25+1):(iday*25) # for getting mean monthly over the 25 hrs as in fortran version
            soilinit = u"K"(mean(ustrip(TAIRs25[sub2]))u"°C") # make initial soil temps equal to mean annual temperature
            T0 = fill(soilinit, numnodes_a)
            #T_soils[step, :] = T0
            θ_soil0_a = collect(fill(SoilMoist[iday], numnodes_a)) # initial soil moisture
        end
        @inbounds for iter = 1:niter
            for i in 1:length(hours)
                if i < length(hours)
                    step = (j - 1) * (length(hours) - 1) + i
                end
                if i == 1 # make first hour of day equal last hour of previous iteration
                    T_soils[step, :] = T0
                    step = (j - 1) * (length(hours) - 1) + i
                    pool += rainfall
                    pools[step] = pool
                    pool = clamp(pool, 0.0u"kg/m^2", maxpool)
                    T_skys[step] = Tsky
                    λ_b, c_p_b, ρ_b = soil_properties(T0, θ_soil0_a, nodes, soilprops, elev, runmoist, false)
                    λ_bulk[step, :] = λ_b
                    c_p_bulk[step, :] = c_p_b
                    ρ_bulk[step, :] = ρ_b
                    if runmoist && i > 1
                        θ_soils[step, :] = infil_out.θ_soil[sub]
                        ψ_soils[step, :] = infil_out.ψ_soil[sub]
                        rh_soils[step, :] = infil_out.rh_soil[sub]
                    end
                else
                    # Parameters
                    params = MicroParams(;
                        soilprops,
                        dep=depths,
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
                    tspan = ((0.0 + (i - 2) * 60)u"minute", (60.0 + (i - 2) * 60)u"minute")  # 1 hour
                    prob = ODEProblem(soil_energy_balance!, T0, tspan, input)
                    sol = solve(prob, Tsit5(); saveat=60.0u"minute", reltol=1e-6u"K", abstol=1e-8u"K")
                    soiltemps = hcat(sol.u...)
                    # account for any phase transition of water in soil
                    T0 = soiltemps[:, 2]
                    if iter == niter # this makes it the same as the R version but really this should happen every time
                     ∑phase, qphase, T0 = phase_transition(soiltemps[:, 2], soiltemps[:, 1], ∑phase, θ_soil0_a, depths)
                    end
                    if i < length(hours)
                        T_soils[step, :] = T0
                    end
                    if runmoist
                        # compute scalar profiles
                        profile_out = get_profile(
                            refhyt=refhyt,
                            ruf=ruf,
                            zh=zh,
                            d0=d0,
                            TAREF=TAIRs[step],
                            VREF=VELs[step],
                            rh=RHs[step],
                            D0cm=u"°C"(T0[1]),  # top layer temp
                            ZEN=ZENRs[step],
                            heights=heights,
                            elev=elev,
                            warn=true
                        )

                        # convection
                        qconv = profile_out.QCONV

                        # evaporation
                        P_atmos = get_pressure(elev)
                        rh_loc = min(0.99, profile_out.RHs[2] / 100)
                        hc = max(abs(qconv / (T0[1] - u"K"(TAIRs[step]))), 0.5u"W/m^2/K")
                        wet_air_out = wet_air(u"K"(TAIRs[step]); rh=RHs[step], P_atmos=P_atmos)
                        c_p_air = wet_air_out.c_p
                        ρ_air = wet_air_out.ρ_air
                        hd = (hc / (c_p_air * ρ_air)) * (0.71 / 0.60)^0.666
                        qevap, gwsurf = evap(tsurf=u"K"(T0[1]), tair=u"K"(TAIRs[step]), rh=RHs[step], rhsurf=100.0, hd=hd, elev=elev, pctwet=pctwet, sat=true)
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
                    if i < length(hours)
                        pools[step] = pool
                    end
                    longwave_out = get_longwave(;
                        elev,
                        rh=RHs[step],
                        tair=TAIRs[step],
                        tsurf=T0[1],
                        slep,
                        sle,
                        cloud=CLDs[step],
                        viewf,
                        shade,
                    )
                    Tsky = longwave_out.Tsky
                    if i < length(hours)
                        T_skys[step] = Tsky
                    end
                    sub = vcat(findall(isodd, 1:numnodes_b), numnodes_b)
                    θ_soil0_a = θ_soil0_b[sub]
                    λ_b, c_p_b, ρ_b = soil_properties(T0, θ_soil0_a, nodes, soilprops, elev, runmoist, false)
                    λ_bulk[step, :] = λ_b
                    c_p_bulk[step, :] = c_p_b
                    ρ_bulk[step, :] = ρ_b
                    if runmoist
                        if i < length(hours)
                            θ_soils[step, :] = infil_out.θ_soil[sub]
                            ψ_soils[step, :] = infil_out.ψ_soil[sub]
                            rh_soils[step, :] = infil_out.rh_soil[sub]
                        end
                    end
                end
            end
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
        c_p_bulk=c_p_bulk,
        ρ_bulk=ρ_bulk,
        pools=pools
    )
end