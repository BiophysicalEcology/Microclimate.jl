function runmicro(;
    # locations, times, depths and heights
    latitude=43.07305u"°", # latitude
    days=[15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349], # days of year to simulate - TODO leap years
    hours=collect(0.:1:24.), # hour of day for solrad
    reference_height=2u"m", # reference height of weather data (air temperature, wind speed, humidity)
    depths=[0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 100.0, 200.0]u"cm", # soil nodes - keep spacing close near the surface
    heights=[1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 100.0]u"cm", # air nodes for temperature, wind speed and humidity profile
    # terrain
    elevation=226.0u"m", # elevation (m)
    horizon_angles=fill(0.0u"°", 24), # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
    slope=0.0u"°", # slope (degrees, range 0-90)
    aspect=0.0u"°", # aspect (degrees, 0 = North, range 0-360)
    ruf=0.004u"m", # m roughness height
    zh=0u"m", # m heat transfer roughness height
    d0=0u"m", # zero plane displacement correction factor
    # soil thermal parameters 
    λ_m=1.25u"W/m/K", # soil minerals thermal conductivity (W/mC)
    ρ_m=2.560u"Mg/m^3", # soil minerals density (Mg/m3)
    c_p_m=870.0u"J/kg/K", # soil minerals specific heat (J/kg-K)
    ρ_b_dry=2.56u"Mg/m^3", # dry soil bulk density (Mg/m3)
    θ_sat=0.26u"m^3/m^3", # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    # soil moisture model soil parameters
    PE=fill(0.7, length(depths) * 2 - 2)u"J/kg", #air entry potential J/kg
    KS=fill(0.0058, length(depths) * 2 - 2)u"kg*s/m^3", #saturated conductivity, kg s/m3
    BB=fill(1.7, length(depths) * 2 - 2), #soil 'b' parameter
    BD=fill(ρ_b_dry, length(depths) * 2 - 2)u"Mg/m^3", # soil bulk density, Mg/m3
    DD=fill(ρ_m, length(depths) * 2 - 2)u"Mg/m^3", # soil mineral density, Mg/m3
    # soil moisture plant parameters
    L=[0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4, 0.4, 0, 0] * 1e4u"m/m^3", # root density at each node, mm/m3 (from Campell 1985 Soil Physics with Basic, p. 131)
    rw=2.5e+10u"m^3/kg/s", # resistance per unit length of root, m3 kg-1 s-1
    pc=-1500.0u"J/kg", # critical leaf water potential for stomatal closure, J kg-1
    rl=2.0e6u"m^4/kg/s", # resistance per unit length of leaf, m3 kg-1 s-1
    sp=10.0, # stability parameter, -
    r1=0.001u"m", # root radius, m
    # soil moisture simulation controls
    im=1e-6u"kg/m^2/s", # maximum overall mass balance error allowed, kg
    maxcount=500, # maximum iterations of soil moisture algorithm
    timestep=360.0u"s", # time step over which to simulate soil moisture (<= 1 hour)
    maxpool=1.0e4u"kg/m^2", # maximum depth of pooling water
    τA=DEFAULT_τA,
    # daily environmental vectors
    albedos=fill(0.1, length(days)), # substrate albedo (decimal %)
    shades=fill(0.0, length(days)), # % shade cast by vegetation
    pctwets=fill(0.0, length(days)), # % surface wetness
    sles=fill(0.96, length(days)), # - surface emissivity
    daily_rainfall=([28, 28.2, 54.6, 79.7, 81.3, 100.1, 101.3, 102.5, 89.7, 62.4, 54.9, 41.2])u"kg/m^2",
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
    # hourly weather vectors
    TAIRs=nothing,
    RHs=nothing,
    VELs=nothing,
    SOLRs=nothing,
    CLDs=nothing,
    RAINs=nothing,
    ZENhr=nothing,
    IRDhr=nothing,
    # intial conditions
    soilinit=u"K".((fill(7.741667, length(depths)))u"°C"),
    SoilMoist=[0.42, 0.42, 0.42, 0.43, 0.44, 0.44, 0.43, 0.42, 0.41, 0.42, 0.42, 0.43],
    LAIs=fill(0.1, length(days)),
    ndmax=3, # number of iterations per day
    daily=false, # doing consecutive days?
    runmoist=false, # run soil moisture algorithm?
    spinup=false, # spin-up the first day by ndmax iterations?
    iuv=false, # this makes it take ages if true!
)

    ndays = length(days)
    tannul = mean(Unitful.ustrip.(vcat(TMAXX, TMINN)))u"°C" # annual mean temperature for getting monthly deep soil temperature (°C)
    #TODO - running mean when longer than a year
    tannulrun = fill(tannul, ndays) # monthly deep soil temperature (2m) (°C)

    # defining view factor based on horizon angles
    viewfactor = 1 - sum(sin.(horizon_angles)) / length(horizon_angles) # convert horizon angles to radians and calc view factor(s)

    # Soil properties
    # set up a profile of soil properites with depth for each day to be run
    numnodes_a = length(depths) # number of soil nodes for temperature calcs and final output
    numnodes_b = numnodes_a * 2 - 2 # number of soil nodes for soil moisture calcs
    # if runmoist & length(depths) != 10
    #     depths_orig = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 100.0, 200.0]u"cm"
    #     mids_orig = [(depths_orig[i] + depths_orig[i+1]) / 2 for i in 1:length(depths_orig)-1]
    #     expanded_depths_orig = sort(unique([depths_orig; mids_orig]))
    #     mids = [(depths[i] + depths[i+1]) / 2 for i in 1:length(depths)-1]
    #     expanded_depths = sort(unique([depths; mids]))
    #     function spline_to_depths(base_depths, base_values, new_depths)
    #         # strip units for interpolation
    #         x = ustrip.(u"cm", base_depths)
    #         y = base_values ./ unit(base_values)  # plain numbers
    #         itp = LinearInterpolation(x, y, extrapolation_bc=Line())  # or CubicSplineInterpolation
    #         # evaluate at new depths
    #         y_new = itp.(ustrip.(u"cm", new_depths)) * unit(base_values)
    #         return y_new
    #     end
    #     L2 = spline_to_depths(expanded_depths, L, expanded_depths)
    # end
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
        latitude,
        elevation,
        horizon_angles,
        slope,
        aspect,
        albedos,
        iuv,
        τA,
    )
    # limit max zenith angles to 90°
    solrad_out.Zenith[solrad_out.Zenith.>90u"°"] .= 90u"°"
    solrad_out.ZenithSlope[solrad_out.ZenithSlope.>90u"°"] .= 90u"°"

    # vector for removing extra interpolated hour
    skip25 = setdiff(1:length(solrad_out.Zenith), 25:25:length(solrad_out.Zenith))

    if isnothing(TAIRs)
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
    end
    ZENRs = solrad_out.Zenith[skip25] # remove every 25th output
    ZSLs = solrad_out.ZenithSlope[skip25] # remove every 25th output
    # adjust for cloud using Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
    doy  = repeat(days, inner=length(hours))[skip25]
    DIRs = solrad_out.Direct[skip25] # remove every 25th output
    DIFs = solrad_out.Scattered[skip25] # remove every 25th output
    global_solar, diffuse_solar, direct_solar = cloud_adjust_radiation(CLDs / 100., DIFs, DIRs, ZENRs, doy)
    if isnothing(SOLRs)
        SOLRs = global_solar
    else
        global_solar = SOLRs # keep global solar output as original solar radiation input
    end


    # Initial conditions
    if !daily && !isnothing(soilinit)
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
    λ_b, c_p_b, ρ_b = soil_properties(T0, θ_soil0_a, nodes_day[:, 1], soilprops, elevation, runmoist, false)
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
        elevation=elevation,
        rh=RHs[1],
        tair=TAIRs[1],
        tsurf=T0[1],
        slep=sles[1],
        sle=sles[1],
        cloud=CLDs[1],
        viewfactor=viewfactor,
        shade=shades[1]
    )
    Tsky = longwave_out.Tsky
    T_skys[1] = Tsky
    # simulate all days
    pool = 0.0u"kg/m^2"
    heights_water_balance = [0.01] .* u"m"
    niter_moist = ustrip(3600 / timestep)
    ∑phase = zeros(Float64, numnodes_a)u"J"
    infil_out = nothing
    for j in 1:ndays
        #j = 1
        iday = j
        lai = LAIs[iday]
        albedo = albedos[iday]
        shade = shades[iday] # daily shade (%)
        sle = sles[iday] # set up vector of ground emissivities for each day
        slep = sle # - cloud emissivity
        pctwet = pctwets[iday] # set up vector of soil wetness for each day
        tdeep = u"K"(tannulrun[iday]) # annual mean temperature for getting daily deep soil temperature (°C)
        nodes = nodes_day[:, iday]
        rainfall = daily_rainfall[iday]
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
        params = MicroParams(;
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
            slep, # check if this is what it should be - sle vs. slep (set as 1 in PAR in Fortran but then changed to user SLE later)
            pctwet,
            nodes,
            tdeep,
            θ_soil=θ_soil0_a,
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
        input = MicroInputs(;
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
                    if runmoist
                        infil_out, pctwet, pool, θ_soil0_b = get_soil_water_balance(;
                                reference_height,
                                ruf,
                                zh,
                                d0,
                                TAIRs,
                                VELs,
                                RHs,
                                ZENRs,
                                T0,
                                heights=heights_water_balance,
                                elevation,
                                pool,
                                θ_soil0_b,
                                PE,
                                KS,
                                BB,
                                BD,
                                DD,
                                depths,
                                timestep,
                                L,
                                rw,
                                pc,
                                rl,
                                sp,
                                r1,
                                lai,
                                im,
                                maxcount,
                                moistlayers,
                                niter_moist,
                                pctwet,
                                step,
                                maxpool,
                        )
                    end
                    pools[step] = pool
                    pool = clamp(pool, 0.0u"kg/m^2", maxpool)
                    T_skys[step] = Tsky
                    λ_b, c_p_b, ρ_b = soil_properties(T0, θ_soil0_a, nodes, soilprops, elevation, runmoist, false)
                    λ_bulk[step, :] = λ_b
                    c_p_bulk[step, :] = c_p_b
                    ρ_bulk[step, :] = ρ_b
                    if runmoist && iday > 1
                        θ_soils[step, :] = infil_out.θ_soil[sub]
                        ψ_soils[step, :] = infil_out.ψ_soil[sub]
                        rh_soils[step, :] = infil_out.rh_soil[sub]
                    end
                else
                    # Parameters
                    params = MicroParams(;
                        soilprops,
                        depths=depths,
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
                        infil_out, pctwet, pool, θ_soil0_b = get_soil_water_balance(;
                                reference_height,
                                ruf,
                                zh,
                                d0,
                                TAIRs,
                                VELs,
                                RHs,
                                ZENRs,
                                T0,
                                heights = heights_water_balance,
                                elevation,
                                pool,
                                θ_soil0_b,
                                PE,
                                KS,
                                BB,
                                BD,
                                DD,
                                depths,
                                timestep,
                                L,
                                rw,
                                pc,
                                rl,
                                sp,
                                r1,
                                lai,
                                im,
                                maxcount,
                                moistlayers,
                                niter_moist,
                                pctwet,
                                step,
                                maxpool,
                        )
                    end
                    if i < length(hours)
                        pools[step] = pool
                    end
                    longwave_out = get_longwave(;
                        elevation,
                        rh=RHs[step],
                        tair=TAIRs[step],
                        tsurf=T0[1],
                        slep,
                        sle,
                        cloud=CLDs[step],
                        viewfactor,
                        shade,
                    )
                    Tsky = longwave_out.Tsky
                    if i < length(hours)
                        T_skys[step] = Tsky
                    end
                    sub = vcat(findall(isodd, 1:numnodes_b), numnodes_b)
                    θ_soil0_a = θ_soil0_b[sub]
                    λ_b, c_p_b, ρ_b = soil_properties(T0, θ_soil0_a, nodes, soilprops, elevation, runmoist, false)
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
    # compute air temperature, wind speed and relative humidity profiles
    profile_out = map(1:length(TAIRs)) do i
        # compute scalar profiles
        get_profile(
            TAREF=TAIRs[i],
            VREF=VELs[i],
            rh=RHs[i],
            D0cm=u"°C"(T_soils[i, 1]),  # top layer temp
            ZEN=ZENRs[i],
            heights=heights,
            elevation=elevation,
            warn=true
        )
    end
    flip2vectors(x) = (; (k => getfield.(x, k) for k in keys(x[1]))...)
    profiles = flip2vectors(profile_out); # pull out each output as a vector
    air_temperature = reduce(hcat, profiles.TAs)'
    wind_speed = reduce(hcat, profiles.VELs)'
    relative_humidity = reduce(hcat, profiles.RHs)'

    return (;
        air_temperature,
        wind_speed,
        relative_humidity,
        cloud_cover=CLDs,
        global_solar,
        direct_solar,
        diffuse_solar,
        zenith_angle=ZENRs,
        sky_temperature=T_skys,
        soil_temperature=T_soils,
        soil_moisture=θ_soils,
        soil_waterpotential=ψ_soils,
        soil_humidity=rh_soils,
        soil_thermalconductivity=λ_bulk,
        soil_specificheat=c_p_bulk,
        soil_bulkdensity=ρ_bulk,
        surface_water=pools,
        solrad_out=solrad_out,
        profile_out=profile_out,
    )
end