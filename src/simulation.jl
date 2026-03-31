"""
    MicroProblem(; latitude, solar_terrain, micro_terrain, soil_moisture_model, soil_thermal_model, environment_daily, ...)

Microclimate simulation problem specification.

# Key formulation choices (type-dispatched)
- `soil_moisture_model`: `CampbellSoilHydraulics(; ..., mode)` where `mode` is
  `PrescribedSoilMoisture()` (default) or `DynamicSoilMoisture()`
- `time_mode`: `MonthlyRepresentativeMode()` (default) or `HourlyMode(; spinup_first_day=false)`
- `convergence`: `FixedSoilTemperatureIterations(3)` (default) or
  `SoilTemperatureConvergenceTolerance(; tolerance, max_iterations_per_day)`
- `diffuse_fraction_model`: `ErbsDiffuseFraction()` (default)
- `vapour_pressure_equation`: `GoffGratch()` (default), `Teten()`, or `Huang()`
"""
@kwdef struct MicroProblem{D,H,Dep,Ht,Lat,SM,ST,MT,SMM,STM,EMM,EH,ED,IST,ISM,VP,SOS,SOK,TM,Conv,DFM}
    # locations, times, depths and heights
    days::D = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349] # days of year to simulate - TODO leap years - why not use real dates?
    hours::H = collect(0.0:1:23.0) # hour of day for solar_radiation
    depths::Dep = DEFAULT_DEPTHS # soil nodes - keep spacing close near the surface
    heights::Ht = [0.01, 2]u"m" # air nodes for temperature, wind speed and humidity profile, last height is reference height for weather data
    latitude::Lat
    # Physical model objects
    solar_model::SM = SolarProblem()
    solar_terrain::ST
    micro_terrain::MT
    soil_moisture_model::SMM # CampbellSoilHydraulics with mode field
    soil_thermal_model::STM
    environment_minmax::EMM
    environment_hourly::EH = nothing
    environment_daily::ED
    # Initial conditions
    initial_soil_temperature::IST = fill(u"K"(7.741667u"°C"), length(depths))
    initial_soil_moisture::ISM = fill(0.42 * 0.25, length(depths))
    # Formulation choices
    vapour_pressure_equation::VP = GoffGratch() # GoffGratch(), Teten(), or Huang()
    # ODE solver for soil temperature integration. Any SciML algorithm works, e.g.:
    #   Tsit5()           — adaptive 5th-order Runge-Kutta (default, accurate)
    #   RK4()             — fixed-step 4th-order RK; set soil_ode_kwargs=(; dt=6u"minute", adaptive=false)
    #   Euler()           — fixed-step Euler; fastest but least accurate
    soil_ode_solver::SOS = Tsit5()
    soil_ode_kwargs::SOK = (; reltol=1e-6u"K", abstol=1e-8u"K")
    time_mode::TM = MonthlyRepresentativeMode() # MonthlyRepresentativeMode() or HourlyMode()
    convergence::Conv = FixedSoilTemperatureIterations(3) # FixedSoilTemperatureIterations or SoilTemperatureConvergenceTolerance
    diffuse_fraction_model::DFM = ErbsDiffuseFraction()
    hourly_rainfall::Bool = false # use hourly rainfall?
end

function example_microclimate_problem(;
    latitude = 43.07305u"°",
    micro_terrain=default_terrain(),
    soil_moisture_model=example_soil_hydraulics(),
    soil_thermal_model=example_soil_thermal_parameters(),
    environment_minmax=example_monthly_weather(),
    environment_daily=example_daily_environment(),
    kw...
)
    MicroProblem(;
         latitude, micro_terrain, solar_terrain, soil_moisture_model, soil_thermal_model, environment_minmax, environment_daily, kw...
    )
end

# TODO some of these can be actual defaults ?
function example_micro_terrain(;
    elevation = 226.0u"m", # elevation (m)
    roughness_height = 0.004u"m", # heat transfer roughness height
    karman_constant = 0.4, # Kármán constant
    dyer_constant = 16.0, # coefficient from Dyer and Hicks for Φ_m (momentum), γ
)
    MicroTerrain(; elevation, roughness_height, karman_constant, dyer_constant)
end

function example_monthly_weather(;
    reference_temperature_min = [-14.3, -12.1, -5.1, 1.2, 6.9, 12.3, 15.2, 13.6, 8.9, 3, -3.2, -10.6]u"°C",
    reference_temperature_max = [-3.2, 0.1, 6.8, 14.6, 21.3, 26.4, 29, 27.7, 23.3, 16.6, 7.8, -0.4]u"°C",
    reference_wind_min = [4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8] * 0.1u"m/s",
    reference_wind_max = [4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8]u"m/s",
    reference_humidity_min = [50.2, 48.4, 48.7, 40.8, 40, 42.1, 45.5, 47.3, 47.6, 45, 51.3, 52.8] ./ 100.0,
    reference_humidity_max = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100] ./ 100.0,
    cloud_min = [50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1] ./ 100.0,
    cloud_max = [50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1] ./ 100.0,
    minima_times = (temp=0, wind=0, humidity=1, cloud=1), # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
    maxima_times = (temp=1, wind=1, humidiy=0, cloud=0), # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
)
    MonthlyMinMaxWeather(reference_temperature_min, reference_temperature_max, reference_wind_min, reference_wind_max, reference_humidity_min, reference_humidity_max, cloud_min, cloud_max, minima_times, maxima_times)
end

function example_soil_thermal_parameters(;
    de_vries_shape_factor = 0.1, # de Vries shape factor, 0.33 for organic soils, 0.1 for mineral
    soil_mineral_conductivity = 1.25u"W/m/K", # soil minerals thermal conductivity
    soil_mineral_density = 2.560u"Mg/m^3", # soil minerals density
    soil_mineral_heat_capacity = 870.0u"J/kg/K", # soil minerals specific heat
    soil_bulk_density = 2.56u"Mg/m^3", # dry soil bulk density
    # TODO why is this calculated elsewhere but also specified here
    soil_saturation_moisture = 0.26u"m^3/m^3", # volumetric water content at saturation (0.1 bar matric potential)
    recirculation_power = 4.0, # power for recirculation function
    return_flow_threshold = 0.162, # return-flow cutoff soil moisture, m^3/m^3
)
    CampbelldeVriesSoilThermal(
        de_vries_shape_factor, soil_mineral_conductivity, soil_mineral_density, soil_mineral_heat_capacity,
        soil_bulk_density, soil_saturation_moisture, recirculation_power, return_flow_threshold
    )
end

# TODO move real defaults to the struct keywords
function example_soil_hydraulics(depths=DEFAULT_DEPTHS;
    bulk_density,
    mineral_density,
    # soil hydraulic parameters
    air_entry_water_potential = fill(0.7, length(depths))u"J/kg", #air entry potential
    saturated_hydraulic_conductivity = fill(0.0058, length(depths))u"kg*s/m^3", #saturated conductivity
    campbell_b_parameter = fill(1.7, length(depths)), #soil 'b' parameter
    soil_bulk_density2 = fill(bulk_density, length(depths))u"Mg/m^3", # soil bulk density
    # TODO what is this why are they different
    soil_mineral_density2 = fill(mineral_density, length(depths))u"Mg/m^3", # soil mineral density
    # plant parameters
    root_density = [0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4, 0.4, 0, 0] * 1e4u"m/m^3", # root density at each node (from Campell 1985 Soil Physics with Basic, p. 131)
    root_resistance = 2.5e+10u"m^3/kg/s", # resistance per unit length of root
    stomatal_closure_potential = -1500.0u"J/kg", # critical leaf water potential for stomatal closure
    leaf_resistance = 2.0e6u"m^4/kg/s", # resistance per unit length of leaf
    stomatal_stability_parameter = 10.0, # stability parameter, -
    root_radius = 0.001u"m", # root radius, m
    # simulation controls
    moist_error = 1e-6u"kg/m^2/s", # maximum overall mass balance error allowed
    moist_count = 500, # maximum iterations of soil moisture algorithm
    moist_step = 360.0u"s", # time step over which to simulate soil moisture (< =  1 hour)
    maxpool = 1.0e4u"kg/m^2", # maximum depth of pooling water
    mode::AbstractSoilMoistureMode = PrescribedSoilMoisture(),
)
    CampbellSoilHydraulics(;
        air_entry_water_potential, saturated_hydraulic_conductivity, campbell_b_parameter, soil_bulk_density2,
        soil_mineral_density2, root_density, root_resistance, stomatal_closure_potential, leaf_resistance, stomatal_stability_parameter,
        root_radius, moist_error, moist_count, moist_step, maxpool, mode,
    )
end
function example_daily_environmental(;
    shade = fill(0.0, length(days)), # fractional shade cast by vegetation
    soil_wetness = fill(0.0, length(days)), # fractional surface wetness
    surface_emissivity = fill(0.96, length(days)), # - surface emissivity
    cloud_emissivity = fill(0.96, length(days)), # - cloud emissivity
    rainfall = ([28, 28.2, 54.6, 79.7, 81.3, 100.1, 101.3, 102.5, 89.7, 62.4, 54.9, 41.2])u"kg/m^2",
    deep_soil_temperature = fill(7.741666u"°C", length(days)),
    leaf_area_index = fill(0.1, length(days)),
)
    EnvironmentTimeseries(; albedo, shade, soil_wetness, surface_emissivity, cloud_emissivity, rainfall, deep_soil_temperature, leaf_area_index)
end

"""
    MicroCache

Pre-allocated workspace for solving a `MicroProblem` without repeated allocation.
Created by `CommonSolve.init(problem::MicroProblem)` and solved in-place with
`CommonSolve.solve!(cache::MicroCache)`.

Use `reinit!(cache, new_problem)` to swap the problem while keeping all
pre-allocated arrays, then call `solve!(cache)` again.

# Example
```julia
cache = init(problem)
output = solve!(cache)

# Change weather and re-solve without allocating
reinit!(cache, modified_problem)
output = solve!(cache)
```
"""
mutable struct MicroCache{MP<:MicroProblem,O<:MicroResult,SR,B,I,SM,N,ND,SP,PB}
    problem::MP
    output::O
    solar_radiation_out::SR
    buffers::B
    ode_integrator::I
    soil_moisture::SM
    nodes::N
    nodes_day::ND
    ∑phase::SP
    profile_buffers::PB
end

function CommonSolve.init(mp::MicroProblem)
    (; days, hours, depths, heights) = mp
    ndays = length(days)
    nhours = length(hours)
    nsteps = ndays * nhours
    num_nodes = length(depths)

    # Solar radiation (allocated once, recomputed on each solve!)
    solar_radiation_out = solve_solar(mp)

    # Output arrays
    output = MicroResult(nsteps, num_nodes, length(heights), solar_radiation_out)

    # Soil buffers
    nodes_day = zeros(num_nodes, ndays)
    nodes_day[1, 1:ndays] .= num_nodes
    nodes = nodes_day[:, 1]
    buffers = (;
        profile = allocate_profile(heights),
        soil_energy_balance = allocate_soil_energy_balance(num_nodes),
        soil_properties = allocate_soil_properties(nodes, mp.soil_thermal_model),
        phase_transition = allocate_phase_transition(length(depths)),
        soil_water_balance = allocate_soil_water_balance(num_nodes),
    )

    # ODE integrator
    soil_moisture = collect(mp.initial_soil_moisture)
    ∑phase = zeros(typeof(1.0u"J"), num_nodes)
    # T0 here is just a prototype for the ODE integrator allocation;
    # solve_soil! will set the real initial condition.
    T0 = if isnothing(mp.initial_soil_temperature)
        SVector(ntuple(_ -> 280.0u"K", num_nodes))
    else
        SVector(ntuple(i -> mp.initial_soil_temperature[i], num_nodes))
    end
    environment_day = get_day(mp.environment_daily, 1)
    environment_instant = get_instant(environment_day, mp.environment_hourly, output, soil_moisture, 1)
    longwave_sky = precompute_longwave_sky(; micro_terrain=mp.micro_terrain, environment_instant, vapour_pressure_equation=mp.vapour_pressure_equation)
    inputs_proto = SoilEnergyInputs(;
        forcing=forcing_day(solar_radiation_out, output, 1),
        buffers, soil_thermal_model=mp.soil_thermal_model, depths, heights,
        solar_terrain=mp.solar_terrain, micro_terrain=mp.micro_terrain,
        nodes, environment_instant,
        soil_wetness=0.0, vapour_pressure_equation=mp.vapour_pressure_equation,
        longwave_sky,
    )
    ode_integrator = allocate_ode_integrator(T0, inputs_proto, mp.soil_ode_solver, mp.soil_ode_kwargs)

    # Profile buffers for solve_air!
    profile_buffers = allocate_profile(heights)

    return MicroCache(mp, output, solar_radiation_out, buffers, ode_integrator, soil_moisture, nodes, nodes_day, ∑phase, profile_buffers)
end

function CommonSolve.solve!(cache::MicroCache)
    mp = cache.problem
    (; environment_minmax, environment_daily, environment_hourly, days, hours, depths, heights) = mp
    output = cache.output
    solar_radiation_out = cache.solar_radiation_out

    # Recompute solar (modifies arrays in-place where possible)
    solve_solar!(solar_radiation_out, mp)

    # Interpolate minmax weather
    interpolate_minmax!(output, environment_minmax, environment_daily, environment_hourly, solar_radiation_out)
    # Adjust solar_radiation given cloud cover to get diffuse fraction
    (; global_radiation, diffuse_fraction) = adjust_for_cloud_cover(output, solar_radiation_out, days, hours; diffuse_fraction_model=mp.diffuse_fraction_model)
    output.diffuse_fraction .= diffuse_fraction
    # Check if cloud-adjusted solar was originally provided
    if !isnothing(environment_hourly)
        if !isnothing(environment_hourly.global_radiation)
            output.global_radiation .= environment_hourly.global_radiation
        else
            output.global_radiation .= global_radiation
        end
    else
        output.global_radiation .= global_radiation
    end
    # Solve soil temperature and moisture
    solve_soil!(cache)
    # Solve air temperatures, windspeed and humidity
    solve_air!(cache)

    return output
end

function CommonSolve.solve(mp::MicroProblem)
    return solve!(init(mp))
end

"""
    reinit!(cache::MicroCache, mp::MicroProblem)

Replace the problem in `cache` while keeping all pre-allocated arrays.
The new problem must have the same dimensions (number of days, hours, depths, heights).
"""
function reinit!(cache::MicroCache, mp::MicroProblem)
    old = cache.problem
    if length(mp.days) != length(old.days) || length(mp.hours) != length(old.hours) ||
       length(mp.depths) != length(old.depths) || length(mp.heights) != length(old.heights)
        throw(DimensionMismatch(
            "reinit! requires identical dimensions. " *
            "Got days=$(length(mp.days)) vs $(length(old.days)), " *
            "hours=$(length(mp.hours)) vs $(length(old.hours)), " *
            "depths=$(length(mp.depths)) vs $(length(old.depths)), " *
            "heights=$(length(mp.heights)) vs $(length(old.heights))"
        ))
    end
    cache.problem = mp
    return cache
end

function solve_solar(mp::MicroProblem)
    (; solar_model, days, hours, solar_terrain) = mp
    # compute clear sky solar radiation
    solar_radiation_out = solar_radiation(solar_model; days, hours, solar_terrain)
    # limit max zenith angles to 90°
    solar_radiation_out.zenith_angle[solar_radiation_out.zenith_angle .> 90u"°"] .= 90u"°"
    solar_radiation_out.zenith_slope_angle[solar_radiation_out.zenith_slope_angle .> 90u"°"] .= 90u"°"

    return solar_radiation_out
end

# In-place recompute: copy new solar values into existing pre-allocated arrays
function solve_solar!(existing::NamedTuple, mp::MicroProblem)
    fresh = solve_solar(mp)
    for k in keys(existing)
        dest = getfield(existing, k)
        src = getfield(fresh, k)
        if dest isa AbstractArray
            dest .= src
        end
    end
    return existing
end

function interpolate_minmax!(output, environment_minmax, environment_daily, environment_hourly, solar_radiation_out)
    # interpolate daily min/max forcing variables to hourly
    reference_temperature, reference_wind_speed, reference_humidity, cloud_cover = hourly_from_min_max(environment_minmax, solar_radiation_out)
    # TODO just use loops for these this allocates
    reference_humidity[reference_humidity .> 1.0] .= 1.0
    cloud_cover[cloud_cover .> 1.0] .= 1.0

    output.cloud_cover .= cloud_cover
    output.reference_temperature .= reference_temperature
    output.pressure .= environment_hourly.pressure
    output.reference_wind_speed .= reference_wind_speed
    output.reference_humidity .= reference_humidity

    return 
end
function interpolate_minmax!(output, environment_minmax::Nothing, environment_daily, environment_hourly, solar_radiation_out)
    output.cloud_cover .= environment_hourly.cloud_cover
    output.reference_temperature .= environment_hourly.reference_temperature
    output.pressure .= environment_hourly.pressure
    output.reference_wind_speed .= environment_hourly.reference_wind_speed
    output.reference_humidity .= environment_hourly.reference_humidity

    return 
end

function adjust_for_cloud_cover(output, solar_radiation_out, days, hours; diffuse_fraction_model=ErbsDiffuseFraction())
    # adjust for cloud using Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
    day_of_year = repeat(days, inner=length(hours))
    zenith_angle = solar_radiation_out.zenith_angle
    direct_horizontal = solar_radiation_out.direct_horizontal
    diffuse_horizontal = solar_radiation_out.diffuse_horizontal
    cloud = output.cloud_cover
    return (; global_radiation, diffuse_fraction) = cloud_adjust_radiation(output, cloud, diffuse_horizontal, direct_horizontal, zenith_angle, day_of_year; diffuse_fraction_model)
end

# Solves soil temperature and moisture using pre-allocated cache buffers
function solve_soil!(cache::MicroCache)
    mp = cache.problem
    (; days, hours, depths, heights) = mp
    output = cache.output
    solar_radiation_out = cache.solar_radiation_out
    buffers = cache.buffers
    ode_integrator = cache.ode_integrator
    soil_moisture = cache.soil_moisture
    nodes = cache.nodes
    nodes_day = cache.nodes_day
    ∑phase = cache.∑phase

    (; solar_terrain, micro_terrain, soil_thermal_model, soil_moisture_model, environment_minmax, environment_daily, initial_soil_temperature, initial_soil_moisture, hourly_rainfall, vapour_pressure_equation, soil_ode_solver, soil_ode_kwargs, time_mode, convergence) = mp
    moisture_mode = soil_moisture_model.mode
    (; moist_step, campbell_b_parameter, soil_bulk_density2, soil_mineral_density2, air_entry_water_potential) = soil_moisture_model
    init_soil_wetness!(moisture_mode)

    ndays = length(days)
    nhours = length(hours)
    num_nodes = length(depths) # number of soil nodes for temperature and soil moisture calcs

    # initial conditions
    if independent_days(time_mode) && isnothing(initial_soil_temperature)
        # TODO below we add the first value to the other 24 so that it mimics the Fortran 25h day but
        # ultimately we should remove that extra value
        t = mean(u"K", [view(output.reference_temperature, 1:nhours); output.reference_temperature[1]])
        T0 = SVector(ntuple(_ -> t, num_nodes))
    else
        if num_nodes != length(initial_soil_temperature)
            error("Initial soil temperature must match length of 'depths'")
        end
        T0 = SVector(ntuple(i -> initial_soil_temperature[i], num_nodes))
    end

    # Reset pre-allocated scratch arrays
    nodes_day[1, 1:ndays] .= num_nodes
    ∑phase .= 0.0u"J"
    soil_moisture .= initial_soil_moisture

    nodes .= nodes_day[:, 1]
    update_soil_properties!(output, buffers.soil_properties, soil_thermal_model;
        soil_temperature=T0, soil_moisture, atmospheric_pressure=101325.0u"Pa", step=1
    )

    soil_saturation_moisture = 1.0 .- soil_bulk_density2 ./ soil_mineral_density2
    output.soil_water_potential[1, :] .= air_entry_water_potential .* (soil_saturation_moisture ./ soil_moisture) .^ campbell_b_parameter
    output.soil_temperature[1, :] .= T0
    output.soil_moisture[1, :] = soil_moisture

    initialise_soil_humidity!(moisture_mode, output, output.soil_water_potential[1, :], T0)

    environment_day = get_day(environment_daily, 1)
    environment_instant = get_instant(environment_day, mp.environment_hourly, output, soil_moisture, 1)
    longwave_sky = precompute_longwave_sky(; micro_terrain, environment_instant, vapour_pressure_equation)

    # simulate all days
    pool = 0.0u"kg/m^2" # initialise depth of pooling water TODO make this an init option
    niter_moist = ustrip(u"s^-1", 3600 / moist_step) # TODO use a solver for soil moisture calc
    infil_out = nothing
    for j in 1:ndays
        iday = j
        nodes .= nodes_day[:, iday]
        environment_day = get_day(environment_daily, iday)
        forcing = forcing_day(solar_radiation_out, output, iday)
        niter = iterations_for_day(time_mode, convergence, j)
        if independent_days(time_mode)
            ∑phase .= 0.0u"J"
            sub2 = (iday*nhours-nhours+1):(iday*nhours) # for getting mean monthly over the 25 hrs as in fortran version
            if isnothing(initial_soil_temperature)
                t = mean(u"K", [output.reference_temperature[sub2]; output.reference_temperature[sub2][1]])
                T0 = SVector(ntuple(_ -> t, num_nodes))
            else
                T0 = SVector(ntuple(i -> initial_soil_temperature[i], num_nodes))
            end
            reset_day_soil_moisture!(moisture_mode, soil_moisture, initial_soil_moisture, j)
        end
        T0 = setindex(T0, environment_instant.deep_soil_temperature, num_nodes) # set deepest node to boundary condition
        use_multi_iter  = may_iterate(convergence)
        iter            = 0
        is_last_iter    = false
        T0_prev_start   = T0  # T0 at start of the previous outer iteration (for convergence check)
        while !is_last_iter
            iter += 1
            T0_iter_start = T0  # save for next iteration's convergence check
            is_last_iter = is_converged(convergence, iter, niter, T0, T0_prev_start)
            T0_prev_start = T0_iter_start
            ∑phase .= 0.0u"J" # TODO decide whether this should happen and fix in Fortran
            for i in 1:length(hours) # loop through hours of day
                step = (j - 1) * length(hours) + i
                environment_instant = get_instant(environment_day, mp.environment_hourly, output, soil_moisture, step)
                longwave_sky = precompute_longwave_sky(; micro_terrain, environment_instant, vapour_pressure_equation)
                inputs = SoilEnergyInputs(; forcing, buffers, soil_thermal_model, depths, heights, solar_terrain, micro_terrain, nodes, environment_instant, soil_wetness=get_soil_wetness(moisture_mode, environment_instant), vapour_pressure_equation, longwave_sky)
                # Precompute soil properties once per hour using T0; ODE reads buffers directly
                soil_properties!(buffers.soil_properties, soil_thermal_model;
                    soil_temperature=T0, soil_moisture,
                    atmospheric_pressure=environment_instant.atmospheric_pressure,
                    vapour_pressure_equation,
                )
                if i == 1 # make first hour of day equal last hour of previous iteration
                    if use_multi_iter
                        ∑phase .= 0.0u"J" # TODO decide whether this should happen and fix in Fortran
                        T0_before = T0
                        T0 = get_soil_temp_timeline!(ode_integrator, T0, i + 1, inputs)
                        if is_last_iter # TODO this should happen every time but at present it doesn't in Fortran version
                            (; accumulated_latent_heat, phase_change_heat, temperature) = phase_transition!(buffers.phase_transition;
                                temperatures=T0, temperatures_past=T0_before, accumulated_latent_heat=∑phase, soil_moisture, depths
                            )
                            ∑phase = accumulated_latent_heat
                            T0 = temperature
                        end
                    end
                else
                    T0_before = T0
                    T0 = get_soil_temp_timeline!(ode_integrator, T0, i, inputs)
                    if is_last_iter # TODO this should happen every time but at present it doesn't in Fortran version
                        (; accumulated_latent_heat, phase_change_heat, temperature) = phase_transition!(buffers.phase_transition;
                            temperatures=T0, temperatures_past=T0_before, accumulated_latent_heat=∑phase, soil_moisture, depths
                        )
                        ∑phase = accumulated_latent_heat
                        T0 = temperature
                    end
                end
                init_soil_obukhov!(buffers, forcing, micro_terrain, heights, T0, i)
                rain = hourly_rainfall ? mp.environment_hourly.rainfall[step] : environment_instant.rainfall
                pool = clamp(pool + rain, 0.0u"kg/m^2", soil_moisture_model.maxpool)
                (; pool, soil_moisture) = step_soil_moisture!(moisture_mode, buffers, soil_moisture_model, output, step;
                    depths, micro_terrain, environment_instant, T0, niter_moist, pool, soil_moisture, vapour_pressure_equation
                )
                # Write to output
                if is_last_iter
                    output.surface_water[step] = pool
                    output.soil_temperature[step, :] .= T0
                    output.sky_temperature[step] = longwave_sky.sky_temperature
                    environment_instant = get_instant(environment_day, mp.environment_hourly, output, soil_moisture, step)
                    update_soil_properties!(output, buffers.soil_properties, soil_thermal_model;
                        soil_temperature=T0, soil_moisture, atmospheric_pressure=environment_instant.atmospheric_pressure, step, vapour_pressure_equation
                    )
                end
            end
        end
    end
    return output
end

# compute air temperature, wind speed and relative humidity profiles
function solve_air!(cache::MicroCache)
    mp = cache.problem
    output = cache.output
    solar_radiation_out = cache.solar_radiation_out
    (; micro_terrain, vapour_pressure_equation) = mp
    profile_buffers = cache.profile_buffers
    for i in 1:size(output.profile.air_temperature, 1)
        surface_temperature = u"°C"(output.soil_temperature[i, 1])
        environment_instant = (;
            atmospheric_pressure=output.pressure[i],
            reference_temperature=output.reference_temperature[i],
            reference_wind_speed=output.reference_wind_speed[i],
            reference_humidity=output.reference_humidity[i],
            zenith_angle=solar_radiation_out.zenith_angle[i],
        )
        result = atmospheric_surface_profile!(profile_buffers; micro_terrain, environment_instant, surface_temperature, vapour_pressure_equation)
        output.profile.air_temperature[i, :]   .= result.air_temperature
        output.profile.wind_speed[i, :]        .= result.wind_speed
        output.profile.relative_humidity[i, :] .= result.relative_humidity
        output.profile.convective_heat_flux[i]  = result.convective_heat_flux
        output.profile.friction_velocity[i]     = result.friction_velocity
    end
end

function allocate_ode_integrator(T0, inputs_proto, solver, solver_kwargs)
    tspan = (0.0u"minute", 60.0u"minute")
    prob = ODEProblem{false}(soil_energy_balance, T0, tspan, inputs_proto)
    return SciMLBase.init(prob, solver; save_everystep=false, save_start=false, save_end=false, solver_kwargs...)
end

function get_soil_temp_timeline!(integrator, T0, i, inputs)
    t0 = (0.0 + (i - 2) * 60)u"minute"
    tf = (60.0 + (i - 2) * 60)u"minute"
    integrator.p = inputs
    SciMLBase.reinit!(integrator, T0; t0, tf, reset_dt=false)
    SciMLBase.solve!(integrator)
    return integrator.u
end

function update_soil_properties!(output, soil_properties_buffers, soil_thermal_model; step, kw...)
    (; bulk_thermal_conductivity, bulk_heat_capacity, bulk_density) = soil_properties!(soil_properties_buffers, soil_thermal_model; kw...)

    output.soil_thermal_conductivity[step, :] .= bulk_thermal_conductivity
    output.soil_heat_capacity[step, :] .= bulk_heat_capacity
    output.soil_bulk_density[step, :] .= bulk_density

    return output
end

function update_soil_water!(output, infil_out, step)
    output.soil_moisture[step, :] .= infil_out.soil_moisture
    output.soil_water_potential[step, :] .= infil_out.soil_water_potential
    output.soil_humidity[step, :] .= infil_out.soil_humidity

    return output
end

# ── Soil moisture stepping dispatch ───────────────────────────────────────

function step_soil_moisture!(mode::DynamicSoilMoisture, buffers, soil_moisture_model, output, step;
    depths, micro_terrain, environment_instant, T0, niter_moist, pool, soil_moisture, vapour_pressure_equation
)
    (; infil_out, soil_wetness, pool, soil_moisture) = get_soil_water_balance!(buffers, soil_moisture_model;
        depths, micro_terrain, environment_instant, T0, niter_moist, pool,
        soil_wetness=mode.soil_wetness, soil_moisture, vapour_pressure_equation
    )
    mode.soil_wetness = soil_wetness
    update_soil_water!(output, infil_out, step)
    return (; pool, soil_moisture)
end

step_soil_moisture!(::PrescribedSoilMoisture, args...; pool, soil_moisture, kw...) = (; pool, soil_moisture)

# TODO eventually make environment a type,
# and this can just be `getindex` on that type.
function forcing_day(solar_radiation_out, output, iday::Int)
    (; pressure, reference_temperature, reference_wind_speed, reference_humidity, cloud_cover) = output
    global_radiation = output.global_radiation
    (; zenith_angle, zenith_slope_angle) = solar_radiation_out

    nhours = 24
    sub1 = (iday*nhours-nhours+1):(iday*nhours)
    tspan = 0.0:60:(60*(nhours-1))

    # get today's weather
    interpolate_solar = scale(interpolate(view(global_radiation, sub1), BSpline(Linear())), tspan)
    interpolate_zenith = scale(interpolate(view(zenith_angle, sub1), BSpline(Linear())), tspan)
    interpolate_slope_zenith = scale(interpolate(view(zenith_slope_angle, sub1), BSpline(Linear())), tspan)
    interpolate_temperature = scale(interpolate(u"K".(view(reference_temperature, sub1)), BSpline(Linear())), tspan)
    interpolate_wind = scale(interpolate(view(reference_wind_speed, sub1), BSpline(Linear())), tspan)
    interpolate_humidity = scale(interpolate(view(reference_humidity, sub1), BSpline(Linear())), tspan)
    interpolate_cloud = scale(interpolate(view(cloud_cover, sub1), BSpline(Linear())), tspan)
    interpolate_pressure = scale(interpolate(view(pressure, sub1), BSpline(Linear())), tspan)

    return MicroForcing(; interpolate_solar, interpolate_zenith, interpolate_slope_zenith, interpolate_temperature, interpolate_wind, interpolate_humidity, interpolate_cloud, interpolate_pressure)
end


# TODO these functions are a bit silly
function get_day(environment_daily, iday)
    # TODO: standardise all these names
    environment_day = (;
        leaf_area_index = environment_daily.leaf_area_index[iday],
        shade = environment_daily.shade[iday], # daily shade (fractional)
        surface_emissivity = environment_daily.surface_emissivity[iday],
        cloud_emissivity = environment_daily.cloud_emissivity[iday], # - cloud emissivity
        soil_wetness = environment_daily.soil_wetness[iday], # set up vector of soil wetness for each day
        deep_soil_temperature = u"K"(environment_daily.deep_soil_temperature[iday]), # daily deep soil temperature (°C)
        rainfall = environment_daily.rainfall[iday],
    )
end
function get_instant(environment_day, environment_hourly, output, soil_moisture, i)
    return (;
        environment_day...,
        # TODO getting data from output means it being correct depends on
        # order of operations in the solve. We need an itermediate object instead
        atmospheric_pressure = output.pressure[i],
        reference_temperature = output.reference_temperature[i],
        reference_wind_speed = output.reference_wind_speed[i],
        reference_humidity = output.reference_humidity[i],
        zenith_angle = output.solar_radiation.zenith_angle[i],
        cloud_cover = output.cloud_cover[i],
        global_radiation = output.global_radiation[i],
        soil_moisture=soil_moisture,
    )
end

# This handles getting values from a Number or array of numbers, or objects of these
maybegetindex(obj::CampbelldeVriesSoilThermal, i::Int) = CampbelldeVriesSoilThermal(; maybegetindex(ConstructionBase.getproperties(obj), i)...)
maybegetindex(props::NamedTuple, i::Int) = map(p -> maybegetindex(p, i), props)
maybegetindex(val::Number, i::Int) = val
maybegetindex(vals::AbstractArray, i::Int) = vals[i]
