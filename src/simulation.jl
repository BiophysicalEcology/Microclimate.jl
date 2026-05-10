"""
    MicroProblem(; site, parameters, environment_daily, ..., config=MicroConfig())

Microclimate simulation problem specification.

Top-level structure:
- `site::Site` — properties of the place
- `parameters::MicroParameters` — physical-parameter structs (soil thermal, soil hydraulics, snow)
- `environment_minmax`, `environment_daily`, `environment_hourly` — input forcings
- `initial_*` — initial conditions
- `config::MicroConfig` — every formulation choice and solver tuning knob
"""
@kwdef struct MicroProblem{D,H,Dep,Ht,SM,S,P,EMM,EH,ED,IST,ISM,ISD,ISNT,ISND,C}
    # locations, times, depths and heights
    days::D = DEFAULT_DAYS # days of year to simulate - TODO leap years - why not use real dates?
    hours::H = DEFAULT_HOURS # hour of day for solar_radiation
    depths::Dep = DEFAULT_DEPTHS # soil nodes - keep spacing close near the surface
    heights::Ht = [0.01, 2]u"m" # air nodes for temperature, wind speed and humidity profile, last height is reference height for weather data
    # Place
    solar_model::SM = SolarProblem()
    site::S
    # All physical-parameter structs (soil thermal, soil hydraulics, snow)
    parameters::P
    # Inputs
    environment_minmax::EMM
    environment_hourly::EH = nothing
    environment_daily::ED
    # Initial conditions
    initial_soil_temperature::IST = fill(u"K"(7.741667u"°C"), length(depths))
    initial_soil_moisture::ISM = fill(0.42 * 0.25, length(depths))
    initial_snow_depth::ISD = 0.0u"cm"
    initial_snow_temperature::ISNT = u"K"(0.0u"°C")
    initial_snow_density::ISND = nothing  # nothing → use parameters.snow.snow_density
    # All formulation choices and solver tuning
    config::C = MicroConfig()
end

function example_microclimate_problem(;
    days = DEFAULT_DAYS,
    hours = DEFAULT_HOURS,
    site = example_site(),
    parameters = MicroParameters(;
        soil_thermal = example_soil_thermal_parameters(),
        soil_hydraulics = example_soil_hydraulics(),
    ),
    environment_minmax = example_monthly_weather(),
    environment_daily = example_daily_environment(days),
    environment_hourly = example_hourly_environment(days, hours; elevation=site.elevation),
    config = MicroConfig(),
    kw...
)
    MicroProblem(;
        days, hours, site, parameters,
        environment_minmax, environment_daily, environment_hourly,
        config,
        kw...,
    )
end

function example_site(;
    latitude = 43.07305u"°",
    longitude = -89.40123u"°",
    elevation = 226.0u"m",
    slope = 0.0u"°",
    aspect = 0.0u"°",
    horizon_angles = fill(0.0u"°", 24),
    sky_view_fraction = 1.0,
    albedo = 0.15,
    roughness_height = 0.004u"m",
    atmos_pressure = atmospheric_pressure(elevation),
)
    Site(;
        latitude, longitude, elevation, slope, aspect, horizon_angles,
        sky_view_fraction, albedo, roughness_height,
        atmospheric_pressure=atmos_pressure,
    )
end

function example_monthly_weather(;
    reference_temperature_min = [-14.3, -12.1, -5.1, 1.2, 6.9, 12.3, 15.2, 13.6, 8.9, 3, -3.2, -10.6]u"°C",
    reference_temperature_max = [-3.2, 0.1, 6.8, 14.6, 21.3, 26.4, 29, 27.7, 23.3, 16.6, 7.8, -0.4]u"°C",
    reference_wind_min = [4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8] * 0.1u"m/s",
    reference_wind_max = [4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8]u"m/s",
    reference_humidity_min = [50.2, 48.4, 48.7, 40.8, 40, 42.1, 45.5, 47.3, 47.6, 45, 51.3, 52.8] ./ 100.0,
    reference_humidity_max = fill(1.0, 12),
    cloud_min = [50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1] ./ 100.0,
    cloud_max = [50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1] ./ 100.0,
    # Per-variable hour offsets indexed [air_temp, wind, humidity, cloud].
    # Air temp & wind: minima offset from sunrise, maxima offset from solar noon.
    # Humidity & cloud: minima offset from solar noon, maxima offset from sunrise.
    minima_times = [0, 0, 1, 1],
    maxima_times = [1, 1, 0, 0],
)
    MonthlyMinMaxEnvironment(;
        reference_temperature_min, reference_temperature_max,
        reference_wind_min, reference_wind_max,
        reference_humidity_min, reference_humidity_max,
        cloud_min, cloud_max,
        minima_times, maxima_times,
    )
end

function example_soil_thermal_parameters(;
    de_vries_shape_factor = 0.1, # de Vries shape factor, 0.33 for organic soils, 0.1 for mineral
    mineral_conductivity = 1.25u"W/m/K", # soil minerals thermal conductivity
    mineral_heat_capacity = 870.0u"J/kg/K", # soil minerals specific heat
    recirculation_power = 4.0, # power for recirculation function
    return_flow_threshold = 0.162, # return-flow cutoff soil moisture, m^3/m^3
)
    CampbelldeVriesSoilThermal(;
        de_vries_shape_factor, mineral_conductivity, mineral_heat_capacity,
        recirculation_power, return_flow_threshold,
    )
end

# TODO move real defaults to the struct keywords
function example_soil_hydraulics(depths=DEFAULT_DEPTHS;
    # Scalars get broadcast to a per-depth vector; an AbstractVector is taken as-is.
    bulk_density = 2.56u"Mg/m^3",
    mineral_density = 2.560u"Mg/m^3",
    # soil hydraulic parameters
    air_entry_water_potential = fill(0.7, length(depths))u"J/kg", #air entry potential
    saturated_hydraulic_conductivity = fill(0.0058, length(depths))u"kg*s/m^3", #saturated conductivity
    campbell_b_parameter = fill(1.7, length(depths)), #soil 'b' parameter
    # plant parameters
    root_density = [0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4, 0.4, 0, 0] * 1e4u"m/m^3", # root density at each node (from Campell 1985 Soil Physics with Basic, p. 131)
    root_resistance = 2.5e+10u"m^3/kg/s", # resistance per unit length of root
    stomatal_closure_potential = -1500.0u"J/kg", # critical leaf water potential for stomatal closure
    leaf_resistance = 2.0e6u"m^4/kg/s", # resistance per unit length of leaf
    stomatal_stability_parameter = 10.0, # stability parameter, -
    root_radius = 0.001u"m", # root radius, m
)
    CampbellSoilHydraulics(;
        air_entry_water_potential, saturated_hydraulic_conductivity, campbell_b_parameter,
        bulk_density = bulk_density isa AbstractVector ? bulk_density : fill(bulk_density, length(depths)),
        mineral_density = mineral_density isa AbstractVector ? mineral_density : fill(mineral_density, length(depths)),
        root_density, root_resistance, stomatal_closure_potential, leaf_resistance, stomatal_stability_parameter,
        root_radius,
    )
end
function example_daily_environment(days=DEFAULT_DAYS;
    shade = fill(0.0, length(days)), # fractional shade cast by vegetation
    soil_wetness = fill(0.0, length(days)), # fractional surface wetness
    surface_emissivity = fill(0.96, length(days)), # - surface emissivity
    cloud_emissivity = fill(0.96, length(days)), # - cloud emissivity
    rainfall = ([28, 28.2, 54.6, 79.7, 81.3, 100.1, 101.3, 102.5, 89.7, 62.4, 54.9, 41.2])u"kg/m^2",
    deep_soil_temperature = fill(7.741666u"°C", length(days)),
    leaf_area_index = fill(0.1, length(days)),
)
    DailyTimeseries(;
        shade, soil_wetness, surface_emissivity, cloud_emissivity,
        rainfall, deep_soil_temperature, leaf_area_index,
    )
end

function example_hourly_environment(days=DEFAULT_DAYS, hours=DEFAULT_HOURS;
    elevation = 226.0u"m",
    pressure = fill(atmospheric_pressure(elevation), length(days) * length(hours)),
    reference_temperature = nothing,
    reference_humidity = nothing,
    reference_wind_speed = nothing,
    global_radiation = nothing,
    cloud_cover = nothing,
    rainfall = nothing,
    zenith_angle = nothing,
    longwave_radiation = nothing,
)
    HourlyTimeseries(;
        pressure, reference_temperature, reference_humidity, reference_wind_speed,
        global_radiation, cloud_cover, rainfall, zenith_angle, longwave_radiation,
    )
end

"""
    MicroState

Mutable simulation state that evolves during `solve!`. Lives on `MicroCache.state`.
"""
mutable struct MicroState{NS,SM,N,ND,SP}
    snow_state::NS
    soil_moisture::SM
    nodes::N
    nodes_day::ND
    ∑phase::SP
end

"""
    MicroBuffers

Pre-allocated workspace. Every buffer the hot path touches lands here once
in `init`. Lives on `MicroCache.buffers`.
"""
struct MicroBuffers{SO,SOB,P,PB,SEB,SP,PT,SWB,SS}
    solar_out::SO                  # SolarRadiation output (NamedTuple of arrays)
    solar::SOB                     # SolarRadiation internal buffers (NamedTuple)
    profile::P                     # atmospheric profile scratch used by the moisture solver
    air_profile::PB                # atmospheric profile scratch used by solve_air!
    soil_energy_balance::SEB
    soil_properties::SP
    phase_transition::PT
    soil_water_balance::SWB
    snow_scratch::SS               # NamedTuple for SnowModel; nothing for NoSnow
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
mutable struct MicroCache{MP<:MicroProblem,O<:MicroResult,S<:MicroState,B<:MicroBuffers,F<:MicroForcing,I,Nsoil}
    problem::MP
    output::O
    state::S
    buffers::B
    forcing::F
    ode_integrator::I
    # Type-level soil-node count, used for `SVector(ntuple(f, Val(N_soil)))` in
    # the day loop so T0 reconstruction is type-stable.
    num_soil_nodes_val::Val{Nsoil}
end

function populate_weather!(output, mp, solar_radiation_out)
    (; days, hours) = mp
    interpolate_minmax!(output, mp.environment_minmax, mp.environment_daily, mp.environment_hourly, solar_radiation_out)
    (; global_radiation, diffuse_fraction) = adjust_for_cloud_cover(output, solar_radiation_out, days, hours; diffuse_fraction_model=mp.config.diffuse_fraction_model)
    output.diffuse_fraction .= diffuse_fraction
    if !isnothing(mp.environment_hourly) && !isnothing(mp.environment_hourly.global_radiation)
        output.global_radiation .= mp.environment_hourly.global_radiation
    else
        output.global_radiation .= global_radiation
    end
    return output
end

# NoSnow has no surface overrides — return the base values directly. Dispatching
# on snow_model type means the SnowModel-only override / NamedTuple-merge code
# is never compiled into the NoSnow execution path.
@inline function apply_snow_overrides(::NoSnow, _, _, _, site, moisture_mode, environment_instant, vapour_pressure_equation)
    return (;
        albedo = site.albedo,
        effective_wetness = get_soil_wetness(moisture_mode, environment_instant),
        longwave_sky = precompute_longwave_sky(; site, environment_instant, vapour_pressure_equation),
    )
end

function apply_snow_overrides(snow_model::SnowModel, snow_state, snow_scratch, step, site, moisture_mode, environment_instant, vapour_pressure_equation)
    overrides = snow_surface_overrides(snow_model, snow_state, snow_scratch, step)
    albedo = isnothing(overrides.albedo) ? site.albedo : overrides.albedo
    effective_wetness = isnothing(overrides.soil_wetness) ? get_soil_wetness(moisture_mode, environment_instant) : overrides.soil_wetness
    env_for_longwave = environment_instant
    if !isnothing(overrides.shade)
        env_for_longwave = merge(env_for_longwave, (; shade=overrides.shade))
    end
    if !isnothing(overrides.emissivity)
        env_for_longwave = merge(env_for_longwave, (; surface_emissivity=overrides.emissivity))
    end
    longwave_sky = precompute_longwave_sky(; site, environment_instant=env_for_longwave, vapour_pressure_equation)
    return (; albedo, effective_wetness, longwave_sky)
end

function build_soil_energy_inputs(; forcing, buffers, soil_thermal_model, depths, heights, site, boundary_layer_model,
    n_snow, num_soil_nodes, nodes, environment_instant, effective_wetness, vapour_pressure_equation,
    longwave_sky, albedo, qfreze, snow_model, snow_state, snow_scratch, soil_moisture,
    bulk_density, mineral_density, maximum_surface_temperature,
)
    return SoilEnergyInputs(; forcing, buffers, soil_thermal_model,
        depths, heights, site, boundary_layer_model,
        nodes=n_snow > 0 ? zeros(n_snow + num_soil_nodes) : nodes,
        environment_instant, soil_wetness=effective_wetness,
        vapour_pressure_equation, longwave_sky, albedo, qfreze,
        snow_model, snow_state, snow_scratch, soil_moisture,
        bulk_density, mineral_density, n_snow, maximum_surface_temperature,
    )
end

function sync_inactive_snow_temps(T_snow::SVector{N}, snow_scratch) where N
    sync = T_snow[1]
    return SVector(ntuple(N) do k
        snow_scratch.node_depths[k] > 0.0u"cm" ? T_snow[k] : sync
    end)
end

# Hoisted to avoid closure-capture of `t` when it's also reassigned inside
# `solve_soil!` (which would force Julia to box `t` for the whole function).
@inline _filled_svec(val, ::Val{N}) where N = SVector(ntuple(_ -> val, Val(N)))

# Hoisted to avoid closure-capture of `T_snow`. Dispatches on `first_active`
# type so the Union{Nothing,Int} from `findfirst` resolves to one method per
# branch — no Union propagation through the calling function.
@inline _sync_snow_to_active(T_snow::SVector{N,T}, ::Nothing, sync_temp) where {N,T} =
    _filled_svec(sync_temp, Val(N))
@inline _sync_snow_to_active(T_snow::SVector{N,T}, first_active::Int, sync_temp) where {N,T} =
    SVector(ntuple(k -> k < first_active ? sync_temp : T_snow[k], Val(N)))

# Dispatched on T_snow's type so the call site infers a concrete SVector
# return rather than `Union{SVector{N1,T}, SVector{N2,T}}`.
combine_ode_state(::Nothing, T0, _) = T0
combine_ode_state(T_snow::SVector, T0, _) = vcat(T_snow, T0)

# Mean of `arr[start:stop]` plus one repeat of the final element (Fortran
# SINEC's TINS uses 25 nodes — 24 hourly values + the last repeated). Avoids
# `[arr[r]; arr[r][1]]` concatenation which heap-allocates intermediate Vectors.
@inline function _tins_mean_K(arr, start::Int, stop::Int)
    n = stop - start + 1
    @inbounds s = u"K"(arr[start])
    @inbounds for i in (start+1):stop
        s += u"K"(arr[i])
    end
    @inbounds s += u"K"(arr[stop])  # 25th sample = repeat of last hour
    return s / (n + 1)
end

# Dispatched split of the integrator state back into (T_snow, T_soil). For
# NoSnow returns (nothing, T_result). For SnowModel{N_snow}, peels off the
# first N_snow entries as T_snow and the rest as T_soil. All sizes are
# compile-time constants from snow_model and T_result's type parameters.
@inline split_ode_state(::NoSnow, T_result) = (nothing, T_result)
@inline function split_ode_state(::SnowModel{N_snow}, T_result::SVector{N_total,T}) where {N_snow, N_total, T}
    T_snow = SVector(ntuple(k -> T_result[k], Val(N_snow)))
    T_soil = SVector(ntuple(k -> T_result[N_snow + k], Val(N_total - N_snow)))
    return T_snow, T_soil
end

function CommonSolve.init(mp::MicroProblem)
    (; days, hours, depths, heights) = mp
    snow_model = mp.parameters.snow
    n_snow = n_snow_nodes(snow_model)
    ndays = length(days)
    nhours = length(hours)
    nsteps = ndays * nhours
    num_soil_nodes = length(depths)
    num_ode_nodes = n_snow + num_soil_nodes

    # Solar radiation: pre-allocate output and scratch buffers once, reuse on
    # each solve! via solar_radiation!.
    solar_radiation_out, solar_buffers = allocate_solar(mp)
    solve_solar!(solar_radiation_out, solar_buffers, mp)

    # Output arrays (always soil-node sized)
    output = MicroResult(nsteps, num_soil_nodes, length(heights), solar_radiation_out, n_snow)

    # Snow state and scratch
    snow_state = initial_snow_state(snow_model, mp.initial_snow_depth, mp.initial_snow_density)
    snow_scratch = allocate_snow_scratch(snow_model, nsteps, num_ode_nodes, depths)

    # Fill output weather arrays with placeholder values so the prototype ODE integrator
    # can be allocated without NaNs. Real values are set at the start of solve!.
    fill!(output.pressure, 101325.0u"Pa")
    fill!(output.reference_temperature, 293.0u"K")
    fill!(output.reference_wind_speed, 1.0u"m/s")
    fill!(output.reference_humidity, 0.5)
    fill!(output.cloud_cover, 0.0)
    fill!(output.global_radiation, 300.0u"W/m^2")

    # All up-front allocation lands in one MicroBuffers struct.
    nodes_day = zeros(num_soil_nodes, ndays)
    nodes_day[1, 1:ndays] .= num_soil_nodes
    nodes = nodes_day[:, 1]
    buffers = MicroBuffers(
        solar_radiation_out,
        solar_buffers,
        allocate_profile(heights),
        allocate_profile(heights),
        allocate_soil_energy_balance(num_ode_nodes),
        allocate_soil_properties(zeros(num_ode_nodes), mp.parameters.soil_thermal, mp.parameters.soil_hydraulics),
        allocate_phase_transition(num_soil_nodes),
        allocate_soil_water_balance(num_soil_nodes),
        snow_scratch,
    )

    # ODE integrator — sized for the combined (snow + soil) vector
    soil_moisture = collect(mp.initial_soil_moisture)
    ∑phase = zeros(typeof(1.0u"J"), num_soil_nodes)
    T0_soil = if isnothing(mp.initial_soil_temperature)
        SVector(ntuple(_ -> 280.0u"K", num_soil_nodes))
    else
        SVector(ntuple(i -> mp.initial_soil_temperature[i], num_soil_nodes))
    end
    T0_snow = init_snow_temperatures(snow_model, mp.initial_snow_temperature)
    T0_ode = combine_ode_state(T0_snow, T0_soil, n_snow)

    # Build prototype depths for ODE (combined snow + soil)
    if n_snow > 0
        depths_proto = collect(combined_depths(snow_model, snow_state, snow_scratch, depths))
    else
        depths_proto = collect(depths)
    end

    populate_weather!(output, mp, solar_radiation_out)

    # Pre-allocate the day-forcing interpolations once. The 8 ScaledInterpolations
    # wrap 25-element coefficient buffers that get refilled per day-iter via
    # `update_forcing_day!` — no per-day construction.
    forcing = allocate_forcing(output, solar_radiation_out)
    update_forcing_day!(forcing, solar_radiation_out, output, 1)

    environment_day = get_day(mp.environment_daily, 1)
    environment_instant = get_instant(environment_day, mp.environment_hourly, output, soil_moisture, 1)
    longwave_sky = precompute_longwave_sky(; site=mp.site, environment_instant, vapour_pressure_equation=mp.config.vapour_pressure_equation)
    inputs_proto = SoilEnergyInputs(;
        forcing,
        buffers, soil_thermal_model=mp.parameters.soil_thermal, depths=depths_proto, heights,
        site=mp.site, boundary_layer_model=mp.config.boundary_layer_model,
        nodes=zeros(num_ode_nodes), environment_instant,
        soil_wetness=0.0, vapour_pressure_equation=mp.config.vapour_pressure_equation,
        longwave_sky, albedo=mp.site.albedo,
        soil_moisture,
        bulk_density=mp.parameters.soil_hydraulics.bulk_density,
        mineral_density=mp.parameters.soil_hydraulics.mineral_density,
        n_snow=n_snow,
        snow_model, snow_state=initial_snow_state(snow_model), snow_scratch,
        maximum_surface_temperature=mp.config.maximum_surface_temperature,
    )
    ode_integrator = allocate_ode_integrator(T0_ode, inputs_proto, mp.config.soil_ode_solver, mp.config.soil_ode_kwargs)

    state = MicroState(snow_state, soil_moisture, nodes, nodes_day, ∑phase)

    return MicroCache(mp, output, state, buffers, forcing, ode_integrator, Val(num_soil_nodes))
end

function CommonSolve.solve!(cache::MicroCache)
    mp = cache.problem
    (; environment_minmax, environment_daily, environment_hourly, days, hours, depths, heights) = mp
    output = cache.output
    solar_radiation_out = cache.buffers.solar_out

    # Recompute solar in place into the pre-allocated cache buffers.
    solve_solar!(solar_radiation_out, cache.buffers.solar, mp)

    populate_weather!(output, mp, solar_radiation_out)
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

function allocate_solar(mp::MicroProblem)
    (; solar_model, days, hours) = mp
    nmax = solar_model.wavelength_count
    nsteps = length(days) * length(hours)
    out = SolarRadiation.allocate_output_arrays(nsteps, length(days), nmax)
    buffers = SolarRadiation.allocate_buffers(nmax, solar_model.diffuse_model)
    return out, buffers
end

# In-place recompute: writes solar values directly into pre-allocated arrays
function solve_solar!(out::NamedTuple, buffers::NamedTuple, mp::MicroProblem)
    (; solar_model, days, hours, site) = mp
    SolarRadiation.solar_radiation!(out, buffers, solar_model;
        days, hours, solar_terrain=SolarRadiation.SolarTerrain(site),
    )
    # Cap zenith angles at 90° (sun-below-horizon sentinel). In-place broadcast
    # avoids the boolean-mask indexing path which allocates the mask vector.
    @. out.zenith_angle = min(out.zenith_angle, 90u"°")
    @. out.zenith_slope_angle = min(out.zenith_slope_angle, 90u"°")
    return out
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
    # adjust for cloud using Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992).
    # `solar_radiation_out.day_of_year` is already a per-step Vector — reuse it
    # instead of `repeat(days, inner=length(hours))` which allocates.
    zenith_angle = solar_radiation_out.zenith_angle
    direct_horizontal = solar_radiation_out.direct_horizontal
    diffuse_horizontal = solar_radiation_out.diffuse_horizontal
    cloud = output.cloud_cover
    doy = solar_radiation_out.day_of_year
    return cloud_adjust_radiation(output, cloud, diffuse_horizontal, direct_horizontal, zenith_angle, doy; diffuse_fraction_model)
end

# Solves soil temperature and moisture using pre-allocated cache buffers
function solve_soil!(cache::MicroCache)
    mp = cache.problem
    (; days, hours, depths, heights) = mp
    snow_model = mp.parameters.snow
    n_snow = n_snow_nodes(snow_model)
    snow_state = cache.state.snow_state
    soil_moisture = cache.state.soil_moisture
    nodes = cache.state.nodes
    nodes_day = cache.state.nodes_day
    ∑phase = cache.state.∑phase
    output = cache.output
    solar_radiation_out = cache.buffers.solar_out
    snow_scratch = cache.buffers.snow_scratch
    buffers = cache.buffers
    ode_integrator = cache.ode_integrator
    forcing = cache.forcing

    (; site, environment_minmax, environment_daily, initial_soil_temperature, initial_soil_moisture) = mp
    (; soil_thermal, soil_hydraulics) = mp.parameters
    soil_thermal_model = soil_thermal
    (; boundary_layer_model, vapour_pressure_equation, time_mode, convergence,
       hourly_rainfall, moist_step, maxpool) = mp.config
    moisture_mode = mp.config.soil_moisture_mode
    (; campbell_b_parameter, bulk_density, mineral_density, air_entry_water_potential) = soil_hydraulics
    init_soil_wetness!(moisture_mode)

    ndays = length(days)
    nhours = length(hours)
    num_soil_nodes = length(depths)
    # Type-level soil-node count for `ntuple(f, Val(N))` constructions below.
    Nsoil = cache.num_soil_nodes_val

    # initial conditions — T0 is always SVector{N_soil}
    # Default fallback (initial_soil_temperature=nothing): uniform daily mean of the
    # first day's reference air temperature, matching Fortran SINEC.f's TINS computation
    # (TIN = mean of the 25-element TAIRRY array; SINEC.f:124-138).
    if isnothing(initial_soil_temperature)
        t_init = _tins_mean_K(output.reference_temperature, 1, nhours)
        T0 = _filled_svec(t_init, Nsoil)
    else
        if num_soil_nodes != length(initial_soil_temperature)
            error("Initial soil temperature must match length of 'depths'")
        end
        T0 = SVector(ntuple(i -> initial_soil_temperature[i], Nsoil))
    end

    # Snow temperature — concrete SVector{N_snow} for SnowModel, `nothing` for
    # NoSnow. Dispatching on `snow_model` type avoids the Union{Nothing,SVector}
    # that a runtime ternary would produce.
    T_snow = init_snow_temperatures(snow_model, mp.initial_snow_temperature)

    # Reset pre-allocated scratch arrays
    nodes_day[1, 1:ndays] .= num_soil_nodes
    ∑phase .= 0.0u"J" # TODO decide whether this should happen and fix in Fortran
    soil_moisture .= initial_soil_moisture
    snow_state = initial_snow_state(snow_model, mp.initial_snow_depth, mp.initial_snow_density)
    reset_snow_scratch!(snow_model, snow_scratch)

    nodes .= nodes_day[:, 1]
    soil_prop_view = if n_snow > 0
        (; bulk_thermal_conductivity = view(buffers.soil_properties.bulk_thermal_conductivity, n_snow+1:n_snow+num_soil_nodes),
           bulk_heat_capacity = view(buffers.soil_properties.bulk_heat_capacity, n_snow+1:n_snow+num_soil_nodes),
           bulk_density = view(buffers.soil_properties.bulk_density, n_snow+1:n_snow+num_soil_nodes))
    else
        buffers.soil_properties
    end
    update_soil_properties!(output, soil_prop_view, soil_thermal_model;
        soil_temperature=T0, soil_moisture, bulk_density, mineral_density,
        atmospheric_pressure=101325.0u"Pa", step=1
    )

    soil_saturation_moisture = 1.0 .- bulk_density ./ mineral_density
    output.soil_water_potential[1, :] .= air_entry_water_potential .* (soil_saturation_moisture ./ soil_moisture) .^ campbell_b_parameter
    output.soil_temperature[1, :] .= T0
    output.soil_moisture[1, :] = soil_moisture

    initialise_soil_humidity!(moisture_mode, output, output.soil_water_potential[1, :], T0)

    environment_day = get_day(environment_daily, 1)
    environment_instant = get_instant(environment_day, mp.environment_hourly, output, soil_moisture, 1)
    longwave_sky = precompute_longwave_sky(; site, environment_instant, vapour_pressure_equation)

    # simulate all days
    pool = 0.0u"kg/m^2" # TODO make this an initialisation option
    qfreze = 0.0u"W/m^2"  # Fortran COMMON/melt/QFREZE: persists across hours TODO make it Q_freeze
    niter_moist = ustrip(u"s^-1", 3600 / moist_step)
    infil_out = nothing
    for j in 1:ndays
        iday = j
        nodes .= nodes_day[:, iday]
        environment_day = get_day(environment_daily, iday)
        update_forcing_day!(forcing, solar_radiation_out, output, iday)
        niter = iterations_for_day(time_mode, convergence, j)
        # Find solar midnight: noon is the hour with minimum zenith angle (sun highest),
        # midnight is 12 hours later. This correctly handles longitude offsets.
        # Solar noon = hour with minimum zenith angle. Manual scan avoids the
        # `arr[range]` + `ustrip.(...)` broadcast allocations.
        day_offset = (j - 1) * nhours
        noon_i = 1
        @inbounds best_zen = solar_radiation_out.zenith_angle[day_offset + 1]
        @inbounds for i in 2:nhours
            z = solar_radiation_out.zenith_angle[day_offset + i]
            if z < best_zen
                best_zen = z
                noon_i = i
            end
        end
        midnight_i = mod(noon_i - 1 + 12, nhours) + 1
        # Per-day state resets, split by what they touch so each TimeMode can opt in
        # independently. Phase reset matches Fortran MICROCLIMATE.f:1206-1209
        # (`qphase(:)=0; sumphase2(:)=0` every monthly-mode pass). Moisture/snow
        # resets are NonConsecutive-only — Fortran preserves these across passes
        # in monthly mode (verified: skipped OSUB returns at L353 before any
        # zeroing of moisture/snow accumulators).
        if reset_phase_per_iter(time_mode)
            ∑phase .= 0.0u"J"
        end
        if reset_moisture_per_day(time_mode)
            reset_day_soil_moisture!(moisture_mode, soil_moisture, initial_soil_moisture, j)
            snow_state = initial_snow_state(snow_model)
            reset_snow_scratch!(snow_model, snow_scratch)
        end
        if reset_snow_per_day(time_mode)
            snow_state = initial_snow_state(snow_model)
            reset_snow_scratch!(snow_model, snow_scratch)
        end
        # Soil-temperature reset (TINS-style uniform fill from this day's mean reference
        # air temperature). Matches Fortran MICROCLIMATE.f:909-916: `T(1:ii) = TINS(1:ii, doy)`
        # fires when `IFINAL == 1`. NonConsecutiveDayMode resets every day; ConsecutiveDayMode
        # never resets after pre-loop init.
        if is_reset_day(time_mode, j)
            sub2_start = (iday-1)*nhours + 1
            sub2_stop  = iday*nhours
            if isnothing(initial_soil_temperature)
                t_reset = _tins_mean_K(output.reference_temperature, sub2_start, sub2_stop)
                T0 = _filled_svec(t_reset, Nsoil)
            else
                T0 = SVector(ntuple(i -> initial_soil_temperature[i], Nsoil))
            end
            T_snow = init_snow_temperatures(snow_model, u"K"(0.0u"°C"))
        end
        T0 = setindex(T0, environment_instant.deep_soil_temperature, num_soil_nodes)
        use_multi_iter  = may_iterate(convergence)
        iter            = 0
        is_last_iter    = false
        T0_prev_start   = T0
        # Fortran MICROCLIMATE.f line 1203: T(I)=WORK(I+520) resets T to start-of-day
        # values after every SFODE call. Save initial conditions for reset.
        T0_day_start    = T0
        T_snow_day_start = T_snow
        ∑phase_day_start = copy(∑phase)  # carry frozen-water state from previous day
        while !is_last_iter
            iter += 1
            T0_iter_start = T0
            is_last_iter = is_converged(convergence, iter, niter, T0, T0_prev_start)
            T0_prev_start = T0_iter_start
            # ConsecutiveDayMode with spinup_first_day: each iter restarts day 1
            # from T0_day_start. NonConsecutiveDayMode (Fortran monthly): T continues
            # across spinup passes via the integrator's Y array, so don't reset.
            if iter_resets_T(time_mode)
                T0 = T0_day_start
                T_snow = T_snow_day_start
            end
            ∑phase .= ∑phase_day_start  # restore carry-over (zeros under reset_phase_per_iter)
            T0_before = T0
            T_snow_before = T_snow
            T0_output = T0  # phase-transition-clamped temperatures for output (see NOTE below)

            # ── Pre-fill day j's row 1 with the start-of-day state (NMR convention).
            # NMR records state at minute (i-1)*60 of each day, so row 1 of each day
            # is the *initial* state at minute 0. Hour i's integration result is then
            # written to row i+1 below (state at minute i*60 = start of hour i+1).
            day_init_step = (j - 1) * length(hours) + 1
            if is_last_iter
                # T0 here is the state at the start of the *output* pass — for
                # ConsecutiveDayMode that equals T0_day_start; for NonConsecutiveDayMode
                # it's the post-prior-spinup-pass state, which is what Fortran writes
                # at OUT(4) = T(1) at TIME=0 of the output pass.
                _write_row!(output.soil_temperature, day_init_step, T0)
                output.surface_water[day_init_step] = pool
                _write_row!(output.soil_moisture, day_init_step, soil_moisture)
                @inbounds for i in eachindex(soil_moisture)
                    output.soil_water_potential[day_init_step, i] =
                        air_entry_water_potential[i] * (soil_saturation_moisture[i] / soil_moisture[i]) ^ campbell_b_parameter[i]
                end
                initialise_soil_humidity!(moisture_mode, output,
                    view(output.soil_water_potential, day_init_step, :), T0)
                # Recompute longwave_sky from the forcing at minute 0 of day j so the
                # sky_temperature we record matches NMR's TIME=0 row exactly. The longwave_sky
                # variable in scope here may be stale from a previous iteration's last hour.
                env_init = get_instant(environment_day, mp.environment_hourly, output, soil_moisture, day_init_step)
                (; longwave_sky_init) = let
                    overrides_init = snow_surface_overrides(snow_model, snow_state, snow_scratch, day_init_step)
                    env_for_lw = env_init
                    if !isnothing(overrides_init.shade)
                        env_for_lw = merge(env_for_lw, (; shade=overrides_init.shade))
                    end
                    if !isnothing(overrides_init.emissivity)
                        env_for_lw = merge(env_for_lw, (; surface_emissivity=overrides_init.emissivity))
                    end
                    (; longwave_sky_init=precompute_longwave_sky(; site, environment_instant=env_for_lw, vapour_pressure_equation))
                end
                output.sky_temperature[day_init_step] = longwave_sky_init.sky_temperature
                update_soil_properties!(output, soil_prop_view, soil_thermal_model;
                    soil_temperature=T0, soil_moisture, bulk_density, mineral_density,
                    atmospheric_pressure=output.pressure[day_init_step],
                    step=day_init_step, vapour_pressure_equation,
                )
                if n_snow > 0
                    _write_row!(output.snow_temperature, day_init_step, T_snow)
                    output.snow_depth[day_init_step] = snow_state.current_depth
                    output.snow_density[day_init_step] = snow_state.density
                    output.snow_fall[day_init_step] = 0.0u"cm/hr"
                end
            end

            # ── Set up first hour's inputs and initialize integrator for full day ──
            step = (j - 1) * length(hours) + 1
            environment_instant = get_instant(environment_day, mp.environment_hourly, output, soil_moisture, step)
            T0 = setindex(T0, environment_instant.deep_soil_temperature, num_soil_nodes)
            (; albedo, effective_wetness, longwave_sky) = apply_snow_overrides(
                snow_model, snow_state, snow_scratch, step, site, moisture_mode, environment_instant, vapour_pressure_equation)
            if n_snow > 0
                compute_effective_depths!(snow_model, snow_state, snow_scratch, depths)
                ode_depths = snow_scratch.effective_depths
                first_active = findfirst(d -> d > 0u"cm", snow_scratch.node_depths)
                sync_temp = isnothing(first_active) ? T0[1] : T_snow[first_active]
                T_snow = _sync_snow_to_active(T_snow, first_active, sync_temp)
            else
                ode_depths = depths
            end
            T0_ode = combine_ode_state(T_snow, T0, n_snow)
            inputs = build_soil_energy_inputs(; forcing, buffers, soil_thermal_model,
                depths=ode_depths, heights, site, boundary_layer_model,
                n_snow, num_soil_nodes, nodes, environment_instant, effective_wetness,
                vapour_pressure_equation, longwave_sky, albedo, qfreze,
                snow_model, snow_state, snow_scratch, soil_moisture,
                bulk_density, mineral_density,
                maximum_surface_temperature=mp.config.maximum_surface_temperature,
            )

            # Initialize integrator for full day (0-1440 min), matching Fortran SFODE.
            # tspan is fixed at integrator init so we omit t0/tf — passing them
            # forces SciMLBase to rebuild the stage cache via `initialize!`,
            # which heap-allocates a `Memory{SVector}` per call.
            SciMLBase.reinit!(ode_integrator, T0_ode)
            ode_integrator.p = inputs

            for i in 1:length(hours)
                step = (j - 1) * length(hours) + i
                T0_before = T0
                T_snow_before = T_snow

                # ── Step integrator to this hour boundary ──
                T_result = step_to_hour!(ode_integrator, i)
                T_snow, T0 = split_ode_state(snow_model, T_result)

                # ── Phase transition (final iteration only) ── # TODO this should happen every time but at present it doesn't in Fortran version so keeping the same for now
                if is_last_iter
                    # `apply_phase_transition` mutates `∑phase` (passed as
                    # `accumulated_latent_heat`) in place and returns the new T.
                    T0 = apply_phase_transition(snow_model, T0, T0_before, buffers.phase_transition, ∑phase, soil_moisture, depths)
                    T0_output = T0
                end

                # ── Snow post-ODE: mass balance, clamping, node activation ──
                # Only on final iteration (Fortran OSUB returns before snow code on non-final iterations)
                snowfall_hour = 0.0u"cm"
                thermal_melt = 0.0u"cm"
                if n_snow > 0 && is_last_iter
                    # Evaluate evaporation at end-of-hour state (Fortran OSUB.f line 744)
                    Q_evap, _ = surface_convection_evaporation(;
                        surface_temperature=T0[1],
                        air_temperature=environment_instant.reference_temperature,
                        wind_speed=environment_instant.reference_wind_speed,
                        relative_humidity=environment_instant.reference_humidity,
                        atmospheric_pressure=environment_instant.atmospheric_pressure,
                        roughness_height=site.roughness_height,
                        reference_height=last(heights),
                        karman_constant=boundary_layer_model.karman_constant,
                        soil_wetness=effective_wetness,
                        vapour_pressure_equation,
                    )
                    snow_state, snowfall_hour, thermal_melt, rain_melt_snow = update_snow(snow_model, snow_state, snow_scratch, T_snow, T_snow_before, environment_instant, step,
                        T0[1], T0_before[1], Q_evap;
                        hourly_rainfall, shade=environment_instant.shade, day_of_year=days[j])
                    # Fortran OSUB.f lines 938-943: QFREZE = melted*(1-snowmelt)*79.7/60/10000
                    # Only applied when DOY > 1 (Fortran guard). Fortran units: cal/(min·cm²).
                    # melt depth (cm) × water density (kg/m³ → kg/m²/cm) × Lf (J/kg) / 1 hr → W/m²
                    melt_factor = snow_model.snow_melt_factor
                    if days[j] > 1 && melt_factor <= 1.0
                        qfreze = uconvert(u"W/m^2", thermal_melt * water_density * LATENT_HEAT_FUSION / 1u"hr" * (1.0 - melt_factor))
                    else
                        qfreze = 0.0u"W/m^2"
                    end
                    # Fortran SNOWLAYER.f lines 67-71, 93-97: reset deactivated node temps to T(1)
                    prev_active = snow_state.active_nodes
                    snow_state = activate_snow_nodes(snow_model, snow_state, snow_scratch, T_snow, step)
                    if snow_state.active_nodes < prev_active
                        T_snow = sync_inactive_snow_temps(T_snow, snow_scratch)
                    end                    
                    snow_state, T_snow, soil_surf_temp = snow_phase_transition(snow_model, snow_state, snow_scratch, T_snow, T_snow_before, environment_instant.atmospheric_pressure, T0[1])
                    if !isnothing(soil_surf_temp)
                        T0 = setindex(T0, soil_surf_temp, 1)
                    end
                    T_snow, clamp_soil = clamp_snow_temperatures(T_snow, snow_scratch)
                    if clamp_soil && T0[1] > u"K"(0.0u"°C")
                        T0 = setindex(T0, u"K"(0.0u"°C"), 1)
                    end
                    T_snow, freeze_soil = freeze_new_snow(snow_model, T_snow, snow_scratch, step)
                    if freeze_soil && T0[1] > u"K"(0.0u"°C")
                        T0 = setindex(T0, u"K"(0.0u"°C"), 1)
                    end
                    snow_state = adjust_snow_near_nodes(snow_model, snow_state, snow_scratch, T_snow, step)
                    prev_active2 = snow_state.active_nodes
                    snow_state = activate_snow_nodes(snow_model, snow_state, snow_scratch, T_snow, step)
                    # Fortran SNOWLAYER.f: also reset temps after the second activation
                    if snow_state.active_nodes < prev_active2
                        T_snow = sync_inactive_snow_temps(T_snow, snow_scratch)
                    end
                end

                # ── Feed modified state back into integrator for next hour ──
                # Matches Fortran OSUB.f line 1087: t=tt (copies modified temps back to ODE state)
                if i < length(hours)
                    next_step = step + 1
                    environment_instant = get_instant(environment_day, mp.environment_hourly, output, soil_moisture, next_step)
                    T0 = setindex(T0, environment_instant.deep_soil_temperature, num_soil_nodes)
                    (; albedo, effective_wetness, longwave_sky) = apply_snow_overrides(
                        snow_model, snow_state, snow_scratch, next_step, site, moisture_mode, environment_instant, vapour_pressure_equation)
                    if n_snow > 0
                        compute_effective_depths!(snow_model, snow_state, snow_scratch, depths)
                        ode_depths = snow_scratch.effective_depths
                    else
                        ode_depths = depths
                    end
                    ode_integrator.p = build_soil_energy_inputs(; forcing, buffers, soil_thermal_model,
                        depths=ode_depths, heights, site, boundary_layer_model,
                        n_snow, num_soil_nodes, nodes, environment_instant, effective_wetness,
                        vapour_pressure_equation, longwave_sky, albedo, qfreze,
                        snow_model, snow_state, snow_scratch, soil_moisture,
                        bulk_density, mineral_density,
                        maximum_surface_temperature=mp.config.maximum_surface_temperature,
                    )
                end
                T0_ode = combine_ode_state(T_snow, T0, n_snow)
                # Direct field write bypasses any custom setproperty! interposing
                # an SVector copy. The integrator's `.u` field is concretely
                # typed `SVector{N,T}` so this is a plain memory write.
                setfield!(ode_integrator, :u, T0_ode)
                SciMLBase.u_modified!(ode_integrator, true)

                init_soil_obukhov!(buffers, forcing, site, boundary_layer_model, heights, T0, i)
                # Fortran OSUB.f: soil moisture runs only on final iteration (after line 353 guard)
                if is_last_iter
                    rain = if hourly_rainfall
                        mp.environment_hourly.rainfall[step]
                    elseif i == midnight_i # daily rainfall falls at solar midnight
                        environment_instant.rainfall
                    else
                        0.0u"kg/m^2"
                    end
                    # Fortran OSUB.f: apply rainmult to rainfall entering condep
                    if n_snow > 0
                        rain = rain * snow_model.rain_multiplier
                    end
                    # Pool update:
                    # Cold snow (air_temp < snow_threshold): rain → snowfall (handled in update_snow); only thermal melt enters pool.
                    # Warm snow or no snow: rain passes through pack to soil surface, rain_melt snow water also drains to surface.
                    #   Fortran OSUB.f line 1116 zeros rainfallb when cursnow > 0 for non-hourly rain (mass conservation
                    #   violation — rain water discarded). Julia instead routes rain to pool, and also tracks rain_melt water
                    #   (snow melted by rain energy; OSUB.f line 950 subtracts rainmelt from cursnow but never adds it to condep
                    #   — another mass loss). Both rain and rain_melt water physically reach the soil surface.
                    snow_present = n_snow > 0 && snow_scratch.snow_depth_hourly[step] > 0.0u"cm"
                    snow_temp_threshold = n_snow > 0 ? snow_model.snow_temperature_threshold : 0.0u"°C"
                    melted_water = n_snow > 0 ? uconvert(u"kg/m^2", thermal_melt * 1.0u"g/cm^3") : 0.0u"kg/m^2"
                    rain_melt_water = n_snow > 0 ? uconvert(u"kg/m^2", rain_melt_snow * snow_state.density) : 0.0u"kg/m^2"
                    if snow_present && u"°C"(environment_instant.reference_temperature) < snow_temp_threshold
                        # Cold snow: rain absorbed into snowpack (handled in update_snow as snowfall); only thermal melt enters pool
                        pool = clamp(pool + melted_water, 0.0u"kg/m^2", maxpool)
                    else
                        # Warm snow or no snow: rain + rain_melt water + thermal melt enter pool
                        pool = clamp(pool + rain + rain_melt_water + melted_water, 0.0u"kg/m^2", maxpool)
                    end
                    # Soil moisture physics; output write happens below at output_step
                    (; pool, soil_moisture, infil_out) = step_soil_moisture!(moisture_mode, buffers, soil_hydraulics;
                        depths, site, boundary_layer_model, environment_instant, T0, niter_moist, pool, soil_moisture,
                        moist_step, moist_error=mp.config.moist_error, moist_count=mp.config.moist_count, maxpool,
                        vapour_pressure_equation, snow_present,
                    )
                end
                # Write hour i's result to step + 1 (NMR convention: state at minute i*60
                # belongs in row (j-1)*nhours + i + 1 = start of hour i+1). Skip when
                # step+1 would land at row 1 of a future day that is going to reset T —
                # the next day's pre-fill will overwrite that slot anyway.
                output_step = step + 1
                in_bounds = output_step <= ndays * length(hours)
                crosses_day = (i == length(hours))
                next_day_resets = crosses_day && j < ndays && is_reset_day(time_mode, j + 1)
                write_output = is_last_iter && in_bounds && !next_day_resets
                if write_output
                    output.surface_water[output_step] = pool
                    _write_row!(output.soil_temperature, output_step, T0_output)
                    output.sky_temperature[output_step] = longwave_sky.sky_temperature
                    if n_snow > 0
                        output.snow_fall[output_step] = uconvert(u"cm/hr", snowfall_hour / 1.0u"hr")
                        output.snow_depth[output_step] = uconvert(u"cm", snow_scratch.snow_depth_hourly[step])
                        output.snow_density[output_step] = uconvert(u"g/cm^3", snow_state.density)
                        _write_row!(output.snow_temperature, output_step, T_snow)
                    end
                    if !isnothing(infil_out)
                        update_soil_water!(output, infil_out, output_step)
                    end
                    environment_instant = get_instant(environment_day, mp.environment_hourly, output, soil_moisture, output_step)
                    update_soil_properties!(output, soil_prop_view, soil_thermal_model;
                        soil_temperature=T0, soil_moisture, bulk_density, mineral_density,
                        atmospheric_pressure=environment_instant.atmospheric_pressure, step=output_step, vapour_pressure_equation
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
    solar_radiation_out = cache.buffers.solar_out
    site = mp.site
    (; boundary_layer_model, vapour_pressure_equation) = mp.config
    profile_buffers = cache.buffers.air_profile
    for i in 1:size(output.profile.air_temperature, 1)
        surface_temperature = u"°C"(output.soil_temperature[i, 1])
        environment_instant = (;
            atmospheric_pressure=output.pressure[i],
            reference_temperature=output.reference_temperature[i],
            reference_wind_speed=output.reference_wind_speed[i],
            reference_humidity=output.reference_humidity[i],
            zenith_angle=solar_radiation_out.zenith_angle[i],
        )
        result = atmospheric_surface_profile!(profile_buffers; site, boundary_layer_model, environment_instant, surface_temperature, vapour_pressure_equation)
        _write_row!(output.profile.air_temperature,   i, result.air_temperature)
        _write_row!(output.profile.wind_speed,        i, result.wind_speed)
        _write_row!(output.profile.relative_humidity, i, result.relative_humidity)
        output.profile.convective_heat_flux[i]  = result.convective_heat_flux
        output.profile.friction_velocity[i]     = result.friction_velocity
    end
end

function allocate_ode_integrator(T0, inputs_proto, solver, solver_kwargs)
    # Full-day integration (0-1440 min), matching Fortran SFODE which integrates
    # continuously. Hourly bookkeeping is done via step!/add_tstop! in the day loop.
    tspan = (0.0u"minute", 1440.0u"minute")
    prob = ODEProblem{false}(soil_energy_balance, T0, tspan, inputs_proto)
    return SciMLBase.init(prob, solver; save_everystep=false, save_start=false, save_end=false, solver_kwargs...)
end

# Step the integrator to an exact hour boundary and return the state.
function step_to_hour!(integrator, hour_index)
    target_t = hour_index * 60.0u"minute"
    SciMLBase.add_tstop!(integrator, target_t)
    while integrator.t < target_t - 1e-10u"minute"
        SciMLBase.step!(integrator)
    end
    return integrator.u
end

function update_soil_properties!(output, soil_properties_buffers, soil_thermal_model; step, kw...)
    (; bulk_thermal_conductivity, bulk_heat_capacity, bulk_density) = soil_properties!(soil_properties_buffers, soil_thermal_model; kw...)

    _write_row!(output.soil_thermal_conductivity, step, bulk_thermal_conductivity)
    _write_row!(output.soil_heat_capacity,        step, bulk_heat_capacity)
    _write_row!(output.soil_bulk_density,         step, bulk_density)

    return output
end

# Explicit-loop row write into a Matrix. Avoids the SubArray + Broadcasted
# heap allocations from `arr[step, :] .= vec` since matrices are column-major.
@inline function _write_row!(arr::AbstractMatrix, step::Int, src::AbstractVector)
    @inbounds for i in eachindex(src)
        arr[step, i] = src[i]
    end
    return arr
end

function update_soil_water!(output, infil_out, step)
    _write_row!(output.soil_moisture, step, infil_out.soil_moisture)
    _write_row!(output.soil_water_potential, step, infil_out.soil_water_potential)
    _write_row!(output.soil_humidity, step, infil_out.soil_humidity)
    return output
end

# ── Soil moisture stepping dispatch ───────────────────────────────────────

function step_soil_moisture!(mode::DynamicSoilMoisture, buffers, soil_hydraulics;
    depths, site, boundary_layer_model, environment_instant, T0, niter_moist, pool, soil_moisture,
    moist_step, moist_error, moist_count, maxpool,
    vapour_pressure_equation, snow_present=false,
)
    (; infil_out, soil_wetness, pool, soil_moisture) = get_soil_water_balance!(buffers, soil_hydraulics;
        depths, site, boundary_layer_model, environment_instant, T0, niter_moist, pool,
        soil_wetness=mode.soil_wetness, soil_moisture,
        moist_step, moist_error, moist_count, maxpool,
        vapour_pressure_equation, snow_present,
    )
    mode.soil_wetness = soil_wetness
    return (; pool, soil_moisture, infil_out)
end

step_soil_moisture!(::PrescribedSoilMoisture, args...; pool, soil_moisture, kw...) = (; pool, soil_moisture, infil_out=nothing)

"""
Build a single `MicroForcing` whose 8 `ScaledInterpolation`s wrap freshly-allocated
25-element coefficient buffers. Built ONCE in `init` and reused across every
day-iter via `update_forcing_day!` which only mutates the underlying `itp.coefs`.

Fortran MICROCLIMATE.f lines 843-852: 25 interpolation nodes (0,60,...,1440).
The 25th node repeats the 24th hour's value (TD(36)=TAIRhr1(DOYF)).
"""
function allocate_forcing(output, solar_radiation_out)
    tspan = 0.0:60:(60*24)  # 0, 60, ..., 1440
    _scaled(buf) = scale(interpolate(buf, BSpline(Linear())), tspan)
    interpolate_solar         = _scaled(zeros(eltype(output.global_radiation), 25))
    interpolate_zenith        = _scaled(zeros(eltype(solar_radiation_out.zenith_angle), 25))
    interpolate_slope_zenith  = _scaled(zeros(eltype(solar_radiation_out.zenith_slope_angle), 25))
    interpolate_temperature   = _scaled(zeros(typeof(1.0u"K"), 25))
    interpolate_wind          = _scaled(zeros(eltype(output.reference_wind_speed), 25))
    interpolate_humidity      = _scaled(zeros(eltype(output.reference_humidity), 25))
    interpolate_cloud         = _scaled(zeros(eltype(output.cloud_cover), 25))
    interpolate_pressure      = _scaled(zeros(eltype(output.pressure), 25))
    return MicroForcing(;
        interpolate_solar, interpolate_zenith, interpolate_slope_zenith,
        interpolate_temperature, interpolate_wind, interpolate_humidity,
        interpolate_cloud, interpolate_pressure,
    )
end

@inline function _copy_day_into!(coefs::AbstractVector, src::AbstractVector, sub1::AbstractRange)
    @inbounds for (k, i) in enumerate(sub1)
        coefs[k] = src[i]
    end
    @inbounds coefs[25] = src[last(sub1)]  # 25th node repeats hour 24
    return nothing
end

@inline function _copy_day_into!(coefs::AbstractVector, src::AbstractVector, sub1::AbstractRange, ::typeof(uconvert), unit)
    @inbounds for (k, i) in enumerate(sub1)
        coefs[k] = uconvert(unit, src[i])
    end
    @inbounds coefs[25] = uconvert(unit, src[last(sub1)])
    return nothing
end

"""
Update the in-place `MicroForcing` for day `iday` by writing into each
`itp.coefs` buffer. No allocation.
"""
function update_forcing_day!(forcing::MicroForcing, solar_radiation_out, output, iday::Int)
    nhours = 24
    sub1 = (iday*nhours - nhours + 1):(iday*nhours)
    # `forcing.interpolate_*.itp.coefs` is the 25-element buffer underneath the
    # ScaledInterpolation, mutated in place.
    _copy_day_into!(forcing.interpolate_solar.itp.coefs,        output.global_radiation,            sub1)
    _copy_day_into!(forcing.interpolate_zenith.itp.coefs,       solar_radiation_out.zenith_angle,   sub1)
    _copy_day_into!(forcing.interpolate_slope_zenith.itp.coefs, solar_radiation_out.zenith_slope_angle, sub1)
    _copy_day_into!(forcing.interpolate_temperature.itp.coefs,  output.reference_temperature,       sub1, uconvert, u"K")
    _copy_day_into!(forcing.interpolate_wind.itp.coefs,         output.reference_wind_speed,        sub1)
    _copy_day_into!(forcing.interpolate_humidity.itp.coefs,     output.reference_humidity,          sub1)
    _copy_day_into!(forcing.interpolate_cloud.itp.coefs,        output.cloud_cover,                 sub1)
    _copy_day_into!(forcing.interpolate_pressure.itp.coefs,     output.pressure,                    sub1)
    return forcing
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
