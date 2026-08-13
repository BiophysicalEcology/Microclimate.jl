# Output quantities are written straight into the result buffers; any other
# quantity in `forcings` (e.g. vapour pressure feeding a derived humidity) is an
# intermediate and needs a scratch buffer. Allocate one per intermediate.
allocate_interpolation_scratch(::Nothing) = nothing
function allocate_interpolation_scratch(env)
    fs = env.forcings
    names = filter(k -> !_is_output_quantity(Val(k)), keys(fs))
    return NamedTuple{names}(map(k -> _quantity_buffer(getproperty(fs, k)), names))
end

# Quantity names the solver writes directly into its result buffers.
@inline _is_output_quantity(::Val) = false
@inline _is_output_quantity(::Val{:reference_temperature}) = true
@inline _is_output_quantity(::Val{:reference_humidity}) = true
@inline _is_output_quantity(::Val{:reference_wind_speed}) = true
@inline _is_output_quantity(::Val{:cloud_cover}) = true

_quantity_buffer(f::DielForcing) = zeros(eltype(first(f.values)), 24 * n_days(f))

function populate_weather!(output, mp, solar_radiation_out, interpolation_scratch)
    (; hours) = mp.model
    days = mp.days
    shortwave_model = mp.model.radiation.shortwave_model
    (; environment_minmax, environment_daily, environment_hourly) = mp.inputs
    interpolate_minmax!(output, environment_minmax, environment_daily, environment_hourly, solar_radiation_out, interpolation_scratch)
    (; global_radiation, diffuse_fraction) = adjust_for_cloud_cover!(
        output, solar_radiation_out, days, hours, shortwave_model,
    )
    output.diffuse_fraction .= diffuse_fraction
    if !isnothing(environment_hourly) && !isnothing(environment_hourly.global_radiation)
        output.global_radiation .= environment_hourly.global_radiation
    else
        output.global_radiation .= global_radiation
    end
    return output
end

# Apply snow surface overrides (albedo, wetness, emissivity) and compute the
# longwave sky terms in one place. Returns concrete Float64 / NamedTuple types
# regardless of whether snow is active so that consumers (in particular the
# integrator's `p` field) see stable types.
#
# PR #102 / Fortran OSUB.f:354 (`IF IFINAL .LT. ND GO TO 200`) skips the snow
# block on non-final spinup iters — overrides are applied only on the final
# (output) iter. We mirror that via the `is_last_iter` flag.
@inline function apply_snow_overrides(::NoSnow, snow_state, snow_scratch, step, site, moisture_mode, environment_instant, longwave_model, vapour_pressure_equation; is_last_iter::Bool=true)
    longwave_sky = precompute_longwave_sky(longwave_model;
        site, environment_instant, vapour_pressure_equation)
    return (;
        albedo = _baseline_albedo(site, environment_instant),
        effective_wetness = get_soil_wetness(moisture_mode, environment_instant),
        longwave_sky,
    )
end

@inline function apply_snow_overrides(snow_model::SnowModel, snow_state, snow_scratch, step,
    site, moisture_mode, environment_instant, longwave_model, vapour_pressure_equation;
    is_last_iter::Bool=true,
)
    # Fortran OSUB.f line 1023-1041: snow surface overrides fire only when the
    # previous hour had snow on the ground (prev_depth > 0). Additionally, on
    # spinup iters (IFINAL < ND) snow code is skipped entirely.
    snow_active = is_last_iter && (
        (step > 1 ? snow_scratch.snow_depth_hourly[step - 1] : snow_state.current_depth) > 0.0u"cm"
    )
    if snow_active
        albedo = snow_albedo(snow_state.days_since_snow)
        effective_wetness = 1.0
        longwave_sky = precompute_longwave_sky(longwave_model;
            site, environment_instant, vapour_pressure_equation,
            surface_emissivity = 0.98,
        )
    else
        albedo = _baseline_albedo(site, environment_instant)
        effective_wetness = get_soil_wetness(moisture_mode, environment_instant)
        longwave_sky = precompute_longwave_sky(longwave_model;
            site, environment_instant, vapour_pressure_equation,
        )
    end
    return (; albedo, effective_wetness, longwave_sky)
end

# Resolve the non-snow baseline albedo: a populated `environment_daily.albedo`
# (surfaced here as `environment_instant.albedo`) wins; otherwise fall back
# to the `Site.albedo` scalar. The inner helper has two methods dispatched
# on whether the env value is `Nothing`, so the return type stays concrete
# per call site (no `Union{Float64, Nothing}` leakage).
@inline _baseline_albedo(site, environment_instant) =
    _pick_albedo(site.albedo, environment_instant.albedo)
@inline _pick_albedo(site_albedo, ::Nothing)   = site_albedo
@inline _pick_albedo(_site_albedo, env_albedo) = env_albedo

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
    (; hours, depths, heights,
       soil_properties_model, soil_hydraulic_model, snow_model,
       vapour_pressure_equation, boundary_layer_model,
       radiation, evaporation_model, soil_energy_model,
       config) = mp.model
    days = mp.days
    longwave_model = radiation.longwave_model
    (; site, soil_profile, environment_minmax, environment_daily, environment_hourly,
       initial_soil_temperature, initial_soil_moisture,
       initial_snow_depth, initial_snow_temperature, initial_snow_density) = mp.inputs
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
    snow_state = initial_snow_state(snow_model, initial_snow_depth, initial_snow_density)
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
    buffers = MicroBuffers(
        solar_radiation_out,
        solar_buffers,
        allocate_profile(boundary_layer_model, heights),
        allocate_profile(boundary_layer_model, heights),
        allocate_soil_energy_balance(soil_energy_model, num_ode_nodes),
        allocate_soil_properties(zeros(num_ode_nodes), soil_properties_model, soil_profile),
        allocate_phase_transition(soil_energy_model.freezing_model, num_soil_nodes),
        allocate_soil_water_balance(soil_hydraulic_model, num_soil_nodes),
        snow_scratch,
        allocate_interpolation_scratch(environment_minmax),
    )

    # ODE integrator — sized for the combined (snow + soil) vector
    soil_moisture = collect(initial_soil_moisture)
    ∑phase = zeros(typeof(1.0u"J"), num_soil_nodes)
    T0_soil = if isnothing(initial_soil_temperature)
        SVector(ntuple(_ -> 280.0u"K", num_soil_nodes))
    else
        SVector(ntuple(i -> initial_soil_temperature[i], num_soil_nodes))
    end
    T0_snow = init_snow_temperatures(snow_model, initial_snow_temperature)
    T0_ode = combine_ode_state(T0_snow, T0_soil, n_snow)

    # Build prototype depths for ODE (combined snow + soil)
    if n_snow > 0
        depths_proto = collect(combined_depths(snow_model, snow_state, snow_scratch, depths))
    else
        depths_proto = collect(depths)
    end

    populate_weather!(output, mp, solar_radiation_out, buffers.interpolation)

    # Pre-allocate the day-forcing interpolations once. The 8 ScaledInterpolations
    # wrap 25-element coefficient buffers that get refilled per day-iter via
    # `update_forcing_day!` — no per-day construction.
    forcing = allocate_forcing(output, solar_radiation_out)
    update_forcing_day!(forcing, solar_radiation_out, output, 1)

    environment_day = get_day(environment_daily, 1)
    environment_instant = get_instant(environment_day, environment_hourly, output, soil_moisture, 1)
    longwave_sky = precompute_longwave_sky(longwave_model;
        site, environment_instant, vapour_pressure_equation)
    inputs_proto = SoilEnergyInputs(;
        model=soil_energy_model,
        evaporation_model,
        atmospheric_radiation_model=longwave_model.atmospheric_radiation_model,
        forcing,
        buffers, soil_properties_model, soil_profile,
        depths=depths_proto, heights,
        site, boundary_layer_model,
        environment_instant,
        soil_wetness=0.0, vapour_pressure_equation,
        longwave_sky, albedo=site.albedo,
        soil_moisture,
        snow_model, snow_state=initial_snow_state(snow_model),
    )
    ode_integrator = allocate_ode_integrator(soil_energy_model, T0_ode, inputs_proto)

    state = MicroState(soil_moisture, ∑phase, similar(∑phase))

    return MicroCache(mp, output, state, buffers, forcing, ode_integrator, Val(num_soil_nodes))
end

function CommonSolve.solve!(cache::MicroCache)
    mp = cache.problem
    output = cache.output
    solar_radiation_out = cache.buffers.solar_out

    # Recompute solar in place into the pre-allocated cache buffers.
    solve_solar!(solar_radiation_out, cache.buffers.solar, mp)

    populate_weather!(output, mp, solar_radiation_out, cache.buffers.interpolation)
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
    reinit!(cache::MicroCache, inputs::MicroInputs)

Swap the inputs on `cache.problem` while keeping the model and every
pre-allocated array. The model is constant across runs, so dimensions are
guaranteed to match — invalid input sizes will surface as a `BoundsError`
or unit-mismatch on first access during `solve!`.
"""
function reinit!(cache::MicroCache, inputs::MicroInputs)
    cache.problem = MicroProblem(cache.problem.days, cache.problem.time_mode, cache.problem.model, inputs)
    return cache
end

function allocate_solar(mp::MicroProblem)
    (; hours) = mp.model
    days = mp.days
    solar_radiation_model = mp.model.radiation.solar_radiation_model
    nmax = solar_radiation_model.wavelength_count
    nsteps = length(days) * length(hours)
    out = SolarRadiation.allocate_output_arrays(nsteps, length(days), nmax)
    buffers = SolarRadiation.allocate_buffers(nmax, solar_radiation_model.diffuse_model)
    return out, buffers
end

# In-place recompute: writes solar values directly into pre-allocated arrays
function solve_solar!(out::NamedTuple, buffers::NamedTuple, mp::MicroProblem)
    (; hours) = mp.model
    days = mp.days
    solar_radiation_model = mp.model.radiation.solar_radiation_model
    site = mp.inputs.site
    SolarRadiation.solar_radiation!(out, buffers, solar_radiation_model;
        days, hours, solar_terrain=SolarRadiation.SolarTerrain(site),
    )
    # Cap zenith angles at 90° (sun-below-horizon sentinel). In-place broadcast
    # avoids the boolean-mask indexing path which allocates the mask vector.
    @. out.zenith_angle = min(out.zenith_angle, 90u"°")
    @. out.zenith_slope_angle = min(out.zenith_slope_angle, 90u"°")
    return out
end

function interpolate_minmax!(output, environment_minmax, environment_daily, environment_hourly, solar_radiation_out, interpolation_scratch)
    # Expand each forcing's diel curve to hourly, writing directly into the
    # matching output buffer.
    # Fill every quantity buffer: output quantities into the result, intermediates
    # into scratch. Derived quantities (e.g. humidity from vapour pressure) read
    # the forcings filled before them.
    buffer_for(vk) = _is_output_quantity(vk) ? getproperty(output, _name(vk)) :
                                               getproperty(interpolation_scratch, _name(vk))
    populate_quantities!(buffer_for, environment_minmax.forcings, solar_radiation_out, consecutive_days(environment_minmax))

    # Clamp humidity and cloud cover to [_, 1.0] without allocating a mask Vector.
    @inbounds for i in eachindex(output.reference_humidity)
        if output.reference_humidity[i] > 1.0
            output.reference_humidity[i] = 1.0
        end
    end
    @inbounds for i in eachindex(output.cloud_cover)
        if output.cloud_cover[i] > 1.0
            output.cloud_cover[i] = 1.0
        end
    end

    output.pressure .= environment_hourly.pressure
    return nothing
end
function interpolate_minmax!(output, environment_minmax::Nothing, environment_daily, environment_hourly, solar_radiation_out, _)
    output.cloud_cover .= environment_hourly.cloud_cover
    output.reference_temperature .= environment_hourly.reference_temperature
    output.pressure .= environment_hourly.pressure
    output.reference_wind_speed .= environment_hourly.reference_wind_speed
    output.reference_humidity .= environment_hourly.reference_humidity

    return 
end

function adjust_for_cloud_cover!(output, solar_radiation_out, days, hours,
    shortwave_model::AbstractShortwaveModel=AngstromMaxwellShortwave(),
)
    # `solar_radiation_out.day_of_year` is already a per-step Vector — reuse it
    # instead of `repeat(days, inner=length(hours))` which allocates.
    zenith_angle = solar_radiation_out.zenith_angle
    direct_horizontal = solar_radiation_out.direct_horizontal
    diffuse_horizontal = solar_radiation_out.diffuse_horizontal
    cloud = output.cloud_cover
    doy = solar_radiation_out.day_of_year
    return shortwave_radiation!(shortwave_model, output, cloud, diffuse_horizontal, direct_horizontal, zenith_angle, doy)
end

# Solves soil temperature and moisture using pre-allocated cache buffers
function solve_soil!(cache::MicroCache)
    mp = cache.problem
    (; hours, depths, heights,
       soil_properties_model, soil_hydraulic_model, snow_model,
       vapour_pressure_equation, boundary_layer_model,
       radiation, evaporation_model, soil_energy_model,
       config) = mp.model
    days = mp.days
    longwave_model = radiation.longwave_model
    soil_freezing_model = soil_energy_model.freezing_model
    (; site, soil_profile, environment_daily, environment_hourly,
       initial_soil_temperature, initial_soil_moisture,
       initial_snow_depth, initial_snow_temperature, initial_snow_density) = mp.inputs
    n_snow = n_snow_nodes(snow_model)
    soil_moisture = cache.state.soil_moisture
    ∑phase = cache.state.∑phase
    ∑phase_day_start = cache.state.∑phase_day_start
    output = cache.output
    buffers = cache.buffers
    solar_radiation_out = buffers.solar_out
    ode_integrator = cache.ode_integrator
    forcing = cache.forcing

    (; convergence, rainfall_schedule, max_surface_pool) = config
    time_mode = mp.time_mode
    moisture_mode = config.soil_moisture_strategy
    (; campbell_b_parameter, air_entry_water_potential) = soil_profile.hydraulics
    (; bulk_density, mineral_density, mineral_conductivity, mineral_heat_capacity) = soil_profile
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
    T_snow = init_snow_temperatures(snow_model, initial_snow_temperature)

    # Reset pre-allocated scratch arrays
    ∑phase .= 0.0u"J" # TODO decide whether this should happen and fix in Fortran
    soil_moisture .= initial_soil_moisture
    snow_state = initial_snow_state(snow_model, initial_snow_depth, initial_snow_density)
    reset_snow_scratch!(snow_model, buffers.snow)

    soil_prop_view = if n_snow > 0
        (; bulk_thermal_conductivity = view(buffers.soil_properties.bulk_thermal_conductivity, n_snow+1:n_snow+num_soil_nodes),
           bulk_heat_capacity = view(buffers.soil_properties.bulk_heat_capacity, n_snow+1:n_snow+num_soil_nodes),
           bulk_density = view(buffers.soil_properties.bulk_density, n_snow+1:n_snow+num_soil_nodes))
    else
        buffers.soil_properties
    end
    update_soil_properties!(output, soil_prop_view, soil_properties_model;
        soil_temperature=T0, soil_moisture, bulk_density, mineral_density,
        mineral_conductivity, mineral_heat_capacity,
        atmospheric_pressure=101325.0u"Pa", step=1
    )

    # Write soil saturation moisture and initial soil water potential without
    # allocating intermediate broadcast Vectors. `saturation_water_content` is
    # the existing scratch buffer used by `soil_water_balance!`; its first
    # `num_soil_nodes` entries hold the same per-depth value (the +1 entry is
    # the deep boundary, set elsewhere).
    soil_saturation_moisture = buffers.soil_water_balance.saturation_water_content
    @inbounds for i in eachindex(soil_moisture)
        soil_saturation_moisture[i] = 1.0 - bulk_density[i] / mineral_density[i]
        output.soil_water_potential[1, i] =
            air_entry_water_potential[i] * (soil_saturation_moisture[i] / soil_moisture[i]) ^ campbell_b_parameter[i]
    end
    _write_row!(output.soil_temperature, 1, T0)
    _write_row!(output.soil_moisture, 1, soil_moisture)

    initialise_soil_humidity!(moisture_mode, output,
        view(output.soil_water_potential, 1, :), T0)

    environment_day = get_day(environment_daily, 1)
    environment_instant = get_instant(environment_day, environment_hourly, output, soil_moisture, 1)
    longwave_sky = precompute_longwave_sky(longwave_model; site, environment_instant, vapour_pressure_equation)

    # simulate all days
    pool = 0.0u"kg/m^2" # TODO make this an initialisation option
    Q_freeze = 0.0u"W/m^2"  # Fortran COMMON/melt/QFREZE: persists across hours
    infil_out = nothing
    for j in 1:ndays
        iday = j
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
            # PR #102: Fortran MICROCLIMATE.f:1206-1209 resets `qphase(:)=0` for
            # *all* nodes, including snow nodes — the within-day latent-heat
            # budget has no meaning across non-consecutive representative days.
            if !isnothing(snow_state)
                snow_state = setproperties(snow_state, (; sum_phase=0.0u"J/m^2"))
            end
        end
        if reset_moisture_per_day(time_mode)
            reset_day_soil_moisture!(moisture_mode, soil_moisture, initial_soil_moisture, j)
        end
        if reset_snow_per_day(time_mode)
            snow_state = initial_snow_state(snow_model)
            reset_snow_scratch!(snow_model, buffers.snow)
        end
        # Soil-temperature reset (TINS-style uniform fill from this day's mean reference
        # air temperature). Matches Fortran MICROCLIMATE.f:909-916: `T(1:ii) = TINS(1:ii, doy)`
        # fires when `IFINAL == 1`. NonConsecutiveDayMode resets every day; ConsecutiveDayMode
        # never resets after pre-loop init.
        if is_reset_day(time_mode, j)
            sub2_start = (iday-1)*nhours + 1
            sub2_stop  = iday*nhours
            # Compute TINS unconditionally so we can also use it for the snow-
            # node temperature reset. Matches Fortran SINEC.f line 130 which
            # writes TINS for every monthly day regardless of initial_soil_temperature.
            t_reset = _tins_mean_K(output.reference_temperature, sub2_start, sub2_stop)
            if isnothing(initial_soil_temperature)
                T0 = _filled_svec(t_reset, Nsoil)
            else
                T0 = SVector(ntuple(i -> initial_soil_temperature[i], Nsoil))
            end
            # PR #102: snow nodes also reset to the day's TINS mean (not 0°C),
            # so the spinup starts from the same near-air-temperature initial
            # condition that Fortran uses.
            T_snow = init_snow_temperatures(snow_model, t_reset)
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
        # Snapshot day-start ∑phase into pre-allocated buffer (no copy alloc).
        @inbounds for i in eachindex(∑phase)
            ∑phase_day_start[i] = ∑phase[i]
        end
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
            # PR #102: also reset snow sum_phase at start of each iter when the
            # mode resets phase per-iter — Fortran resets `qphase` for all nodes.
            if !isnothing(snow_state) && reset_phase_per_iter(time_mode)
                snow_state = setproperties(snow_state, (; sum_phase=0.0u"J/m^2"))
            end
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
                env_init = get_instant(environment_day, environment_hourly, output, soil_moisture, day_init_step)
                # Same Fortran-aligned snow override (snow → emissivity 0.98) but
                # without the type-unstable conditional `merge` — use the
                # apply_snow_overrides helper from above.
                (; longwave_sky) = apply_snow_overrides(snow_model, snow_state, buffers.snow,
                    day_init_step, site, moisture_mode, env_init,
                    longwave_model, vapour_pressure_equation)
                output.sky_temperature[day_init_step] = longwave_sky.sky_temperature
                update_soil_properties!(output, soil_prop_view, soil_properties_model;
                    soil_temperature=T0, soil_moisture, bulk_density, mineral_density,
                    mineral_conductivity, mineral_heat_capacity,
                    atmospheric_pressure=output.pressure[day_init_step],
                    step=day_init_step, vapour_pressure_equation,
                )
                # PR #102: snow is now recorded at `step` (the hour-loop index)
                # not at `day_init_step`. Matches Fortran OSUB convention where
                # snow is computed and output at the same TIME point it is
                # updated. Skip pre-filling snow here — the per-hour loop writes it.
            end

            # ── Set up first hour's inputs and initialize integrator for full day ──
            step = (j - 1) * length(hours) + 1
            environment_instant = get_instant(environment_day, environment_hourly, output, soil_moisture, step)
            T0 = setindex(T0, environment_instant.deep_soil_temperature, num_soil_nodes)
            (; albedo, effective_wetness, longwave_sky) = apply_snow_overrides(
                snow_model, snow_state, buffers.snow, step, site, moisture_mode, environment_instant,
                longwave_model, vapour_pressure_equation;
                is_last_iter)
            if n_snow > 0
                # PR #102: on non-final iters (spinup) Fortran still maintains
                # `snowhr(methour)` and `snode` via the SNOWLAYER call inside
                # OSUB before the `IFINAL < ND` return — keep the depth array
                # populated so the ODE sees the same node geometry every iter.
                if !is_last_iter
                    buffers.snow.snow_depth_hourly[step] = snow_state.current_depth
                    snow_state = activate_snow_nodes!(snow_model, snow_state, buffers.snow, T_snow, step)
                end
                compute_effective_depths!(snow_model, snow_state, buffers.snow, depths)
                ode_depths = buffers.snow.effective_depths
                first_active = findfirst(d -> d > 0u"cm", buffers.snow.node_depths)
                sync_temp = isnothing(first_active) ? T0[1] : T_snow[first_active]
                T_snow = _sync_snow_to_active(T_snow, first_active, sync_temp)
            else
                ode_depths = depths
            end
            T0_ode = combine_ode_state(T_snow, T0, n_snow)
            update_soil_energy_inputs!(ode_integrator.p;
                depths=ode_depths, environment_instant,
                soil_wetness=effective_wetness, longwave_sky, albedo,
                Q_freeze, snow_state,
            )

            # Reset the integrator for a fresh day-pass (0-1440 min), matching
            # Fortran SFODE. The default `reinit!` allocates ~528 B per call:
            # ~440 B from `initialize!` rebuilding the Tsit5 stage cache
            # `Memory{SVector}` and ~88 B from the `tstops`/`saveat` push paths.
            # We do the minimal subset of `reinit!`'s work in place — the Tsit5
            # `ConstantCache` stage values are recomputed every step anyway, so
            # skipping the cache reinit is safe.
            _maybe_reinit_integrator!(ode_integrator, T0_ode)

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
                    T0 = apply_phase_transition(snow_model, soil_freezing_model, T0, T0_before, buffers.phase_transition, ∑phase, soil_moisture, depths)
                    T0_output = T0
                end

                # ── Snow post-ODE: mass balance, clamping, node activation ──
                # Only on final iteration (Fortran OSUB returns before snow code on non-final iterations)
                snowfall_hour = 0.0u"cm"
                thermal_melt = 0.0u"cm"
                if n_snow > 0 && is_last_iter
                    # Evaluate evaporation at end-of-hour state (Fortran OSUB.f line 744)
                    Q_evap, _ = surface_convection_evaporation(evaporation_model;
                        boundary_layer_model,
                        surface_temperature=T0[1],
                        air_temperature=environment_instant.reference_temperature,
                        wind_speed=environment_instant.reference_wind_speed,
                        relative_humidity=environment_instant.reference_humidity,
                        atmospheric_pressure=environment_instant.atmospheric_pressure,
                        roughness_height=site.roughness_height,
                        reference_height=last(heights),
                        soil_wetness=effective_wetness,
                        vapour_pressure_equation,
                    )
                    snow_state, snowfall_hour, thermal_melt, rain_melt_snow = update_snow!(snow_model, snow_state, buffers.snow, T_snow, T_snow_before, environment_instant, step,
                        T0[1], T0_before[1], Q_evap;
                        rainfall_schedule, shade=environment_instant.shade, day_of_year=days[j])
                    # Fortran OSUB.f lines 938-943: QFREZE = melted*(1-snowmelt)*79.7/60/10000
                    # Only applied when DOY > 1 (Fortran guard). Fortran units: cal/(min·cm²).
                    # melt depth (cm) × water density (kg/m³ → kg/m²/cm) × Lf (J/kg) / 1 hr → W/m²
                    melt_factor = snow_model.snow_melt_factor
                    if days[j] > 1 && melt_factor <= 1.0
                        Q_freeze = uconvert(u"W/m^2", thermal_melt * water_density * LATENT_HEAT_FUSION / 1u"hr" * (1.0 - melt_factor))
                    else
                        Q_freeze = 0.0u"W/m^2"
                    end
                    # Fortran SNOWLAYER.f lines 67-71, 93-97: reset deactivated node temps to T(1)
                    prev_active = snow_state.active_nodes
                    snow_state = activate_snow_nodes!(snow_model, snow_state, buffers.snow, T_snow, step)
                    if snow_state.active_nodes < prev_active
                        T_snow = sync_inactive_snow_temps(T_snow, buffers.snow)
                    end
                    snow_state, T_snow, soil_surf_temp = snow_phase_transition(snow_model, snow_state, buffers.snow, T_snow, T_snow_before, environment_instant.atmospheric_pressure, T0[1])
                    # NicheMapR-parity: skip the soil-surface clamp cascade
                    # from snow_phase_transition / clamp_snow_temperatures /
                    # freeze_new_snow. Keeping T0[1] at its post-ODE value lets
                    # the next-hour snow_thermal_melt fire at the right hours
                    # and matches NicheMapR's monthly-snow output (mean depth
                    # diff 0.48 cm over 12 days vs 2.15 cm with cascade on).
                    T_snow, _ = clamp_snow_temperatures(T_snow, buffers.snow)
                    # NicheMapR-parity: see comment above on the soil-surface
                    # clamp cascade.
                    T_snow, _ = freeze_new_snow(snow_model, T_snow, buffers.snow, step)
                    # NicheMapR-parity: see comment above on the soil-surface
                    # clamp cascade.
                    snow_state = adjust_snow_near_nodes!(snow_model, snow_state, buffers.snow, T_snow, step)
                    prev_active2 = snow_state.active_nodes
                    snow_state = activate_snow_nodes!(snow_model, snow_state, buffers.snow, T_snow, step)
                    # Fortran SNOWLAYER.f: also reset temps after the second activation
                    if snow_state.active_nodes < prev_active2
                        T_snow = sync_inactive_snow_temps(T_snow, buffers.snow)
                    end
                end

                # ── Feed modified state back into integrator for next hour ──
                # Matches Fortran OSUB.f line 1087: t=tt (copies modified temps back to ODE state)
                if i < length(hours)
                    next_step = step + 1
                    environment_instant = get_instant(environment_day, environment_hourly, output, soil_moisture, next_step)
                    T0 = setindex(T0, environment_instant.deep_soil_temperature, num_soil_nodes)
                    (; albedo, effective_wetness, longwave_sky) = apply_snow_overrides(
                        snow_model, snow_state, buffers.snow, next_step, site, moisture_mode, environment_instant,
                        longwave_model, vapour_pressure_equation;
                        is_last_iter)
                    if n_snow > 0
                        if !is_last_iter
                            buffers.snow.snow_depth_hourly[next_step] = snow_state.current_depth
                            snow_state = activate_snow_nodes!(snow_model, snow_state, buffers.snow, T_snow, next_step)
                        end
                        compute_effective_depths!(snow_model, snow_state, buffers.snow, depths)
                        ode_depths = buffers.snow.effective_depths
                    else
                        ode_depths = depths
                    end
                    update_soil_energy_inputs!(ode_integrator.p;
                        depths=ode_depths, environment_instant,
                        soil_wetness=effective_wetness, longwave_sky, albedo,
                        Q_freeze, snow_state,
                    )
                end
                T0_ode = combine_ode_state(T_snow, T0, n_snow)
                # Direct field write bypasses any custom setproperty! interposing
                # an SVector copy. The integrator's `.u` field is concretely
                # typed `SVector{N,T}` so this is a plain memory write.
                setfield!(ode_integrator, :u, T0_ode)
                SciMLBase.u_modified!(ode_integrator, true)

                init_surface_fluxes!(boundary_layer_model, buffers, forcing, site, heights, T0, i)
                # Fortran OSUB.f: soil moisture runs only on final iteration (after line 353 guard)
                if is_last_iter
                    rain = current_rainfall(rainfall_schedule;
                        environment_instant, environment_hourly=environment_hourly,
                        step, i, midnight_i,
                    )
                    # Fortran OSUB.f: apply rainmult to rainfall entering condep
                    if n_snow > 0
                        rain = rain * snow_model.rain_multiplier
                    end
                    # Pool update:
                    # Cold snow (air_temp < snow_threshold): rain → snowfall (handled in update_snow!); only thermal melt enters pool.
                    # Warm snow or no snow: rain passes through pack to soil surface, rain_melt snow water also drains to surface.
                    #   Fortran OSUB.f line 1116 zeros rainfallb when cursnow > 0 for non-hourly rain (mass conservation
                    #   violation — rain water discarded). Julia instead routes rain to pool, and also tracks rain_melt water
                    #   (snow melted by rain energy; OSUB.f line 950 subtracts rainmelt from cursnow but never adds it to condep
                    #   — another mass loss). Both rain and rain_melt water physically reach the soil surface.
                    snow_present = n_snow > 0 && buffers.snow.snow_depth_hourly[step] > 0.0u"cm"
                    snow_temp_threshold = n_snow > 0 ? snow_model.snow_temperature_threshold : 0.0u"°C"
                    melted_water = n_snow > 0 ? uconvert(u"kg/m^2", thermal_melt * 1.0u"g/cm^3") : 0.0u"kg/m^2"
                    rain_melt_water = n_snow > 0 ? uconvert(u"kg/m^2", rain_melt_snow * snow_state.density) : 0.0u"kg/m^2"
                    if snow_present && u"°C"(environment_instant.reference_temperature) < snow_temp_threshold
                        # Cold snow: rain absorbed into snowpack (handled in update_snow! as snowfall); only thermal melt enters pool
                        pool = clamp(pool + melted_water, 0.0u"kg/m^2", max_surface_pool)
                    else
                        # Warm snow or no snow: rain + rain_melt water + thermal melt enter pool
                        pool = clamp(pool + rain + rain_melt_water + melted_water, 0.0u"kg/m^2", max_surface_pool)
                    end
                    # Soil moisture physics; output write happens below at output_step
                    (; pool, soil_moisture, infil_out) = step_soil_moisture!(moisture_mode, buffers, soil_hydraulic_model;
                        soil_profile, depths, site, boundary_layer_model, environment_instant, T0, pool, soil_moisture,
                        max_surface_pool, evaporation_model, vapour_pressure_equation, snow_present,
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
                    update_soil_water!(output, soil_moisture, infil_out, output_step)
                    environment_instant = get_instant(environment_day, environment_hourly, output, soil_moisture, output_step)
                    update_soil_properties!(output, soil_prop_view, soil_properties_model;
                        soil_temperature=T0, soil_moisture, bulk_density, mineral_density,
                        mineral_conductivity, mineral_heat_capacity,
                        atmospheric_pressure=environment_instant.atmospheric_pressure, step=output_step, vapour_pressure_equation
                    )
                end
                # PR #102: snow output writes at `step` (not `step+1`), matching
                # Fortran OSUB convention where snow is computed and recorded at
                # the same TIME point it is updated. Row 1 of each day is the
                # post-midnight-snowfall state.
                if is_last_iter && n_snow > 0 && step <= ndays * length(hours)
                    output.snow_fall[step]    = uconvert(u"cm/hr", snowfall_hour / 1.0u"hr")
                    output.snow_depth[step]   = uconvert(u"cm", buffers.snow.snow_depth_hourly[step])
                    # Corrected NicheMapR convention: stored density at row N uses
                    # the post-update depth and snow_age at hour N. Internally
                    # `snow_state.density` is the pre-update density (used for
                    # snowfall conversion and as `prevden` in the next hour's
                    # compaction ratio). Recompute for output only.
                    output.snow_density[step] = uconvert(u"g/cm^3",
                        evolve_snow_density(snow_model, snow_state))
                    _write_row!(output.snow_temperature, step, T_snow)
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
    site = mp.inputs.site
    (; boundary_layer_model, vapour_pressure_equation, snow_model) = mp.model
    profile_buffers = cache.buffers.air_profile
    n_snow = n_snow_nodes(snow_model)
    for i in 1:size(output.profile.air_temperature, 1)
        # PR #102: when snow is present, the boundary-layer surface is the
        # top of the snowpack (snow node 1) not the soil surface.
        has_snow = n_snow > 0 && output.snow_depth[i] > 0.0u"cm"
        surface_temperature = has_snow ?
            u"°C"(output.snow_temperature[i, 1]) :
            u"°C"(output.soil_temperature[i, 1])
        environment_instant = (;
            atmospheric_pressure=output.pressure[i],
            reference_temperature=output.reference_temperature[i],
            reference_wind_speed=output.reference_wind_speed[i],
            reference_humidity=output.reference_humidity[i],
            zenith_angle=solar_radiation_out.zenith_angle[i],
        )
        result = atmospheric_surface_profile!(boundary_layer_model, profile_buffers; site, environment_instant, surface_temperature, vapour_pressure_equation)
        _write_row!(output.profile.air_temperature,   i, result.air_temperature)
        _write_row!(output.profile.wind_speed,        i, result.wind_speed)
        _write_row!(output.profile.relative_humidity, i, result.relative_humidity)
        output.profile.convective_heat_flux[i]  = result.convective_heat_flux
        output.profile.friction_velocity[i]     = result.friction_velocity
    end
end

function allocate_ode_integrator(model::SoilHeatTransport1D, T0, inputs_proto)
    # Full-day integration (0-1440 min), matching Fortran SFODE which integrates
    # continuously. Hourly bookkeeping is done via step!/add_tstop! in the day loop.
    tspan = (0.0u"minute", 1440.0u"minute")
    prob = ODEProblem{false}(soil_energy_balance, T0, tspan, inputs_proto)
    return SciMLBase.init(prob, model.ode_solver;
        save_everystep=false, save_start=false, save_end=false, model.ode_kwargs...)
end

# Reset integrator state for a fresh pass through the (0..1440 min) day.
# Multi-step solvers (Adams family) need their f-history rebuilt via
# `reinit!`'s cache pass; single-step Runge–Kutta solvers (Tsit5, RK4, …)
# recompute their stage cache every step anyway, so we skip `reinit_cache`
# and avoid the ~440 B/call allocation. `integrator.alg` has a concrete
# type at the call site, so the compiler constant-folds the branch.
@inline function _maybe_reinit_integrator!(integrator, u_new)
    if ismultistep(integrator.alg)
        SciMLBase.reinit!(integrator, u_new;
            t0=zero(integrator.t), erase_sol=false, reset_dt=true,
        )
    else
        SciMLBase.reinit!(integrator, u_new;
            t0=zero(integrator.t), erase_sol=false, reset_dt=true,
            reinit_cache=false,
        )
    end
    return integrator
end

# Step the integrator to an exact hour boundary and return the state.
# Fortran DSUB.f:551 clamps T(I) in-place at every derivative call. Our
# `soil_energy_balance` only clamps the local copy used for derivatives, so
# fixed-step solvers (e.g. ABM43) can let the integrator state drift outside
# physical bounds. We re-clamp the state after every internal step here so
# any solver — adaptive or fixed-step — stays in [-81°C, max_surface_T].
function step_to_hour!(integrator, hour_index)
    target_t = hour_index * 60.0u"minute"
    SciMLBase.add_tstop!(integrator, target_t)
    T_lo = u"K"(-81.0u"°C")
    T_hi = u"K"(integrator.p.model.maximum_surface_temperature)
    while integrator.t < target_t - 1e-10u"minute"
        SciMLBase.step!(integrator)
        u_clamped = map(t -> clamp(t, T_lo, T_hi), integrator.u)
        if u_clamped !== integrator.u
            setfield!(integrator, :u, u_clamped)
            setfield!(integrator, :uprev, u_clamped)
            SciMLBase.u_modified!(integrator, true)
        end
    end
    return integrator.u
end

function update_soil_properties!(output, soil_properties_buffers, soil_properties_model; step, kw...)
    (; bulk_thermal_conductivity, bulk_heat_capacity, bulk_density) = soil_properties!(soil_properties_buffers, soil_properties_model; kw...)

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

# Explicit-loop column read from a Matrix into a Vector. Avoids the
# `dst .= arr[:, col]` SubArray + Broadcasted heap allocations.
@inline function _read_col!(dst::AbstractVector, arr::AbstractMatrix, col::Int)
    @inbounds for i in eachindex(dst)
        dst[i] = arr[i, col]
    end
    return dst
end

function update_soil_water!(output, soil_moisture, infil_out, step)
    _write_row!(output.soil_moisture, step, soil_moisture)
    if !isnothing(infil_out)
        _write_row!(output.soil_water_potential, step, infil_out.soil_water_potential)
        _write_row!(output.soil_humidity, step, infil_out.soil_humidity)
    end
    return output
end

# ── Soil moisture stepping dispatch ───────────────────────────────────────

function step_soil_moisture!(mode::DynamicSoilMoisture, buffers, soil_hydraulic_model;
    soil_profile, depths, site, boundary_layer_model, environment_instant, T0, pool, soil_moisture,
    max_surface_pool, evaporation_model, vapour_pressure_equation, snow_present=false,
)
    (; moisture_tolerance, moisture_max_iterations, moisture_timestep) = mode
    niter_moist = ustrip(u"s^-1", 3600 / moisture_timestep)
    (; infil_out, soil_wetness, pool, soil_moisture) = soil_water_balance!(buffers, soil_hydraulic_model;
        soil_profile, depths, site, boundary_layer_model, environment_instant, T0, niter_moist, pool,
        soil_wetness=mode.soil_wetness, soil_moisture,
        moisture_timestep, moisture_tolerance, moisture_max_iterations, max_surface_pool,
        evaporation_model, vapour_pressure_equation, snow_present,
    )
    mode.soil_wetness = soil_wetness
    return (; pool, soil_moisture, infil_out)
end

step_soil_moisture!(::PrescribedSoilMoisture, args...; pool, soil_moisture, kw...) = (; pool, soil_moisture, infil_out=nothing)

"""
Build a single `Forcing` whose 8 `ScaledInterpolation`s wrap freshly-allocated
25-element coefficient buffers. Built ONCE in `init` and reused across every
day-iter via `update_forcing_day!` which only mutates the underlying `itp.coefs`.

Fortran MICROCLIMATE.f lines 843-852: 25 interpolation nodes (0,60,...,1440).
The 25th node repeats the 24th hour's value (TD(36)=TAIRhr1(DOYF)).
"""
function allocate_forcing(output, solar_radiation_out)
    tspan = 0.0:60:(60*24)  # 0, 60, ..., 1440
    _scaled(buf) = scale(interpolate(buf, BSpline(Interpolations.Linear())), tspan)
    return Forcing(
        _scaled(zeros(eltype(output.global_radiation), 25)),                  # solar
        _scaled(zeros(eltype(solar_radiation_out.zenith_angle), 25)),         # zenith
        _scaled(zeros(eltype(solar_radiation_out.zenith_slope_angle), 25)),   # slope_zenith
        _scaled(zeros(typeof(1.0u"K"), 25)),                                  # temperature
        _scaled(zeros(eltype(output.reference_wind_speed), 25)),              # wind
        _scaled(zeros(eltype(output.reference_humidity), 25)),                # humidity
        _scaled(zeros(eltype(output.cloud_cover), 25)),                       # cloud
        _scaled(zeros(eltype(output.pressure), 25)),                          # pressure
        _scaled(zeros(eltype(output.diffuse_fraction), 25)),                  # diffuse_fraction
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
Update the in-place `Forcing` for day `iday` by writing into each
`itp.coefs` buffer. No allocation.
"""
function update_forcing_day!(forcing::Forcing, solar_radiation_out, output, iday::Int)
    nhours = 24
    sub1 = (iday*nhours - nhours + 1):(iday*nhours)
    # `forcing.<field>.itp.coefs` is the 25-element buffer underneath the
    # ScaledInterpolation, mutated in place.
    _copy_day_into!(forcing.solar.itp.coefs,            output.global_radiation,                sub1)
    _copy_day_into!(forcing.zenith.itp.coefs,           solar_radiation_out.zenith_angle,       sub1)
    _copy_day_into!(forcing.slope_zenith.itp.coefs,     solar_radiation_out.zenith_slope_angle, sub1)
    _copy_day_into!(forcing.temperature.itp.coefs,      output.reference_temperature,           sub1, uconvert, u"K")
    _copy_day_into!(forcing.wind.itp.coefs,             output.reference_wind_speed,            sub1)
    _copy_day_into!(forcing.humidity.itp.coefs,         output.reference_humidity,              sub1)
    _copy_day_into!(forcing.cloud.itp.coefs,            output.cloud_cover,                     sub1)
    _copy_day_into!(forcing.pressure.itp.coefs,         output.pressure,                        sub1)
    _copy_day_into!(forcing.diffuse_fraction.itp.coefs, output.diffuse_fraction,                sub1)
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
        # `nothing` here → `apply_snow_overrides` falls back to `site.albedo`.
        # A populated `environment_daily.albedo` overrides per-day.
        albedo = _maybe_indexed(environment_daily.albedo, iday),
    )
end

# Type-stable accessor for optional daily fields. Branches on whether the
# field is `nothing` so the inferred return type is concrete per call site.
@inline _maybe_indexed(::Nothing, _) = nothing
@inline _maybe_indexed(v, iday) = v[iday]
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

# Generic per-layer field selection used by `soil_properties!`. The
# variant-specific method (`maybegetindex(::CampbelldeVriesSoilProperties, i)`)
# lives in `soil_properties/campbell_devries.jl`.
maybegetindex(props::NamedTuple, i::Int) = map(p -> maybegetindex(p, i), props)
maybegetindex(val::Number, i::Int) = val
maybegetindex(vals::AbstractArray, i::Int) = vals[i]
