# ── Snow model types ──────────────────────────────────────────────────────

struct NoSnow end

struct SnowModel{N}
    snow_temperature_threshold::typeof(1.0u"°C")
    snow_density::typeof(1.0u"g/cm^3")
    snow_melt_factor::Float64
    undercatch::Float64
    rain_multiplier::Float64
    rain_melt_factor::Float64
    density_function::NTuple{4,Float64}  # empirical coefficients (dimensionless)
    snow_conductivity::typeof(1.0u"W/m/K")  # 0 W/m/K = use Aggarwal
    canopy_interception::Float64
    grass_shade::Bool
    min_snow_depth::typeof(1.0u"cm")
    snow_node_thresholds::NTuple{N,typeof(1.0u"cm")}
    melt_threshold::typeof(1.0u"°C")
end

function SnowModel(;
    snow_temperature_threshold=1.5u"°C",
    snow_density=0.375u"g/cm^3",
    snow_melt_factor=1.0,
    undercatch=1.0,
    rain_multiplier=1.0,
    rain_melt_factor=0.0125,
    density_function=(0.0, 0.0, 0.0, 0.0),
    snow_conductivity=0.0u"W/m/K",
    canopy_interception=0.0,
    grass_shade=false,
    min_snow_depth=2.0u"cm",
    snow_node_thresholds=DEFAULT_SNOW_NODE_THRESHOLDS .* u"cm",
    melt_threshold=0.4u"°C",
)
    N = length(snow_node_thresholds)
    SnowModel{N}(
        snow_temperature_threshold, snow_density, snow_melt_factor,
        undercatch, rain_multiplier, rain_melt_factor, density_function,
        snow_conductivity, canopy_interception, grass_shade, min_snow_depth,
        NTuple{N,typeof(1.0u"cm")}(snow_node_thresholds), melt_threshold,
    )
end

n_snow_nodes(::NoSnow) = 0
n_snow_nodes(::SnowModel{N}) where N = N

# ── Immutable snow state ─────────────────────────────────────────────────

struct SnowState
    current_depth::typeof(1.0u"cm")
    snow_age::Float64
    days_since_snow::Float64
    density::typeof(1.0u"g/cm^3")
    previous_density::typeof(1.0u"g/cm^3")
    cumulative_melt::typeof(1.0u"cm")
    extra_snow::typeof(1.0u"kg/m^2")
    active_nodes::Int
    sum_phase::typeof(1.0u"J/m^2")
end

initial_snow_state(::NoSnow) = nothing
function initial_snow_state(sm::SnowModel)
    SnowState(0.0u"cm", 0.0, 0.3, sm.snow_density, 1.0u"g/cm^3",
              0.0u"cm", 0.0u"kg/m^2", 0, 0.0u"J/m^2")
end

# ── Pre-allocated scratch arrays ─────────────────────────────────────────

allocate_snow_scratch(::NoSnow, args...) = nothing

function allocate_snow_scratch(::SnowModel{N}, nsteps, num_ode_nodes, depths) where N
    return (;
        node_depths           = fill(0.0u"cm", N),
        snow_depth_hourly     = fill(0.0u"cm", nsteps),
        mean_temperature      = fill(0.0u"°C", N),
        mean_temperature_past = fill(0.0u"°C", N),
        phase_heat            = fill(0.0u"J/m^2", N),
        layer_mass            = fill(0.0u"kg/m^2", N),
        effective_depths      = fill(zero(eltype(depths)), num_ode_nodes),
    )
end

# ── Scratch reset (for independent days) ─────────────────────────────────

reset_snow_scratch!(::NoSnow, _) = nothing
function reset_snow_scratch!(::SnowModel, scratch)
    scratch.node_depths .= 0.0u"cm"
    scratch.snow_depth_hourly .= 0.0u"cm"
    scratch.phase_heat .= 0.0u"J/m^2"
    scratch.layer_mass .= 0.0u"kg/m^2"
    return nothing
end

# ── Snow node activation ─────────────────────────────────────────────────

activate_snow_nodes(::NoSnow, args...) = nothing
function activate_snow_nodes(snow_model::SnowModel{N}, state::SnowState, scratch, snow_temperature, step) where N
    (; node_depths, snow_depth_hourly) = scratch
    cursnow = snow_depth_hourly[step]
    active_nodes = state.active_nodes

    if cursnow > 300.0u"cm"
        active_nodes = 0
    end

    if cursnow < snow_model.min_snow_depth
        node_depths .= 0.0u"cm"
        snow_depth_hourly[step] = 0.0u"cm"
        return setproperties(state, (; active_nodes=0, current_depth=0.0u"cm"))
    end

    maxsnode = 0
    for i in 1:N
        if cursnow > snow_model.snow_node_thresholds[i]
            maxsnode = i
        end
    end
    maxsnode = min(maxsnode, N - 1)

    for i in 1:maxsnode
        node_depths[i + (N - maxsnode)] = snow_model.snow_node_thresholds[i]
    end

    if maxsnode < active_nodes
        for i in 1:(N - maxsnode)
            node_depths[i] = 0.0u"cm"
        end
    end

    return setproperties(state, (; current_depth=cursnow, active_nodes=maxsnode))
end

# ── Effective depth vectors for the ODE ──────────────────────────────────

function combined_depths(snow_model::SnowModel{N}, state::SnowState, scratch, soil_depths) where N
    cursnow = state.current_depth
    snow_offset = cursnow >= snow_model.min_snow_depth ? cursnow : 0.0u"cm"
    return vcat(
        SVector(ntuple(i -> uconvert(unit(eltype(soil_depths)), scratch.node_depths[i]), Val(N))),
        SVector(ntuple(i -> soil_depths[i] + uconvert(unit(eltype(soil_depths)), snow_offset), Val(length(soil_depths))))
    )
end

compute_effective_depths!(::NoSnow, _, _, _) = nothing
function compute_effective_depths!(snow_model::SnowModel{N}, state::SnowState, scratch, soil_depths) where N
    (; effective_depths, node_depths) = scratch
    cursnow = state.current_depth
    snow_offset = cursnow >= snow_model.min_snow_depth ? cursnow : 0.0u"cm"

    for i in 1:N
        effective_depths[i] = uconvert(unit(eltype(soil_depths)), node_depths[i])
    end
    for i in eachindex(soil_depths)
        effective_depths[N + i] = soil_depths[i] + uconvert(unit(eltype(soil_depths)), snow_offset)
    end
    return nothing
end

# ── Snow thermal properties ──────────────────────────────────────────────

function write_snow_properties!(snow_model::SnowModel{N}, state::SnowState, scratch, snow_temperature,
    property_buffers, atmospheric_pressure, vapour_pressure_equation=GoffGratch(),
) where N
    snowdens = evolve_snow_density(snow_model, state)
    state = setproperties(state, (; density=snowdens))

    (; node_depths) = scratch

    if state.current_depth >= snow_model.min_snow_depth
        # Mixing ratio at 0°C
        e_sat = wet_air_properties(273.15u"K", 1.0, atmospheric_pressure; vapour_pressure_equation).vapour_pressure
        mixing_ratio = 0.622 * e_sat / (atmospheric_pressure - e_sat)
        # Specific heat: ice + humid air weighted by density fraction
        density_fraction = ustrip(u"g/cm^3", snowdens)
        specific_heat = (2100.0 * density_fraction + (1005.0 + 1820.0 * mixing_ratio) * (1.0 - density_fraction))u"J/kg/K"

        # Conductivity (Aggarwal 2009 or user-defined)
        conductivity = if snow_model.snow_conductivity > 0.0u"W/m/K"
            snow_model.snow_conductivity
        else
            ρ = ustrip(u"kg/m^3", snowdens)
            (0.00395 + 0.00084 * ρ - 1.7756e-6 * ρ^2 + 3.80635e-9 * ρ^3)u"W/m/K"
        end

        for i in 1:N
            property_buffers.bulk_thermal_conductivity[i] = uconvert(unit(property_buffers.bulk_thermal_conductivity[N+1]), conductivity)
            temp = u"°C"(snow_temperature[i])
            if temp > -0.45u"°C" && temp <= 0.4u"°C"
                property_buffers.bulk_heat_capacity[i] = uconvert(unit(property_buffers.bulk_heat_capacity[N+1]), specific_heat + uconvert(u"J/kg/K", LATENT_HEAT_FUSION / 1.0u"K"))
            else
                property_buffers.bulk_heat_capacity[i] = uconvert(unit(property_buffers.bulk_heat_capacity[N+1]), specific_heat)
            end
            is_active = node_depths[i] > 0.0u"cm" || (i < N && node_depths[min(N, i + 1)] > 0.0u"cm")
            property_buffers.bulk_density[i] = uconvert(unit(property_buffers.bulk_density[N+1]), is_active ? snowdens : 0.0u"g/cm^3")
        end
    else
        for i in 1:N
            property_buffers.bulk_thermal_conductivity[i] = property_buffers.bulk_thermal_conductivity[N + 1]
            property_buffers.bulk_heat_capacity[i] = property_buffers.bulk_heat_capacity[N + 1]
            property_buffers.bulk_density[i] = uconvert(unit(property_buffers.bulk_density[N+1]), 0.0u"g/cm^3")
        end
    end
    return state
end

# ── Snow temperature operations ──────────────────────────────────────────

function clamp_snow_temperatures(snow_temperature::SVector{N}, scratch) where N
    freeze = u"K"(0.0u"°C")
    return SVector(ntuple(Val(N)) do i
        scratch.node_depths[i] > 0.0u"cm" && snow_temperature[i] > freeze ? freeze : snow_temperature[i]
    end)
end

function freeze_new_snow(snow_temperature::SVector{N}, scratch, step) where N
    if step > 1 && scratch.snow_depth_hourly[step] > 0.0u"cm" && scratch.snow_depth_hourly[step - 1] < 1e-8u"cm"
        freeze = u"K"(0.0u"°C")
        return SVector(ntuple(Val(N)) do i
            snow_temperature[i] > freeze ? freeze : snow_temperature[i]
        end)
    end
    return snow_temperature
end

# ── Snow albedo ──────────────────────────────────────────────────────────

function snow_albedo(days_since_snow)
    (-9.874 * log(max(days_since_snow, 0.3)) + 78.3434) / 100.0
end

# ── Snow density evolution ───────────────────────────────────────────────
# The density function coefficients are empirical and operate on unitless
# depth (cm) and age values, producing a density fraction (g/cm³ numerically).
# Returns the (possibly updated) density value.

function evolve_snow_density(snow_model::SnowModel, state::SnowState)
    densfun = snow_model.density_function
    densfun[1] > 0 || return state.density
    cursnow = ustrip(u"cm", state.current_depth)
    snowage = state.snow_age
    new_density = if densfun[3] > 0
        (densfun[1] - densfun[2]) * (1.0 - exp(-densfun[3] * cursnow - densfun[4] * snowage)) + densfun[2]
    else
        min(0.9167, densfun[1] * snowage + densfun[2])
    end
    return new_density * u"g/cm^3"
end

# ── Snow surface overrides ───────────────────────────────────────────────

snow_surface_overrides(::NoSnow, _, _, _) = (;
    albedo=nothing, emissivity=nothing, soil_wetness=nothing, shade=nothing
)
function snow_surface_overrides(snow_model::SnowModel, state::SnowState, scratch, step)
    if step > 1 && scratch.snow_depth_hourly[step - 1] > 0.0u"cm"
        return (;
            albedo    = snow_albedo(state.days_since_snow),
            emissivity = nothing,
            soil_wetness = 1.0,
            shade = snow_model.grass_shade ? 0.0 : nothing,
        )
    end
    return (; albedo=nothing, emissivity=nothing, soil_wetness=nothing, shade=nothing)
end

# ── Snow layer geometry ──────────────────────────────────────────────────

function snow_layer_geometry(snow_model::SnowModel{N}, state::SnowState, scratch, layer_index) where N
    (; node_depths) = scratch
    maxsnode = state.active_nodes
    node_index = clamp(layer_index + N - 1 - maxsnode, 1, N)
    thickness = if layer_index != maxsnode + 1
        upper = max(1, layer_index + N - maxsnode)
        node_depths[upper] - node_depths[node_index]
    else
        max(0.0u"cm", state.current_depth - snow_model.snow_node_thresholds[max(1, maxsnode)])
    end
    return (; node_index, thickness)
end

# ── Snow specific heat ───────────────────────────────────────────────────
# Empirical formula: weighted average of ice and humid air specific heats.

function snow_specific_heat(snow_density, atmospheric_pressure)
    e_sat = wet_air_properties(273.15u"K", 1.0, atmospheric_pressure).vapour_pressure
    mixing_ratio = 0.622 * e_sat / (atmospheric_pressure - e_sat)
    density_fraction = ustrip(u"g/cm^3", snow_density)
    return (2100.0 * density_fraction + (1005.0 + 1820.0 * mixing_ratio) * (1.0 - density_fraction))u"J/kg/K"
end

# ── Thermal melt ─────────────────────────────────────────────────────────

function snow_thermal_melt(snow_model::SnowModel{N}, state::SnowState, scratch,
    snow_temperature, snow_temperature_before, atmospheric_pressure,
) where N
    maxsnode = state.active_nodes
    snowdens = state.density
    (; mean_temperature, mean_temperature_past) = scratch

    for i in 1:N
        mean_temperature[i] = u"°C"((snow_temperature[i] + snow_temperature[min(N, i + 1)]) / 2)
        mean_temperature_past[i] = u"°C"((snow_temperature_before[i] + snow_temperature_before[min(N, i + 1)]) / 2)
    end

    specific_heat = snow_specific_heat(snowdens, atmospheric_pressure)

    total_melt = 0.0u"kg/m^2"
    for j in 1:(maxsnode + 1)
        (; node_index, thickness) = snow_layer_geometry(snow_model, state, scratch, j)
        layer_mass = uconvert(u"kg/m^2", thickness * snowdens)

        if mean_temperature[node_index] > snow_model.melt_threshold
            ΔT = mean_temperature[node_index] - snow_model.melt_threshold
            sensible_energy = uconvert(u"J/m^2", ΔT * specific_heat * layer_mass)
            latent_energy = uconvert(u"J/m^2", LATENT_HEAT_FUSION * layer_mass)
            total_melt += (sensible_energy + latent_energy) / LATENT_HEAT_FUSION
        end
    end

    return uconvert(u"cm", total_melt / snowdens)
end

# ── Snow phase transition (freezing within snowpack) ─────────────────────

snow_phase_transition(::NoSnow, args...) = nothing
function snow_phase_transition(snow_model::SnowModel{N}, state::SnowState, scratch,
    snow_temperature, snow_temperature_before, atmospheric_pressure,
) where N
    maxsnode = state.active_nodes
    (; mean_temperature, mean_temperature_past, phase_heat, layer_mass) = scratch

    state.current_depth < 200.0u"cm" || return state

    specific_heat = snow_specific_heat(state.density, atmospheric_pressure)

    phase_heat .= 0.0u"J/m^2"
    layer_mass .= 0.0u"kg/m^2"
    for j in 1:(maxsnode + 1)
        j < N && j != maxsnode + 1 || continue
        (; node_index, thickness) = snow_layer_geometry(snow_model, state, scratch, j)
        mass = uconvert(u"kg/m^2", thickness * state.density)

        if mean_temperature_past[node_index] >= 0.0u"°C" && mean_temperature[node_index] <= 0.0u"°C"
            layer_mass[node_index] = max(0.0u"kg/m^2", mass)
            phase_heat[node_index] = uconvert(u"J/m^2",
                (mean_temperature_past[node_index] - mean_temperature[node_index]) * specific_heat * layer_mass[node_index])
        end
    end

    new_sum_phase = state.sum_phase + sum(phase_heat)
    total_mass = sum(layer_mass)
    if total_mass > 0.0u"kg/m^2" && new_sum_phase > uconvert(u"J/m^2", LATENT_HEAT_FUSION * total_mass)
        new_sum_phase = 0.0u"J/m^2"
    end
    return setproperties(state, (; sum_phase=new_sum_phase))
end

# ── Snow mass balance ────────────────────────────────────────────────────

update_snow(::NoSnow, args...; kw...) = (nothing, 0.0u"cm")

function update_snow(snow_model::SnowModel{N}, state::SnowState, scratch,
    snow_temperature, snow_temperature_before, environment_instant, step;
    hourly_rainfall=false, shade=0.0,
) where N
    (; snow_depth_hourly) = scratch

    air_temperature = environment_instant.reference_temperature
    rainfall = environment_instant.rainfall
    is_first_step = (step == 1)
    is_midnight = (mod(step - 1, 24) == 0)

    # Density evolution
    prev_dens = state.density
    new_dens = evolve_snow_density(snow_model, state)
    density_ratio = new_dens / prev_dens
    if step > 1 && density_ratio != 0 && isfinite(ustrip(density_ratio))
        snow_depth_hourly[step - 1] /= ustrip(density_ratio)
    end

    # ── Snowfall ──
    snowfall = 0.0u"cm"
    days_since_snow = state.days_since_snow
    if u"°C"(air_temperature) <= snow_model.snow_temperature_threshold && rainfall > 0.0u"kg/m^2"
        if is_midnight || hourly_rainfall
            snowfall = uconvert(u"cm", rainfall * snow_model.rain_multiplier * snow_model.undercatch / new_dens)
            if shade > 0
                snowfall *= (1.0 - snow_model.canopy_interception)
            end
            days_since_snow = 0.3
        else
            days_since_snow += 1.0 / 25.0
        end
    else
        days_since_snow += 1.0 / 25.0
    end
    days_since_snow = max(days_since_snow, 0.3)

    # ── Sublimation ──
    sublimation = 0.0u"cm"
    snow_surface_temperature = snow_temperature[1]
    if snow_surface_temperature > u"K"(0.0u"°C")
        atmospheric_pressure = environment_instant.atmospheric_pressure
        relative_humidity = environment_instant.reference_humidity
        wet_air_out = wet_air_properties(u"K"(air_temperature), relative_humidity, atmospheric_pressure)
        heat_transfer_coefficient = 5.0u"W/m^2/K"
        mass_transfer_coefficient = (heat_transfer_coefficient / (wet_air_out.specific_heat * wet_air_out.density)) * (0.71 / 0.60)^0.666
        _, evaporation_mass_flux = evaporation(;
            surface_temperature=u"K"(snow_surface_temperature),
            air_temperature=u"K"(air_temperature),
            relative_humidity,
            surface_relative_humidity=1.0,
            mass_transfer_coefficient,
            atmospheric_pressure,
            soil_wetness=1.0,
            saturated=false,
        )
        evaporation_mass_flux = max(0.0u"g/s/m^2", evaporation_mass_flux)
        sublimation = uconvert(u"cm", evaporation_mass_flux * 1.0u"hr" / new_dens)
    end

    # ── Rain melt (Anderson model) ──
    rain_melt = 0.0u"cm"
    if u"°C"(air_temperature) >= snow_model.snow_temperature_threshold && rainfall > 0.0u"kg/m^2"
        if is_midnight || hourly_rainfall
            water_depth = uconvert(u"cm", rainfall / water_density)
            density_fraction = ustrip(u"g/cm^3", new_dens)
            air_temp_C = ustrip(u"°C", air_temperature)
            rain_melt = snow_model.rain_melt_factor * water_depth * air_temp_C * 24.0 / density_fraction
            if hourly_rainfall
                rain_melt /= 24.0
            end
            rain_melt = max(0.0u"cm", rain_melt)
        end
    end

    # ── Thermal melt ──
    # snow_thermal_melt needs state with updated density
    state_for_melt = setproperties(state, (; density=new_dens))
    thermal_melt = if is_first_step
        0.0u"cm"
    else
        snow_thermal_melt(snow_model, state_for_melt, scratch, snow_temperature, snow_temperature_before,
            environment_instant.atmospheric_pressure) * snow_model.snow_melt_factor
    end
    cumulative_melt = state.cumulative_melt + thermal_melt

    # ── Net accumulation ──
    carried_over = uconvert(u"cm", state.extra_snow / new_dens)
    net_accumulation = max(0.0u"cm", snowfall - sublimation) + carried_over
    extra_snow = 0.0u"kg/m^2"

    # ── Snow depth update ──
    previous_depth = step > 1 ? snow_depth_hourly[step - 1] : 0.0u"cm"
    new_depth = max(0.0u"cm", previous_depth + net_accumulation - rain_melt - thermal_melt)

    if new_depth < snow_model.min_snow_depth && (step <= 1 || snow_depth_hourly[step - 1] < 1e-8u"cm")
        extra_snow = uconvert(u"kg/m^2", new_depth * new_dens)
        new_depth = 0.0u"cm"
    end

    snow_depth_hourly[step] = new_depth

    # ── Snow age tracking ──
    snow_age = if is_first_step
        0.0
    elseif snow_depth_hourly[step] > 0.0u"cm" && snow_depth_hourly[step - 1] > 0.0u"cm"
        state.snow_age + 1.0 / 25.0
    else
        0.0
    end

    new_state = SnowState(new_depth, snow_age, days_since_snow, new_dens, new_dens,
                          cumulative_melt, extra_snow, state.active_nodes, state.sum_phase)
    return (new_state, snowfall)
end

# ── Adjust snow depth near node boundaries ────────────────────────────────

adjust_snow_near_nodes(::NoSnow, args...) = nothing
function adjust_snow_near_nodes(snow_model::SnowModel, state::SnowState, scratch, step)
    (; snow_depth_hourly) = scratch

    snowout = snow_depth_hourly[step]
    maxsnode = state.active_nodes
    threshold = snow_model.snow_node_thresholds[max(1, maxsnode)]
    extra_snow = state.extra_snow

    if step > 1 && snowout > 0.0u"cm"
        if (snowout - threshold < 0.5u"cm") && (snowout <= snow_depth_hourly[step - 1])
            extra_snow += uconvert(u"kg/m^2", (snowout - threshold) * state.density)
            snow_depth_hourly[step] = max(0.0u"cm", threshold - 0.1u"cm")
        end
    end

    snowout = snow_depth_hourly[step]
    if snowout > 0.0u"cm"
        threshold = snow_model.snow_node_thresholds[max(1, state.active_nodes)]
        if snowout < snow_model.min_snow_depth + 0.5u"cm" || snowout < threshold + 0.5u"cm"
            extra_snow = uconvert(u"kg/m^2", 0.5u"cm" * state.density)
            snow_depth_hourly[step] = threshold + 0.5u"cm"
        end
    end
    current_depth = snow_depth_hourly[step]
    return setproperties(state, (; extra_snow, current_depth))
end

# ── Phase transition dispatch ────────────────────────────────────────────

function apply_phase_transition(::NoSnow, soil_temperature, soil_temperature_past, buffers, accumulated_latent_heat, soil_moisture, depths)
    phase_transition!(buffers;
        temperatures=soil_temperature, temperatures_past=soil_temperature_past,
        accumulated_latent_heat, soil_moisture, depths,
    )
end

function apply_phase_transition(::SnowModel, soil_temperature, soil_temperature_past, buffers, accumulated_latent_heat, soil_moisture, depths)
    phase_transition!(buffers;
        temperatures=soil_temperature, temperatures_past=soil_temperature_past,
        accumulated_latent_heat, soil_moisture, depths,
    )
end
