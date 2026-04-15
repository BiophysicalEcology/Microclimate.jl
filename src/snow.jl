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
    # Evolve density for property calculations only — do NOT update state.density here.
    # The density update belongs in update_snow, which needs the old density for the
    # compaction ratio (densrat = snowdens/prevden, Fortran OSUB.f line 754).
    snowdens = evolve_snow_density(snow_model, state)

    (; node_depths) = scratch

    if state.current_depth >= snow_model.min_snow_depth
        specific_heat = snow_specific_heat(snowdens, atmospheric_pressure; vapour_pressure_equation)

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

"""
Clamp active snow node temperatures to ≤ 0°C, and if the last snow node (N) is active
and above 0°C, also clamp the first soil node.  Matches Fortran OSUB.f lines 976-989.
Returns `(snow_temperature, clamp_soil_surface::Bool)`.
"""
function clamp_snow_temperatures(snow_temperature::SVector{N}, scratch) where N
    freeze = u"K"(0.0u"°C")
    clamp_soil = scratch.node_depths[N] > 0.0u"cm" && snow_temperature[N] > freeze
    snow_temperature = SVector(ntuple(Val(N)) do i
        scratch.node_depths[i] > 0.0u"cm" && snow_temperature[i] > freeze ? freeze : snow_temperature[i]
    end)
    return (snow_temperature, clamp_soil)
end

"""
Freeze recently fallen snow nodes.  Fortran OSUB.f lines 990-1013: computes exactly
which nodes the new snow covers (`maxsnode3`), freezes only those.  If zero nodes
are covered, freezes the soil surface node instead.

Returns `(snow_temperature, freeze_soil_surface::Bool)`.
"""
function freeze_new_snow(snow_model::SnowModel{N}, snow_temperature::SVector{N}, scratch, step) where N
    if step > 1 && scratch.snow_depth_hourly[step] > 0.0u"cm" && scratch.snow_depth_hourly[step - 1] < 1e-8u"cm"
        freeze = u"K"(0.0u"°C")
        depth = scratch.snow_depth_hourly[step]
        # Count how many node thresholds the new snow exceeds
        maxsnode3 = 0
        for i in 1:N
            if depth > snow_model.snow_node_thresholds[i]
                maxsnode3 = i
            end
        end
        if maxsnode3 == 0
            # Snow is shallower than all thresholds — freeze soil surface (node 9 in Fortran)
            return (snow_temperature, true)
        else
            # Freeze only the specific snow nodes that the new snow covers
            # Fortran: do i=9-maxsnode3, 8-maxsnode3 → Julia indices: (N+1-maxsnode3):(N-maxsnode3)
            # But since N=8, this is (9-maxsnode3):(8-maxsnode3), matching Fortran's node mapping
            lo = N + 1 - maxsnode3
            hi = N - maxsnode3
            snow_temperature = SVector(ntuple(Val(N)) do i
                (i >= lo && i <= hi && snow_temperature[i] > freeze) ? freeze : snow_temperature[i]
            end)
            return (snow_temperature, false)
        end
    end
    return (snow_temperature, false)
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
            emissivity = 0.98,  # Fortran OSUB.f line 1028: SLE = 0.98
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

function snow_specific_heat(snow_density, atmospheric_pressure; vapour_pressure_equation=GoffGratch())
    wap = wet_air_properties(273.15u"K", 1.0, atmospheric_pressure; vapour_pressure_equation)
    mixing_ratio = wap.mixing_ratio
    density_fraction = ustrip(u"g/cm^3", snow_density)
    # Fortran OSUB.f line 847: RW/1.+RW = 2*RW due to operator precedence
    return (2100.0 * density_fraction + (1005.0 + 1820.0 * (mixing_ratio / 1.0 + mixing_ratio)) * (1.0 - density_fraction))u"J/kg/K"
end

# ── Thermal melt ─────────────────────────────────────────────────────────

function snow_thermal_melt(snow_model::SnowModel{N}, state::SnowState, scratch,
    snow_temperature, snow_temperature_before, atmospheric_pressure,
    soil_surface_temperature, soil_surface_temperature_before,
) where N
    maxsnode = state.active_nodes
    snowdens = state.density
    (; mean_temperature, mean_temperature_past) = scratch

    for i in 1:N
        t_below = i < N ? snow_temperature[i + 1] : soil_surface_temperature
        t_below_past = i < N ? snow_temperature_before[i + 1] : soil_surface_temperature_before
        mean_temperature[i] = u"°C"((snow_temperature[i] + t_below) / 2)
        mean_temperature_past[i] = u"°C"((snow_temperature_before[i] + t_below_past) / 2)
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

    # Fortran converts g/m² to cm via *0.0001, i.e. divides by ice density (1 g/cm³)
    return uconvert(u"cm", total_melt / 1.0u"g/cm^3")
end

# ── Snow phase transition (freezing within snowpack) ─────────────────────

snow_phase_transition(::NoSnow, args...) = (nothing, nothing, nothing)
function snow_phase_transition(snow_model::SnowModel{N}, state::SnowState, scratch,
    snow_temperature::SVector{N}, snow_temperature_before, atmospheric_pressure,
    soil_surface_temperature,
) where N
    maxsnode = state.active_nodes
    (; mean_temperature, mean_temperature_past, phase_heat, layer_mass) = scratch

    if state.current_depth >= 200.0u"cm"
        return (state, snow_temperature, soil_surface_temperature)
    end

    specific_heat = snow_specific_heat(state.density, atmospheric_pressure)
    freeze = u"K"(0.0u"°C")

    # Track which nodes to clamp to 0°C (Fortran: t(cnd)=0, t(cnd+1)=0)
    clamp_to_zero = fill(false, N)
    clamp_soil_surface = false

    phase_heat .= 0.0u"J/m^2"
    layer_mass .= 0.0u"kg/m^2"
    for j in 1:(maxsnode + 1)
        j <= N || continue
        (; node_index, thickness) = snow_layer_geometry(snow_model, state, scratch, j)
        mass = uconvert(u"kg/m^2", thickness * state.density)

        if mean_temperature_past[node_index] >= 0.0u"°C" && mean_temperature[node_index] <= 0.0u"°C"
            layer_mass[node_index] = max(0.0u"kg/m^2", mass)
            phase_heat[node_index] = uconvert(u"J/m^2",
                (mean_temperature_past[node_index] - mean_temperature[node_index]) * specific_heat * layer_mass[node_index])
            # Fortran OSUB.f lines 876-895: set t(cnd)=0, t(cnd+1)=0
            clamp_to_zero[node_index] = true
            if node_index < N
                clamp_to_zero[node_index + 1] = true
            else
                # cnd+1 is the first soil node
                clamp_soil_surface = true
            end
        end
    end

    # Apply temperature clamping to 0°C
    snow_temperature = SVector(ntuple(Val(N)) do i
        clamp_to_zero[i] ? freeze : snow_temperature[i]
    end)
    if clamp_soil_surface
        soil_surface_temperature = freeze
    end

    new_sum_phase = state.sum_phase + sum(phase_heat)
    total_mass = sum(layer_mass)

    # Fortran OSUB.f lines 914-931: if sumphase exceeds latent heat budget,
    # set nodes AND their neighbors to -0.5°C, reset individual qphase entries,
    # and reset phase accumulator.
    if total_mass > 0.0u"kg/m^2" && new_sum_phase > uconvert(u"J/m^2", LATENT_HEAT_FUSION * total_mass)
        freeze_overshoot = u"K"(-0.5u"°C")
        # Mark nodes that overshoot AND their neighbors (Fortran: t(cnd)=-0.5, t(cnd+1)=-0.5)
        overshoot_node = fill(false, N)
        for i in 1:N
            if snow_temperature[i] < freeze - 1e-8u"K" && phase_heat[i] > 0.0u"J/m^2"
                overshoot_node[i] = true
                phase_heat[i] = 0.0u"J/m^2"  # Fortran: qphase(cnd)=0
                if i < N
                    overshoot_node[i + 1] = true
                end
            end
        end
        snow_temperature = SVector(ntuple(Val(N)) do i
            overshoot_node[i] ? freeze_overshoot : snow_temperature[i]
        end)
        new_sum_phase = 0.0u"J/m^2"
    end

    return (setproperties(state, (; sum_phase=new_sum_phase)), snow_temperature, soil_surface_temperature)
end

# ── Snow mass balance ────────────────────────────────────────────────────

update_snow(::NoSnow, args...; kw...) = (nothing, 0.0u"cm", 0.0u"cm")

function update_snow(snow_model::SnowModel{N}, state::SnowState, scratch,
    snow_temperature, snow_temperature_before, environment_instant, step,
    soil_surface_temperature, soil_surface_temperature_before,
    Q_evaporation;
    hourly_rainfall=false, shade=0.0, day_of_year=1,
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
    # Fortran OSUB.f lines 782-817: use Q_evaporation from the energy balance ODE,
    # convert to mass via latent heat of vaporisation, zero if soil surface ≤ 0°C.
    sublimation = 0.0u"cm"
    if u"°C"(soil_surface_temperature) > 0.0u"°C" && Q_evaporation > 0.0u"W/m^2"
        latent_heat = enthalpy_of_vaporisation(soil_surface_temperature)
        gwsurf = Q_evaporation / latent_heat  # kg/s/m²
        sublimation = uconvert(u"cm", gwsurf * 1.0u"hr" / new_dens)
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
    # Fortran OSUB.f lines 836-837: only compute thermal melt if prevsnow >= minsnow
    state_for_melt = setproperties(state, (; density=new_dens))
    previous_depth = step > 1 ? snow_depth_hourly[step - 1] : 0.0u"cm"
    thermal_melt = if is_first_step || previous_depth < snow_model.min_snow_depth
        0.0u"cm"
    else
        raw_melt = snow_thermal_melt(snow_model, state_for_melt, scratch, snow_temperature, snow_temperature_before,
            environment_instant.atmospheric_pressure,
            soil_surface_temperature, soil_surface_temperature_before)
        # Fortran OSUB.f lines 938-940: melt factor only applied when DOY > 1
        day_of_year > 1 ? raw_melt * snow_model.snow_melt_factor : raw_melt
    end
    cumulative_melt = state.cumulative_melt + thermal_melt
    # Fortran OSUB.f lines 957-962: at end of day (siout(1)==1380), clamp cummelted
    # to previous snow depth, then immediately reset to zero. The clamp is dead code
    # in Fortran too (overwritten on the next line). Possibly a bug — the clamped
    # value may have been intended for output.
    if mod(step, 24) == 0
        # cumulative_melt = min(cumulative_melt, previous_depth)
        cumulative_melt = 0.0u"cm"
    end

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
    return (new_state, snowfall, thermal_melt)
end

# ── Adjust snow depth near node boundaries ────────────────────────────────

adjust_snow_near_nodes(::NoSnow, args...) = nothing
function adjust_snow_near_nodes(snow_model::SnowModel, state::SnowState, scratch, snow_temperature, step)
    (; snow_depth_hourly) = scratch

    snowout = snow_depth_hourly[step]
    maxsnode = state.active_nodes
    threshold = snow_model.snow_node_thresholds[max(1, maxsnode)]
    extra_snow = state.extra_snow

    # First adjustment
    if step > 1 && snowout > 0.0u"cm"
        if (snowout - threshold < 0.5u"cm") && (snowout <= snow_depth_hourly[step - 1])
            extra_snow += uconvert(u"kg/m^2", (snowout - threshold) * state.density)
            snow_depth_hourly[step] = max(0.0u"cm", threshold - 0.1u"cm")
        end
    end

    # Fortran OSUB.f: call snowlayer between adjustments to recompute node thresholds
    state = setproperties(state, (; extra_snow))
    state = activate_snow_nodes(snow_model, state, scratch, snow_temperature, step)

    # Second adjustment — use recomputed threshold
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
