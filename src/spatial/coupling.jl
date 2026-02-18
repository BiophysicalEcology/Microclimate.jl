"""
Coupling between lateral flows and point physics.

Handles feedbacks:
- Surface water (pool) → infiltration → soil moisture
- Cold air depth → surface air temperature → energy balance
- Soil saturation → reduced infiltration capacity
- Wet surface → albedo change, increased evaporation
"""

using Base.Threads: @threads

"""
    update_infiltration_from_pool!(state, soil_params, dt)

Transfer water from surface pool to soil based on infiltration capacity.

Infiltration depends on:
- Current soil moisture (saturated soil can't infiltrate)
- Soil hydraulic conductivity
- Pool depth (more water = more pressure)

Returns total infiltrated depth.
"""
function update_infiltration_from_pool!(
    state::SpatialMicroState,
    soil_params,
    dt::Float64,
)
    nx, ny, _ = size(state)
    total_infiltrated = 0.0

    # Get soil parameters
    # TODO: handle spatially varying soil params
    k_sat = 0.001  # m/hour, saturated hydraulic conductivity placeholder
    θ_sat = 0.45   # saturation moisture

    @threads for j in 1:ny
        for i in 1:nx
            pool = state.surface_water[i, j]
            if pool > 0
                θ = state.soil_moisture[i, j, 1]

                # Infiltration capacity decreases as soil saturates
                saturation_fraction = θ / θ_sat
                infiltration_capacity = k_sat * (1 - saturation_fraction)^2 * dt

                # Actual infiltration limited by pool and capacity
                infiltrated = min(pool, infiltration_capacity)

                # Update state
                state.surface_water[i, j] -= infiltrated
                state.soil_moisture[i, j, 1] += infiltrated / 0.1  # crude: 10cm top layer

                total_infiltrated += infiltrated
            end
        end
    end

    return total_infiltrated
end

"""
    update_surface_temperature_from_cold_air!(state, terrain)

Modify surface temperature based on cold air pool depth.

Cold air pooling in valleys reduces nighttime temperatures.
Effect depends on:
- Cold air depth
- Cold air temperature (assumed ~dewpoint)
- Mixing with ambient air
"""
function update_surface_temperature_from_cold_air!(
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain,
)
    nx, ny, _ = size(state)

    # Temperature of cold air (approximation)
    cold_air_temperature = 273.0  # K, near freezing

    # Mixing coefficient: how much cold air affects surface
    mixing_depth = 1.0  # m, characteristic mixing depth

    @threads for j in 1:ny
        for i in 1:nx
            cold_depth = state.cold_air_depth[i, j]
            if cold_depth > 0
                # Mixing fraction increases with cold air depth
                mixing_fraction = 1 - exp(-cold_depth / mixing_depth)

                # Blend surface temperature toward cold air temperature
                T_surface = state.soil_temperature[i, j, 1]
                T_new = T_surface * (1 - mixing_fraction) + cold_air_temperature * mixing_fraction

                state.soil_temperature[i, j, 1] = T_new
            end
        end
    end
end

"""
    update_albedo_from_wetness!(albedo, state)

Modify surface albedo based on surface wetness.

Wet surfaces have lower albedo (darker).
"""
function update_albedo_from_wetness!(
    albedo::AbstractMatrix,
    state::SpatialMicroState,
    dry_albedo::Float64 = 0.3,
    wet_albedo::Float64 = 0.1,
)
    nx, ny, _ = size(state)

    # Wetness affects albedo: wet = darker
    water_depth_for_full_wet = 0.01  # m

    @threads for j in 1:ny
        for i in 1:nx
            pool = state.surface_water[i, j]
            wetness = min(1.0, pool / water_depth_for_full_wet)

            albedo[i, j] = dry_albedo * (1 - wetness) + wet_albedo * wetness
        end
    end
end

"""
    update_evaporation_from_pool!(state, potential_evap, dt)

Evaporate water from surface pools.

Ponded water evaporates at potential rate (like open water).
Reduces pool depth, adds to atmospheric moisture.
"""
function update_evaporation_from_pool!(
    state::SpatialMicroState,
    potential_evap::AbstractMatrix,  # m/hour
    dt::Float64,
)
    nx, ny, _ = size(state)
    total_evaporated = 0.0

    @threads for j in 1:ny
        for i in 1:nx
            pool = state.surface_water[i, j]
            if pool > 0
                evap = min(pool, potential_evap[i, j] * dt)
                state.surface_water[i, j] -= evap
                total_evaporated += evap
            end
        end
    end

    return total_evaporated
end

"""
    check_soil_saturation!(state)

Check for soil saturation and generate return flow.

When soil moisture exceeds saturation, excess water
returns to surface pool (seepage).
"""
function check_soil_saturation!(state::SpatialMicroState, θ_sat::Float64 = 0.45)
    nx, ny, nz = size(state)

    @threads for j in 1:ny
        for i in 1:nx
            for k in 1:nz
                θ = state.soil_moisture[i, j, k]
                if θ > θ_sat
                    # Excess water returns to surface
                    excess = (θ - θ_sat) * 0.1  # crude: 10cm layer depth
                    state.soil_moisture[i, j, k] = θ_sat
                    state.surface_water[i, j] += excess
                end
            end
        end
    end
end
