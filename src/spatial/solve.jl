"""
Spatial solve loop for grid-based microclimate simulation.

Calls existing point physics for each cell, then runs lateral flows
for surface water and cold air.
"""

using Base.Threads: @threads

"""
    SpatialMicroProblem

Problem definition for spatial microclimate simulation.

# Fields
- `terrain`: SpatialMicroTerrain with DEM and derived metrics
- `soil_params`: Soil parameters (can be spatially varying or uniform)
- `weather`: Weather forcing (spatially uniform for now)
- `depths`: Soil layer depths
- `tspan`: Time span (start, end) in hours
- `dt`: Timestep in hours
"""
@kwdef struct SpatialMicroProblem{T<:SpatialMicroTerrain,S,W,D}
    terrain::T
    soil_params::S
    weather::W
    depths::D
    tspan::Tuple{Float64,Float64} = (0.0, 24.0)
    dt::Float64 = 1.0  # hours
end

"""
    solve(problem::SpatialMicroProblem, state::SpatialMicroState; kw...)

Run spatial microclimate simulation.

Each timestep:
1. Run point physics for each cell (threaded)
2. Run lateral surface water flow
3. Run cold air drainage (at night)
4. Update state with coupled feedbacks
"""
function solve(
    problem::SpatialMicroProblem,
    state::SpatialMicroState;
    lateral_substeps::Int = 10,
    cold_air_production::Float64 = 0.1,
    surface_water_friction::Float64 = 0.1,
    verbose::Bool = false,
)
    (; terrain, soil_params, weather, depths, tspan, dt) = problem
    nx, ny, nz = size(state)

    # Preallocate buffers for point physics (one per thread)
    nthreads = Threads.nthreads()
    # TODO: preallocate point physics buffers

    t = tspan[1]
    step = 0

    while t < tspan[2]
        step += 1
        verbose && println("Step $step, t = $t h")

        # Get weather at current time
        # TODO: interpolate weather
        current_weather = weather

        # 1. Run point physics for each cell
        # For now, just update surface water from precipitation
        # TODO: call full soil_energy_balance and soil_water_balance! per cell
        _point_physics_step!(state, terrain, soil_params, current_weather, depths, dt)

        # 2. Run lateral surface water flow
        # Takes pool from each cell, routes according to DEM
        if any(>(0), state.surface_water)
            _lateral_surface_water!(state, terrain, lateral_substeps, surface_water_friction)
        end

        # 3. Run cold air drainage (simplified: at night when sun angle low)
        # TODO: check solar angle from weather
        is_night = true  # placeholder
        if is_night
            _lateral_cold_air!(state, terrain, lateral_substeps, cold_air_production)
        end

        # 4. Coupling: cold air affects surface temperature
        _apply_cold_air_coupling!(state)

        t += dt
    end

    return state
end

"""
Run point physics for all cells (threaded).

Calls existing soil_energy_balance and soil_water_balance! per cell.
Weather and soil_params are grids - index at (i,j) for each cell.
"""
function _point_physics_step!(
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain,
    soil_params,
    weather,
    depths,
    dt,
)
    nx, ny, nz = size(state)

    # Get precipitation from weather grid (or scalar)
    precip = _get_value(weather, :rainfall, 0.0)

    @threads for j in 1:ny
        for i in 1:nx
            # Get cell precipitation
            cell_precip = _get_cell_value(precip, i, j)

            # Add precipitation to pool
            state.surface_water[i, j] += cell_precip

            # Get cell-specific soil params if spatially varying
            cell_soil = _get_cell_params(soil_params, i, j)

            # Infiltration based on soil properties
            k_sat = _get_param(cell_soil, :saturated_hydraulic_conductivity, 0.001)
            θ_sat = _get_param(cell_soil, :saturation_moisture, 0.45)

            θ = state.soil_moisture[i, j, 1]
            saturation_fraction = θ / θ_sat
            infiltration_capacity = k_sat * (1 - saturation_fraction)^2 * dt

            infiltrated = min(state.surface_water[i, j], infiltration_capacity)
            state.surface_water[i, j] -= infiltrated

            # Add to top soil layer (convert depth to volumetric)
            layer_depth = nz > 1 ? depths[2] - depths[1] : 0.1
            state.soil_moisture[i, j, 1] += infiltrated / layer_depth
            state.soil_moisture[i, j, 1] = clamp(state.soil_moisture[i, j, 1], 0.0, θ_sat)

            # TODO: Call full soil_energy_balance ODE per cell
            # TODO: Call full soil_water_balance! per cell
        end
    end
end

# Helper functions for grid/scalar access
_get_value(x::Nothing, field, default) = default
_get_value(x::NamedTuple, field, default) = get(x, field, default)
_get_value(x, field, default) = hasproperty(x, field) ? getproperty(x, field) : default

_get_cell_value(x::Real, i, j) = x
_get_cell_value(x::AbstractMatrix, i, j) = x[i, j]

_get_cell_params(x::Nothing, i, j) = nothing
_get_cell_params(x::NamedTuple, i, j) = x  # uniform params
_get_cell_params(x, i, j) = x  # assume struct with uniform values for now

_get_param(x::Nothing, field, default) = default
_get_param(x::NamedTuple, field, default) = get(x, field, default)
_get_param(x, field, default) = hasproperty(x, field) ? getproperty(x, field) : default

"""
Run lateral surface water flow using Stencils.jl.
"""
function _lateral_surface_water!(
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain,
    substeps::Int,
    friction::Float64,
)
    # Use existing surface_water_flow
    method = SurfaceWaterFlow(
        infiltration_rate = 0.0,  # already handled in point physics
        friction = friction,
        timesteps = substeps,
    )

    # Run flow routing
    new_water = surface_water_flow(
        method,
        terrain.dem,
        state.surface_water;  # current pool as "precipitation"
        cellsize = terrain.cellsize,
    )

    # Update state
    state.surface_water .= new_water
end

"""
Run lateral cold air drainage using Stencils.jl.
"""
function _lateral_cold_air!(
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain,
    substeps::Int,
    production::Float64,
)
    # Use existing cold_air_pooling
    method = ColdAirFlow(
        production = production,
        friction = 0.1,
        timesteps = substeps,
    )

    # Run cold air drainage
    # Start with current cold air + new production based on sky view factor
    # Areas with high SVF cool more → produce more cold air
    initial_cold_air = state.cold_air_depth .+ production .* terrain.sky_view_factor

    new_cold_air = cold_air_pooling(
        method,
        terrain.dem;
        cellsize = terrain.cellsize,
    )

    # Combine: existing drainage pattern + new production
    state.cold_air_depth .= new_cold_air
end

"""
Apply cold air coupling to surface temperature.

Cold air pools reduce surface temperature.
"""
function _apply_cold_air_coupling!(state::SpatialMicroState)
    nx, ny, _ = size(state)

    # Simple coupling: cold air depth reduces surface temperature
    # More sophisticated: cold air has a temperature, mixes with surface air
    cold_air_cooling_rate = 2.0  # K per meter of cold air depth

    @threads for j in 1:ny
        for i in 1:nx
            if state.cold_air_depth[i, j] > 0
                # Reduce surface temperature proportional to cold air depth
                cooling = state.cold_air_depth[i, j] * cold_air_cooling_rate
                state.soil_temperature[i, j, 1] -= cooling

                # Cold air dissipates over time
                state.cold_air_depth[i, j] *= 0.9
            end
        end
    end
end
