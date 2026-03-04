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
@kwdef struct SpatialMicroProblem{T<:SpatialMicroTerrain,S,W,D,C<:ColdAirCouplingMethod,I<:InfiltrationMethod}
    terrain::T
    soil_params::S
    weather::W
    depths::D
    tspan::Tuple{Float64,Float64} = (0.0, 24.0)
    dt::Float64 = 1.0  # hours
    cold_air_coupling::C = LinearCoolingCoupling()
    infiltration::I = SimpleInfiltration()
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
    (; terrain, soil_params, weather, depths, tspan, dt, cold_air_coupling, infiltration) = problem
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
        _add_precipitation!(state, current_weather)

        # 1b. Infiltration from surface pools to soil
        apply_infiltration!(state, depths, infiltration, dt)

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
        apply_cold_air_coupling!(state, cold_air_coupling)

        t += dt
    end

    return state
end

"""
Add precipitation to surface water pools.

Weather can be a NamedTuple with :rainfall field, or a struct with rainfall property.
Rainfall can be a scalar (uniform) or matrix (spatially varying).
"""
function _add_precipitation!(state::SpatialMicroState, weather)
    nx, ny, _ = size(state)
    precip = _get_value(weather, :rainfall, 0.0)

    @threads for j in 1:ny
        for i in 1:nx
            state.surface_water[i, j] += _get_cell_value(precip, i, j)
        end
    end
end

# Helper functions for grid/scalar access
_get_value(x::Nothing, field, default) = default
_get_value(x::NamedTuple, field, default) = get(x, field, default)
_get_value(x, field, default) = hasproperty(x, field) ? getproperty(x, field) : default

_get_cell_value(x::Real, i, j) = x
_get_cell_value(x::AbstractMatrix, i, j) = x[i, j]

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

