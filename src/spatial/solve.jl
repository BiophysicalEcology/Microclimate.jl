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
- `snow_model`: Snow model (NoSnow() to disable, or DegreeDaySnow, Snow17, UtahEnergyBalance, KearneySnow)
- `snow_redistribution`: Snow redistribution method (nothing to disable)
"""
@kwdef struct SpatialMicroProblem{T<:SpatialMicroTerrain,S,W,D,SM,SR}
    terrain::T
    soil_params::S
    weather::W
    depths::D
    tspan::Tuple{Float64,Float64} = (0.0, 24.0)
    dt::Float64 = 1.0  # hours
    snow_model::SM = NoSnow()
    snow_redistribution::SR = nothing
end

"""
    solve(problem::SpatialMicroProblem, state::SpatialMicroState; kw...)

Run spatial microclimate simulation.

Each timestep:
1. Snow accumulation and melt (if snow model enabled)
2. Run point physics for each cell (threaded)
3. Run lateral surface water flow (including snowmelt)
4. Run snow redistribution (if enabled)
5. Run cold air drainage (at night)
6. Update state with coupled feedbacks (cold air, snow insulation)
"""
function solve(
    problem::SpatialMicroProblem,
    state::SpatialMicroState;
    lateral_substeps::Int = 10,
    cold_air_production::Float64 = 0.1,
    surface_water_friction::Float64 = 0.1,
    wind_speed::Float64 = 0.0,
    wind_direction::Float64 = 0.0,
    verbose::Bool = false,
)
    (; terrain, soil_params, weather, depths, tspan, dt, snow_model, snow_redistribution) = problem
    nx, ny, nz = size(state)

    t = tspan[1]
    step = 0

    while t < tspan[2]
        step += 1
        verbose && println("Step $step, t = $t h")

        # Get weather at current time
        # TODO: interpolate weather
        current_weather = weather

        # Get weather values for this timestep
        precip = _get_value(current_weather, :rainfall, 0.0)
        air_temp = _get_value(current_weather, :air_temperature, 10.0)

        # 1. Snow accumulation and melt → returns water output (rain + melt)
        water_output = _snow_water_output!(state, terrain, snow_model, precip, air_temp, dt)

        # 2. Route water output to surface water, then run point physics
        _add_water_to_surface!(state, water_output)
        _point_physics_step!(state, terrain, soil_params, current_weather, depths, dt)

        # 3. Run lateral surface water flow (includes snowmelt runoff)
        if any(>(0), state.surface_water)
            _lateral_surface_water!(state, terrain, lateral_substeps, surface_water_friction)
        end

        # 4. Run snow redistribution (if enabled and snow present)
        _lateral_snow!(state, terrain, snow_redistribution, wind_speed, wind_direction)

        # 5. Run cold air drainage (simplified: at night when sun angle low)
        # TODO: check solar angle from weather
        is_night = true  # placeholder
        if is_night
            _lateral_cold_air!(state, terrain, lateral_substeps, cold_air_production)
        end

        # 6. Coupling effects
        _apply_cold_air_coupling!(state)

        # Snow-cold air coupling
        _snow_coupling!(state, terrain, snow_model)

        t += dt
    end

    return state
end

"""
Add water to surface water grid.
"""
function _add_water_to_surface!(state::SpatialMicroState, water::AbstractMatrix)
    state.surface_water .+= water
end

function _add_water_to_surface!(state::SpatialMicroState, water::Real)
    state.surface_water .+= water
end

"""
Compute water output from snow physics. Returns matrix of water (rain + melt) per cell.
Snow state is updated in-place, but surface_water is NOT modified here.
"""
# NoSnow: all precipitation passes through
function _snow_water_output!(
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain,
    ::NoSnow,
    precipitation::Union{AbstractMatrix,Real},
    air_temperature::Union{AbstractMatrix,Real},
    dt::Real,
)
    # Return precipitation unchanged
    if precipitation isa Real
        nx, ny, _ = size(state)
        return fill(precipitation, nx, ny)
    else
        return copy(precipitation)
    end
end

function _snow_water_output!(
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain,
    snow_model::SnowModel,
    precipitation::Union{AbstractMatrix,Real},
    air_temperature::Union{AbstractMatrix,Real},
    dt::Real,
)
    nx, ny, _ = size(state)
    water_output = zeros(nx, ny)

    @threads for j in 1:ny
        for i in 1:nx
            precip = precipitation isa Real ? precipitation : precipitation[i, j]
            T_air = air_temperature isa Real ? air_temperature : air_temperature[i, j]

            # Partition precipitation into rain and snow
            (; snowfall, rainfall) = snow_accumulation(snow_model, precip, T_air)

            # Accumulate snow
            if snowfall > 0
                _accumulate_snow!(state, i, j, snowfall, snow_model)
            end

            # Compute melt
            swe = state.snow_water_equivalent[i, j]
            melt = 0.0
            if swe > 0
                melt = _compute_snow_melt(snow_model, state, i, j, T_air, rainfall, dt)

                # Remove melt from snowpack
                state.snow_water_equivalent[i, j] = max(0.0, swe - melt)

                # Update depth
                if state.snow_water_equivalent[i, j] > 0
                    state.snow_depth[i, j] = state.snow_water_equivalent[i, j] * 1000.0 / state.snow_density[i, j]
                else
                    state.snow_depth[i, j] = 0.0
                end

                # Age the snow
                state.snow_age[i, j] += dt / 24.0
            end

            # Water output = rain + melt
            water_output[i, j] = rainfall + melt
        end
    end

    return water_output
end

# Dispatch for different snow models
function _compute_snow_melt(model::GenericSnowModel, state, i, j, air_temp, rainfall, dt)
    snow_melt(model, state.snow_water_equivalent[i, j], air_temp, dt)
end

function _compute_snow_melt(model::KearneySnow, state, i, j, air_temp, rainfall, dt)
    # Simplified: use degree-day like approach for spatial
    # Full implementation would need energy flux
    T_threshold = 0.0
    if air_temp <= T_threshold
        return 0.0
    end

    # Basic temperature melt
    degree_hours = (air_temp - T_threshold) * dt
    melt_rate = 3.0 / 24.0  # mm/°C/hour
    melt = melt_rate * degree_hours / 1000.0  # convert to m

    # Rain-on-snow melt
    if rainfall > 0 && air_temp > 0
        melt += rainfall * air_temp * model.rain_melt_factor
    end

    return min(melt, state.snow_water_equivalent[i, j])
end

function _compute_snow_melt(model::Snow17, state, i, j, air_temp, rainfall, dt)
    # Simplified Snow17 for spatial
    T_melt = model.melt_base_temperature
    if air_temp <= T_melt
        return 0.0
    end

    mf = (model.maximum_melt_factor + model.minimum_melt_factor) / 2  # average melt factor
    mf_scaled = mf * (dt / 6.0)  # scale to timestep

    melt = mf_scaled * (air_temp - T_melt) / 1000.0  # mm to m

    return min(melt, state.snow_water_equivalent[i, j])
end

function _compute_snow_melt(model::UtahEnergyBalance, state, i, j, air_temp, rainfall, dt)
    # UEB needs full energy balance - simplified here
    T_threshold = 0.0
    if air_temp <= T_threshold
        return 0.0
    end

    # Approximate with temperature index
    melt_rate = 4.0 / 24.0  # mm/°C/hour
    melt = melt_rate * (air_temp - T_threshold) * dt / 1000.0

    return min(melt, state.snow_water_equivalent[i, j])
end

"""
Run lateral snow redistribution.
"""
# No redistribution
function _lateral_snow!(
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain,
    ::Nothing,
    wind_speed::Real,
    wind_direction::Real,
)
    # Do nothing
end

function _lateral_snow!(
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain,
    method::SnowRedistributionMethod,
    wind_speed::Real,
    wind_direction::Real,
)
    # Skip if no snow present
    if !any(>(0), state.snow_water_equivalent)
        return
    end
    new_snow = snow_redistribution(
        method,
        terrain.dem,
        state.snow_depth,
        wind_speed,
        wind_direction;
        cellsize = terrain.cellsize,
    )

    # Update snow depth and recalculate SWE
    for j in axes(state.snow_depth, 2), i in axes(state.snow_depth, 1)
        old_depth = state.snow_depth[i, j]
        new_depth = new_snow[i, j]

        if old_depth > 0 && new_depth > 0
            # Scale SWE proportionally
            ratio = new_depth / old_depth
            state.snow_water_equivalent[i, j] *= ratio
        elseif new_depth > 0 && old_depth == 0
            # Snow deposited where there was none
            # Assume density of surrounding snow or fresh snow
            state.snow_water_equivalent[i, j] = new_depth * state.snow_density[i, j] / 1000.0
        else
            state.snow_water_equivalent[i, j] = 0.0
        end

        state.snow_depth[i, j] = new_depth
    end
end

"""
Snow-soil coupling effects.
"""
# NoSnow: do nothing
_snow_coupling!(state::SpatialMicroState, terrain::SpatialMicroTerrain, ::NoSnow) = nothing

# Active snow model: apply coupling effects
function _snow_coupling!(state::SpatialMicroState, terrain::SpatialMicroTerrain, ::SnowModel)
    apply_cold_air_to_snow!(state, terrain)
    update_surface_temperature_for_snow!(state, terrain)
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
