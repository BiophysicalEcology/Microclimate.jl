"""
Spatial solve loop for grid-based microclimate simulation.

Calls existing point physics for each cell, then runs lateral flows
for surface water and cold air.
"""

using Base.Threads: @threads
using Unitful: Unitful, ustrip, @u_str

"""
    SpatialMicroProblem

Problem definition for spatial microclimate simulation.

# Fields
- `terrain`: SpatialMicroTerrain with elevation and derived metrics
- `soil_params`: Soil parameters (can be spatially varying or uniform)
- `weather`: Weather forcing (spatially uniform for now)
- `depths`: Soil layer depths
- `tspan`: Time span (start, end) with time units
- `timestep`: Timestep with time units
- `klam21`: KLAM21 cold air drainage model parameters
- `infiltration`: Infiltration method
"""
@kwdef struct SpatialMicroProblem{T<:SpatialMicroTerrain,S,W,D,TS,DT,K<:KLAM21,I<:InfiltrationMethod}
    terrain::T
    soil_params::S
    weather::W
    depths::D
    tspan::TS = (0.0u"hr", 24.0u"hr")
    timestep::DT = 1.0u"hr"
    klam21::K = KLAM21()
    infiltration::I = SimpleInfiltration()
end

"""
    solve(problem::SpatialMicroProblem, state::SpatialMicroState; kw...)

Run spatial microclimate simulation.

Each timestep:
1. Run point physics for each cell (threaded) - produces P (heat loss rate)
2. Run lateral surface water flow
3. Run KLAM21 cold air drainage (uses P from point physics)
4. Apply cold air coupling to surface temperature
"""
function solve(
    problem::SpatialMicroProblem,
    state::SpatialMicroState;
    lateral_substeps::Int = 10,
    surface_water_friction = 0.1,  # scalar or grid
    surface_water_infiltration = 0.0,  # scalar or grid, during lateral flow
    roughness_length = 0.05u"m",  # z₀ for KLAM21 friction, scalar or grid
    verbose::Bool = false,
)
    (; terrain, soil_params, weather, depths, tspan, timestep, klam21, infiltration) = problem
    nx, ny, nz = size(state)

    t = tspan[1]
    step = 0

    # Preallocate P (heat loss rate) grid with units
    P = fill(0.0u"W/m^2", nx, ny)

    while t < tspan[2]
        step += 1
        verbose && println("Step $step, t = $t")

        # Get weather at current time
        # TODO: interpolate weather properly
        current_weather = weather

        # 1. Run point physics for each cell
        # For now, just update surface water from precipitation
        # TODO: call full soil_energy_balance and soil_water_balance! per cell
        # The energy balance would compute P (heat loss rate) for KLAM21
        _add_precipitation!(state, current_weather)

        # 1b. Infiltration from surface pools to soil
        apply_infiltration!(state, depths, infiltration, timestep)

        # 1c. Compute heat loss rate P for KLAM21
        # P = net longwave radiation loss, depends on surface temp and sky view
        _compute_heat_loss!(P, state, terrain, current_weather)

        # 2. Run lateral surface water flow
        if any(x -> ustrip(u"m", x) > 0, state.surface_water)
            _lateral_surface_water!(state, terrain;
                substeps = lateral_substeps,
                friction = surface_water_friction,
                infiltration_rate = surface_water_infiltration)
        end

        # 3. Run KLAM21 cold air drainage (only at night)
        if _is_nighttime(current_weather)
            klam21_step!(klam21, state.cold_air, P, terrain.dem, roughness_length;
                         timestep = timestep,
                         cellsize = terrain.cellsize)
        end

        # 4. Apply cold air coupling to surface temperature
        apply_cold_air_coupling!(state, KLAM21Coupling())

        t += timestep
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
    precip = _get_value(weather, :rainfall, 0.0u"m")

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

_get_cell_value(x::Union{Real,Unitful.Quantity}, i, j) = x
_get_cell_value(x::AbstractMatrix, i, j) = x[i, j]

"""
Run lateral surface water flow using Stencils.jl.
"""
function _lateral_surface_water!(
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain;
    substeps::Int,
    friction,
    infiltration_rate = 0.0,
)
    # Run flow routing - parameters passed to function, not struct
    new_water = surface_water_flow(
        SurfaceWaterFlow(),
        terrain.dem,
        state.surface_water;  # current pool as "precipitation"
        infiltration_rate = infiltration_rate,  # already handled in point physics
        friction = friction,
        timesteps = substeps,
        cellsize = terrain.cellsize,
    )

    # Update state
    state.surface_water .= new_water
end

# =============================================================================
# Heat Loss Calculation for KLAM21
# =============================================================================

const STEFAN_BOLTZMANN = 5.670374419e-8u"W/m^2/K^4"

"""
    _compute_heat_loss!(P, state, terrain, weather)

Compute heat loss rate P (W/m²) for KLAM21 cold air drainage.

Heat loss is computed from net longwave radiation:
- Surface emits: ε σ T_surface⁴
- Surface receives: SVF × σ T_sky⁴ + (1-SVF) × σ T_surface⁴

Where SVF is sky view factor. Areas with high SVF lose more heat.

Returns P > 0 when surface is cooling (typical nighttime).
"""
function _compute_heat_loss!(
    P::AbstractMatrix,
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain,
    weather,
)
    nx, ny, _ = size(state)
    T_air = _get_value(weather, :air_temperature, 288.0u"K")
    emissivity = 0.96  # typical soil/vegetation

    σ = STEFAN_BOLTZMANN

    @threads for j in 1:ny
        for i in 1:nx
            T_surface = state.soil_temperature[i, j, 1]
            svf = terrain.sky_view_factor[i, j]
            T_air_cell = _get_cell_value(T_air, i, j)

            # Estimate sky temperature (simplified: T_sky ≈ T_air - 20K for clear sky)
            # In reality depends on humidity, clouds
            T_sky = T_air_cell - 20.0u"K"

            # Net longwave radiation (positive = cooling)
            # Surface emits to sky (scaled by SVF) and receives from sky
            emitted = emissivity * σ * T_surface^4
            received_sky = emissivity * σ * T_sky^4 * svf
            received_terrain = emissivity * σ * T_air_cell^4 * (1 - svf)

            # P = net heat loss (positive when cooling)
            net_loss = emitted - received_sky - received_terrain
            P[i, j] = max(0.0u"W/m^2", net_loss)
        end
    end
end

"""
    _is_nighttime(weather)

Check if it's nighttime based on solar radiation or explicit flag.
Returns true if global radiation is zero or negligible.
"""
function _is_nighttime(weather)
    solar = _get_value(weather, :global_radiation, 0.0u"W/m^2")
    threshold = 10.0u"W/m^2"
    if solar isa Unitful.Quantity
        return solar < threshold
    elseif solar isa Real
        return solar < 10.0
    else
        # Grid - check if all values are below threshold
        return all(x -> ustrip(u"W/m^2", x) < 10.0, solar)
    end
end
