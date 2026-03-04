"""
Surface water flow models for microclimate simulation.

Surface water flows downslope following terrain, with losses to infiltration
and evaporation, accumulating in topographic depressions.

Similar to cold air flow but with:
- Precipitation as input (instead of radiative cooling)
- Infiltration as primary sink (soil-dependent)
- Evaporation as secondary sink
"""

using Stencils: Stencils, Moore, StencilArray,
    center, neighbors, offsets, mapstencil!, scatterstencil!, Remove, Wrap
using Geomorphometry: Geomorphometry
using StaticArrays: SVector

abstract type SurfaceWaterMethod end

"""
    SurfaceWaterFlow(; infiltration_rate=0.001, friction=0.1, timesteps=100)

Stencil-based surface water flow simulation.

Simulates surface water runoff that flows downslope, with losses to
infiltration. Uses iterative stencil operations with `SwitchingStencilArray`.

# Parameters
- `infiltration_rate::Float64=0.001`: Infiltration rate (m/timestep)
- `friction::Float64=0.1`: Flow friction coefficient (0-1)
- `timesteps::Int=100`: Number of simulation timesteps

# Returns
Surface water depth field (m).
"""
Base.@kwdef struct SurfaceWaterFlow <: SurfaceWaterMethod
    infiltration_rate::Float64 = 0.001
    friction::Float64 = 0.1
    timesteps::Int = 100
end

"""
    surface_water_flow(method, dem, precipitation; cellsize=(1,1), soil_moisture=nothing)

Compute surface water distribution for a DEM given precipitation input.

# Arguments
- `method::SurfaceWaterMethod`: The method to use
- `dem::AbstractMatrix`: Digital elevation model (m)
- `precipitation::AbstractMatrix`: Precipitation depth (m) - can be uniform or spatially varying

# Keywords
- `cellsize`: Cell dimensions (dx, dy) in meters
- `soil_moisture`: Optional soil moisture grid for variable infiltration

# Returns
Surface water depth field (m).
"""
function surface_water_flow end

function surface_water_flow(
    method::SurfaceWaterFlow,
    dem::AbstractMatrix{T},
    precipitation::Union{AbstractMatrix,Real};
    cellsize = Geomorphometry.cellsize(dem),
    soil_moisture = nothing,
    boundary = Remove(zero(Float32)),
) where T
    precip_grid = precipitation isa Real ? fill(Float32(precipitation), size(dem)) : Float32.(precipitation)
    _surface_water_flow(
        dem, precip_grid, method.infiltration_rate, method.friction,
        method.timesteps, cellsize, soil_moisture, boundary
    )
end

function _surface_water_flow(
    dem::AbstractMatrix{T},
    precipitation::AbstractMatrix,
    infiltration_rate::Float64,
    friction::Float64,
    timesteps::Int,
    cellsize,
    soil_moisture,
    boundary,
) where T
    rows, cols = size(dem)

    # Initialize surface water with precipitation
    water = Float32.(precipitation)
    water_next = similar(water)

    # Precompute distances
    δx, δy = Float64.(abs.(cellsize))
    δxy = sqrt(δx^2 + δy^2)

    flow_coef = Float32(1.0 - friction)
    infil_rate = Float32(infiltration_rate)

    hood = Moore{1}()
    dists = SVector{8,Float32}(δxy, δy, δxy, δx, δx, δxy, δy, δxy)

    for _ in 1:timesteps
        water_sa = StencilArray(water, hood; boundary)
        dem_sa = StencilArray(dem, hood; boundary)

        # Step 1: Compute residuals (water remaining after outflow)
        # Each cell computes how much it keeps after sending water to lower neighbors
        mapstencil!(water_next, water_sa, dem_sa) do water_hood, dem_hood
            c_water = center(water_hood)
            c_dem = center(dem_hood)

            # Apply infiltration
            c_water = max(zero(Float32), c_water - infil_rate)

            # Water surface elevation = terrain + water depth
            c_surface = c_dem + c_water

            water_neighbors = neighbors(water_hood)
            dem_neighbors = neighbors(dem_hood)

            # Compute outflows limited by equilibration
            # Don't send more than half the surface difference to each neighbor
            total_outflow = zero(Float32)
            for k in eachindex(water_hood)
                n_dem = dem_neighbors[k]
                n_water = water_neighbors[k]
                n_surface = n_dem + n_water
                Δz = c_surface - n_surface

                if Δz > zero(Δz)
                    # Limit outflow to half the surface difference (equilibration limit)
                    # Also scale by flow coefficient for gradual flow
                    max_equil_flow = Δz * Float32(0.5) * flow_coef
                    total_outflow += max_equil_flow
                end
            end

            # Limit total outflow to available water
            total_outflow = min(total_outflow, c_water)

            # Return residual (what stays in this cell)
            c_water - total_outflow
        end

        # Step 2: Scatter outflows to neighbors
        # Each cell sends water to lower neighbors based on water surface gradient
        scatterstencil!(+, water_next, water_sa, dem_sa) do water_hood, dem_hood
            c_water = center(water_hood)
            c_dem = center(dem_hood)

            # Same infiltration as above
            c_water = max(zero(Float32), c_water - infil_rate)
            c_surface = c_dem + c_water

            water_neighbors = neighbors(water_hood)
            dem_neighbors = neighbors(dem_hood)

            # Compute outflows limited by equilibration
            outflows_raw = map(eachindex(water_hood)) do k
                n_dem = dem_neighbors[k]
                n_water = water_neighbors[k]
                n_surface = n_dem + n_water
                Δz = c_surface - n_surface

                if Δz > zero(Δz)
                    # Limit outflow to half the surface difference
                    Δz * Float32(0.5) * flow_coef
                else
                    zero(Float32)
                end
            end

            # Scale outflows to not exceed available water
            total = sum(outflows_raw)
            if total > c_water && total > zero(Float32)
                scale = c_water / total
                map(o -> o * scale, outflows_raw)
            else
                outflows_raw
            end
        end

        # Swap buffers
        water, water_next = water_next, water
    end

    return water
end

"""
    surface_water_event(dem, precipitation; kwargs...)

Convenience function for simulating a single precipitation event.

# Arguments
- `dem`: Digital elevation model
- `precipitation`: Total precipitation (m), uniform or grid

# Keywords
- `infiltration_rate=0.001`: Infiltration rate (m/timestep)
- `friction=0.1`: Flow friction
- `duration=100`: Event duration in timesteps
- `cellsize=(1,1)`: Cell size

# Returns
Final surface water depth distribution.
"""
function surface_water_event(
    dem::AbstractMatrix,
    precipitation::Union{AbstractMatrix,Real};
    infiltration_rate = 0.001,
    friction = 0.1,
    duration = 100,
    cellsize = Geomorphometry.cellsize(dem),
)
    method = SurfaceWaterFlow(;
        infiltration_rate,
        friction,
        timesteps = duration,
    )
    surface_water_flow(method, dem, precipitation; cellsize)
end

# =============================================================================
# Infiltration Methods
# =============================================================================

"""
    InfiltrationMethod

Abstract type for methods that transfer surface water to soil.

Infiltration depends on soil properties, current moisture, and ponded water depth.
"""
abstract type InfiltrationMethod end

"""
    SimpleInfiltration(; k_sat=0.001, θ_sat=0.45)

Simple infiltration model based on saturated hydraulic conductivity.

Infiltration capacity decreases quadratically as soil approaches saturation:

    capacity = k_sat * (1 - θ/θ_sat)^2 * dt

Actual infiltration is limited by available ponded water.

# Parameters
- `k_sat::Float64=0.001`: Saturated hydraulic conductivity (m/hour)
- `θ_sat::Float64=0.45`: Saturation volumetric moisture content
"""
Base.@kwdef struct SimpleInfiltration <: InfiltrationMethod
    k_sat::Float64 = 0.001
    θ_sat::Float64 = 0.45
end

"""
    apply_infiltration!(state, depths, method::InfiltrationMethod, dt)

Transfer water from surface pools to soil based on infiltration capacity.

Modifies `state.surface_water` and `state.soil_moisture` in place.
Returns total infiltrated depth across all cells.
"""
function apply_infiltration! end

function apply_infiltration!(
    state::SpatialMicroState,
    depths,
    method::SimpleInfiltration,
    dt::Float64,
)
    nx, ny, nz = size(state)
    (; k_sat, θ_sat) = method
    total_infiltrated = Threads.Atomic{Float64}(0.0)

    # Compute top layer depth
    layer_depth = nz > 1 ? depths[2] - depths[1] : 0.1

    @threads for j in 1:ny
        for i in 1:nx
            pool = state.surface_water[i, j]
            if pool > 0
                θ = state.soil_moisture[i, j, 1]
                saturation_fraction = θ / θ_sat
                infiltration_capacity = k_sat * (1 - saturation_fraction)^2 * dt

                infiltrated = min(pool, infiltration_capacity)
                state.surface_water[i, j] -= infiltrated
                state.soil_moisture[i, j, 1] += infiltrated / layer_depth
                state.soil_moisture[i, j, 1] = clamp(state.soil_moisture[i, j, 1], 0.0, θ_sat)

                Threads.atomic_add!(total_infiltrated, infiltrated)
            end
        end
    end

    return total_infiltrated[]
end

# =============================================================================
# Evaporation from Surface Pools
# =============================================================================

"""
    apply_evaporation!(state, potential_evap, dt)

Evaporate water from surface pools at the potential rate.

Ponded water evaporates like open water (at potential evaporation rate).
Actual evaporation is limited by available pool depth.

# Arguments
- `state::SpatialMicroState`: Spatial state to modify
- `potential_evap::AbstractMatrix`: Potential evaporation rate (m/hour)
- `dt::Float64`: Timestep (hours)

Returns total evaporated depth across all cells.
"""
function apply_evaporation!(
    state::SpatialMicroState,
    potential_evap::AbstractMatrix,
    dt::Float64,
)
    nx, ny, _ = size(state)
    total_evaporated = Threads.Atomic{Float64}(0.0)

    @threads for j in 1:ny
        for i in 1:nx
            pool = state.surface_water[i, j]
            if pool > 0
                evap = min(pool, potential_evap[i, j] * dt)
                state.surface_water[i, j] -= evap
                Threads.atomic_add!(total_evaporated, evap)
            end
        end
    end

    return total_evaporated[]
end

# =============================================================================
# Soil Saturation Check
# =============================================================================

"""
    check_saturation!(state; θ_sat=0.45, layer_depth=0.1)

Check for soil saturation and generate return flow.

When soil moisture exceeds saturation, excess water returns to surface pool.

# Arguments
- `state::SpatialMicroState`: Spatial state to modify
- `θ_sat::Float64=0.45`: Saturation volumetric moisture content
- `layer_depth::Float64=0.1`: Depth of each soil layer (m)
"""
function check_saturation!(
    state::SpatialMicroState;
    θ_sat::Float64 = 0.45,
    layer_depth::Float64 = 0.1,
)
    nx, ny, nz = size(state)

    @threads for j in 1:ny
        for i in 1:nx
            for k in 1:nz
                θ = state.soil_moisture[i, j, k]
                if θ > θ_sat
                    excess = (θ - θ_sat) * layer_depth
                    state.soil_moisture[i, j, k] = θ_sat
                    state.surface_water[i, j] += excess
                end
            end
        end
    end
end
