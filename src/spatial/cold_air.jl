"""
Cold air pooling models for microclimate simulation.

Based on the analogy between cold air drainage and water flow, cold air
(being denser than warm air) flows downslope and accumulates in topographic
depressions during calm, clear nights.

References:
- Lundquist et al. (2008) "Automated algorithm for mapping regions of cold-air pooling"
- Schwab (2000) "Reliefanalytische Verfahren zur Abschaetzung naechtlicher Kaltluftbewegungen"
- Dietrich & Böhner (2008) "Cold Air Production and Flow in a Low Mountain Range"
"""

using Stencils: Stencils, Moore, StencilArray, SwitchingStencilArray,
    center, neighbors, offsets, mapstencil!
using Geomorphometry: Geomorphometry
using StaticArrays: SVector

abstract type ColdAirPoolingMethod end

"""
    ColdAirFlow(; production=1.0, friction=0.1, timesteps=100)

Stencil-based cold air flow simulation.

Simulates cold air as a fluid that is produced by radiative cooling and flows
downslope, accumulating in depressions. Uses iterative stencil operations with
`SwitchingStencilArray` for efficiency.

# Parameters
- `production::Float64=1.0`: Cold air production per timestep (arbitrary units)
- `friction::Float64=0.1`: Friction coefficient (0-1, higher = slower flow)
- `timesteps::Int=100`: Number of simulation timesteps

# Returns
Cold air depth field (same units as production × timesteps).
"""
Base.@kwdef struct ColdAirFlow <: ColdAirPoolingMethod
    production::Float64 = 1.0
    friction::Float64 = 0.1
    timesteps::Int = 100
end

"""
    cold_air_pooling(method, dem; cellsize=(1,1))

Compute cold air pooling for a digital elevation model.

# Arguments
- `method::ColdAirPoolingMethod`: The method to use
- `dem::AbstractMatrix`: Digital elevation model

# Returns
Depends on method - see specific method documentation.
"""
function cold_air_pooling end

function cold_air_pooling(
    method::ColdAirFlow,
    dem::AbstractMatrix{T};
    cellsize = Geomorphometry.cellsize(dem),
) where T
    _cold_air_flow(dem, method.production, method.friction, method.timesteps, cellsize)
end

function _cold_air_flow(
    dem::AbstractMatrix{T},
    production::Float64,
    friction::Float64,
    timesteps::Int,
    cellsize,
) where T
    rows, cols = size(dem)

    # Initialize cold air depth field
    air = zeros(Float32, rows, cols)

    # Precompute distances to neighbors
    δx, δy = Float64.(abs.(cellsize))
    δxy = sqrt(δx^2 + δy^2)

    # Flow coefficient
    flow_coef = Float32(1.0 - friction)
    prod_f32 = Float32(production)

    # Create stencil arrays
    hood = Moore{1}()
    ssa = SwitchingStencilArray(air, hood)
    dem_sa = StencilArray(dem, hood)

    # Precompute distance for each offset in Moore neighborhood
    # Moore offsets are in row-major order: (-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)
    dists = SVector{8,Float32}(δxy, δy, δxy, δx, δx, δxy, δy, δxy)

    for _ in 1:timesteps
        # Add production to source array
        ssa.source .+= prod_f32

        ssa = mapstencil!(ssa, dem_sa) do air_hood, dem_hood
            c_air = center(air_hood)
            c_dem = center(dem_hood)

            air_neighbors = neighbors(air_hood)
            dem_neighbors = neighbors(dem_hood)

            net_flow = zero(Float32)

            for i in eachindex(air_hood)
                n_dem = dem_neighbors[i]
                n_air = air_neighbors[i]
                dist = dists[i]

                Δz = c_dem - n_dem  # positive if we're higher

                if Δz > zero(Δz)
                    # We're higher - cold air flows out
                    slope_factor = Δz / dist
                    outflow = c_air * slope_factor * flow_coef
                    net_flow -= outflow
                elseif Δz < zero(Δz)
                    # Neighbor is higher - cold air flows in
                    slope_factor = -Δz / dist
                    inflow = n_air * slope_factor * flow_coef
                    net_flow += inflow
                end
            end

            # Return new cold air depth (clamp to non-negative)
            max(zero(Float32), c_air + net_flow)
        end
    end

    return Array(ssa)
end
