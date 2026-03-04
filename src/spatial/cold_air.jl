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

# =============================================================================
# Cold Air Coupling Methods
# =============================================================================

"""
    ColdAirCouplingMethod

Abstract type for methods that couple cold air depth to surface temperature.

Cold air pooling affects surface temperature through mixing with ambient air.
Different methods model this coupling with varying physical fidelity.
"""
abstract type ColdAirCouplingMethod end

"""
    ExponentialMixingCoupling(; mixing_depth=1.0, cold_air_temperature=273.0)

Exponential mixing model for cold air coupling.

Surface temperature is blended toward the cold air temperature with a mixing
fraction that increases exponentially with cold air depth:

    mixing_fraction = 1 - exp(-cold_air_depth / mixing_depth)
    T_new = T_surface * (1 - mixing_fraction) + T_cold * mixing_fraction

# Parameters
- `mixing_depth::Float64=1.0`: Characteristic depth for mixing (m)
- `cold_air_temperature::Float64=273.0`: Temperature of cold air pool (K)
"""
Base.@kwdef struct ExponentialMixingCoupling <: ColdAirCouplingMethod
    mixing_depth::Float64 = 1.0
    cold_air_temperature::Float64 = 273.0
end

"""
    LinearCoolingCoupling(; cooling_rate=2.0, dissipation_rate=0.9)

Linear cooling model for cold air coupling.

Surface temperature is reduced proportionally to cold air depth:

    cooling = cold_air_depth * cooling_rate
    T_new = T_surface - cooling

Cold air dissipates each timestep by the dissipation rate.

# Parameters
- `cooling_rate::Float64=2.0`: Cooling per meter of cold air depth (K/m)
- `dissipation_rate::Float64=0.9`: Fraction of cold air retained per timestep
"""
Base.@kwdef struct LinearCoolingCoupling <: ColdAirCouplingMethod
    cooling_rate::Float64 = 2.0
    dissipation_rate::Float64 = 0.9
end

"""
    apply_cold_air_coupling!(state, method::ColdAirCouplingMethod)

Apply cold air coupling to modify surface temperature based on cold air depth.

Modifies `state.soil_temperature[:,:,1]` in place based on `state.cold_air_depth`.
"""
function apply_cold_air_coupling! end

function apply_cold_air_coupling!(
    state::SpatialMicroState,
    method::ExponentialMixingCoupling,
)
    nx, ny, _ = size(state)
    (; mixing_depth, cold_air_temperature) = method

    @threads for j in 1:ny
        for i in 1:nx
            cold_depth = state.cold_air_depth[i, j]
            if cold_depth > 0
                mixing_fraction = 1 - exp(-cold_depth / mixing_depth)
                T_surface = state.soil_temperature[i, j, 1]
                state.soil_temperature[i, j, 1] = T_surface * (1 - mixing_fraction) + cold_air_temperature * mixing_fraction
            end
        end
    end
end

function apply_cold_air_coupling!(
    state::SpatialMicroState,
    method::LinearCoolingCoupling,
)
    nx, ny, _ = size(state)
    (; cooling_rate, dissipation_rate) = method

    @threads for j in 1:ny
        for i in 1:nx
            if state.cold_air_depth[i, j] > 0
                cooling = state.cold_air_depth[i, j] * cooling_rate
                state.soil_temperature[i, j, 1] -= cooling
                state.cold_air_depth[i, j] *= dissipation_rate
            end
        end
    end
end
