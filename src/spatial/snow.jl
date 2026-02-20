"""
Spatial snow processes for grid-based microclimate simulation.

Includes:
- Snow redistribution by wind (transport from exposed to sheltered areas)
- Gravitational snow redistribution (avalanching on steep slopes)
- Integration of point snow physics with lateral processes

Snow redistribution follows similar patterns to cold air pooling but with
different physics: wind transports snow from windward to leeward slopes,
and steep slopes can shed snow to lower elevations.

References:
- Liston, G.E. and M. Sturm (1998). A snow-transport model for complex terrain.
  Journal of Glaciology, 44(148), 498-516.
- Winstral, A. and D. Marks (2002). Simulating wind fields and snow redistribution
  using terrain-based parameters. Hydrological Processes, 16(18), 3585-3603.
"""

using Stencils: Stencils, Moore, StencilArray, SwitchingStencilArray,
    center, neighbors, offsets, mapstencil!
using Geomorphometry: Geomorphometry
using StaticArrays: SVector

abstract type SnowRedistributionMethod end

"""
    WindDrivenSnowRedistribution(; transport_coefficient=0.1, ...)

Wind-driven snow redistribution model.

Snow is eroded from windward slopes and deposited in sheltered areas.
Transport rate depends on wind speed, fetch distance, and snow properties.

# Parameters
- `transport_coefficient::Float64=0.1`: Base transport rate coefficient
- `critical_wind_speed::Float64=5.0`: Minimum wind speed for transport (m/s)
- `deposition_coefficient::Float64=0.3`: Rate of deposition in sheltered areas
- `maximum_erosion_fraction::Float64=0.1`: Max fraction eroded per timestep
- `timesteps::Int=10`: Number of redistribution substeps
"""
Base.@kwdef struct WindDrivenSnowRedistribution <: SnowRedistributionMethod
    transport_coefficient::Float64 = 0.1
    critical_wind_speed::Float64 = 5.0       # m/s
    deposition_coefficient::Float64 = 0.3
    maximum_erosion_fraction::Float64 = 0.1
    timesteps::Int = 10
end

"""
    GravitationalSnowRedistribution(; critical_slope=38.0, ...)

Gravitational snow redistribution (avalanching).

Snow is redistributed from steep slopes to lower elevations when slope
exceeds a critical angle or when snowpack becomes unstable.

# Parameters
- `critical_slope::Float64=38.0`: Slope angle triggering redistribution (degrees)
- `redistribution_fraction::Float64=0.5`: Fraction of excess snow redistributed
- `timesteps::Int=5`: Number of redistribution substeps
"""
Base.@kwdef struct GravitationalSnowRedistribution <: SnowRedistributionMethod
    critical_slope::Float64 = 38.0           # degrees
    redistribution_fraction::Float64 = 0.5
    timesteps::Int = 5
end

"""
    snow_redistribution(method, dem, snow_depth, wind_speed, wind_direction; cellsize)

Redistribute snow spatially based on wind or gravity.

# Arguments
- `method::SnowRedistributionMethod`: Redistribution method
- `dem::AbstractMatrix`: Digital elevation model (m)
- `snow_depth::AbstractMatrix`: Current snow depth field (m)
- `wind_speed::Real`: Wind speed (m/s)
- `wind_direction::Real`: Wind direction (degrees from north, clockwise)

# Returns
New snow depth field after redistribution (m).
"""
function snow_redistribution end

function snow_redistribution(
    method::WindDrivenSnowRedistribution,
    dem::AbstractMatrix{T},
    snow_depth::AbstractMatrix,
    wind_speed::Real,
    wind_direction::Real;
    cellsize = Geomorphometry.cellsize(dem),
) where T
    if wind_speed < method.critical_wind_speed
        return copy(snow_depth)
    end

    _wind_snow_redistribution(
        dem, snow_depth, wind_speed, wind_direction,
        method.transport_coefficient,
        method.deposition_coefficient,
        method.maximum_erosion_fraction,
        method.timesteps,
        cellsize
    )
end

function _wind_snow_redistribution(
    dem::AbstractMatrix{T},
    snow_depth::AbstractMatrix,
    wind_speed::Real,
    wind_direction::Real,
    transport_coef::Float64,
    deposition_coef::Float64,
    max_erosion::Float64,
    timesteps::Int,
    cellsize,
) where T
    rows, cols = size(dem)

    # Initialize snow field
    snow = Float32.(snow_depth)

    # Wind direction vector (direction wind is coming FROM)
    # Convert from meteorological (from N, clockwise) to mathematical
    wind_rad = deg2rad(270 - wind_direction)
    wind_x = cos(wind_rad)
    wind_y = sin(wind_rad)

    # Precompute distances
    δx, δy = Float64.(abs.(cellsize))
    δxy = sqrt(δx^2 + δy^2)

    hood = Moore{1}()
    ssa = SwitchingStencilArray(snow, hood)
    dem_sa = StencilArray(dem, hood)

    # Neighbor directions (matching Moore neighborhood order)
    # (-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1)
    neighbor_dirs = SVector{8,Tuple{Float32,Float32}}(
        (-1/δxy, -1/δxy), (-1/δy, 0), (-1/δxy, 1/δxy),
        (0, -1/δx), (0, 1/δx),
        (1/δxy, -1/δxy), (1/δy, 0), (1/δxy, 1/δxy)
    )
    dists = SVector{8,Float32}(δxy, δy, δxy, δx, δx, δxy, δy, δxy)

    transport_scale = Float32(transport_coef * wind_speed / 10.0)
    max_erosion_f32 = Float32(max_erosion)
    deposition_f32 = Float32(deposition_coef)

    for _ in 1:timesteps
        ssa = mapstencil!(ssa, dem_sa) do snow_hood, dem_hood
            c_snow = center(snow_hood)
            c_dem = center(dem_hood)

            if c_snow <= 0
                return c_snow
            end

            snow_neighbors = neighbors(snow_hood)
            dem_neighbors = neighbors(dem_hood)

            net_transport = zero(Float32)

            for k in eachindex(snow_hood)
                n_dem = dem_neighbors[k]
                n_snow = snow_neighbors[k]

                # Direction to neighbor
                dir_x, dir_y = neighbor_dirs[k]

                # Dot product with wind direction
                # Positive = downwind, negative = upwind
                wind_alignment = dir_x * wind_x + dir_y * wind_y

                # Elevation difference
                Δz = n_dem - c_dem

                if wind_alignment > 0.3  # Downwind neighbor
                    # Sheltering effect: deposition increases with slope
                    shelter_factor = max(0, Δz / dists[k])
                    deposition = c_snow * transport_scale * wind_alignment * deposition_f32 * shelter_factor
                    net_transport -= deposition

                elseif wind_alignment < -0.3  # Upwind neighbor
                    # Erosion from upwind cell toward us
                    if n_snow > 0
                        exposure_factor = max(0, -Δz / dists[k])  # exposed if we're higher
                        erosion = n_snow * transport_scale * (-wind_alignment) * exposure_factor
                        net_transport += erosion
                    end
                end
            end

            # Limit erosion
            if net_transport < 0
                net_transport = max(net_transport, -c_snow * max_erosion_f32)
            end

            max(zero(Float32), c_snow + net_transport)
        end
    end

    return Array(ssa)
end

function snow_redistribution(
    method::GravitationalSnowRedistribution,
    dem::AbstractMatrix{T},
    snow_depth::AbstractMatrix;
    cellsize = Geomorphometry.cellsize(dem),
) where T
    _gravitational_snow_redistribution(
        dem, snow_depth,
        method.critical_slope,
        method.redistribution_fraction,
        method.timesteps,
        cellsize
    )
end

function _gravitational_snow_redistribution(
    dem::AbstractMatrix{T},
    snow_depth::AbstractMatrix,
    critical_slope::Float64,
    redistribution_fraction::Float64,
    timesteps::Int,
    cellsize,
) where T
    rows, cols = size(dem)

    snow = Float32.(snow_depth)

    δx, δy = Float64.(abs.(cellsize))
    δxy = sqrt(δx^2 + δy^2)

    critical_tan = Float32(tand(critical_slope))
    redist_frac = Float32(redistribution_fraction)

    hood = Moore{1}()
    ssa = SwitchingStencilArray(snow, hood)
    dem_sa = StencilArray(dem, hood)

    dists = SVector{8,Float32}(δxy, δy, δxy, δx, δx, δxy, δy, δxy)

    for _ in 1:timesteps
        ssa = mapstencil!(ssa, dem_sa) do snow_hood, dem_hood
            c_snow = center(snow_hood)
            c_dem = center(dem_hood)

            snow_neighbors = neighbors(snow_hood)
            dem_neighbors = neighbors(dem_hood)

            # Find steepest downhill neighbor
            max_slope = zero(Float32)
            total_outflow = zero(Float32)

            for k in eachindex(snow_hood)
                n_dem = dem_neighbors[k]
                Δz = c_dem - n_dem  # positive if we're higher

                if Δz > 0
                    slope = Δz / dists[k]
                    if slope > critical_tan
                        # Redistribute proportional to excess slope
                        excess = (slope - critical_tan) / critical_tan
                        outflow = c_snow * redist_frac * min(excess, 1.0f0)
                        total_outflow += outflow
                    end
                end
            end

            # Limit outflow to available snow
            total_outflow = min(total_outflow, c_snow * 0.5f0)

            # Also receive snow from uphill
            inflow = zero(Float32)
            for k in eachindex(snow_hood)
                n_dem = dem_neighbors[k]
                n_snow = snow_neighbors[k]
                Δz = n_dem - c_dem  # positive if neighbor is higher

                if Δz > 0 && n_snow > 0
                    slope = Δz / dists[k]
                    if slope > critical_tan
                        excess = (slope - critical_tan) / critical_tan
                        inflow += n_snow * redist_frac * min(excess, 1.0f0) / 8.0f0
                    end
                end
            end

            max(zero(Float32), c_snow - total_outflow + inflow)
        end
    end

    return Array(ssa)
end

# ============================================================================
# Snow coupling functions for spatial simulation
# ============================================================================

"""
    update_snow_from_precipitation!(state, terrain, precipitation, air_temperature, model)

Update snow state grid from precipitation, partitioning rain vs snow.

Rain passes through to surface water, snow accumulates.
"""
function update_snow_from_precipitation!(
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain,
    precipitation::Union{AbstractMatrix,Real},
    air_temperature::Union{AbstractMatrix,Real},
    model::SnowModel,
)
    nx, ny = size(terrain)

    Threads.@threads for j in 1:ny
        for i in 1:nx
            precip = precipitation isa Real ? precipitation : precipitation[i, j]
            T_air = air_temperature isa Real ? air_temperature : air_temperature[i, j]

            if precip > 0
                (; snowfall, rainfall) = snow_accumulation(model, precip, T_air)

                # Add snowfall to snowpack
                if snowfall > 0
                    _accumulate_snow!(state, i, j, snowfall, model)
                end

                # Rain goes to surface water (or through snowpack to surface)
                state.surface_water[i, j] += rainfall
            end
        end
    end
end

function _accumulate_snow!(state::SpatialMicroState, i::Int, j::Int, snowfall::Real, model::SnowModel)
    fresh_density = _fresh_snow_density(model)

    # Current state
    swe_old = state.snow_water_equivalent[i, j]
    ρ_old = state.snow_density[i, j]

    # New SWE
    swe_new = swe_old + snowfall
    state.snow_water_equivalent[i, j] = swe_new

    # Update density (mass-weighted average)
    if swe_old > 0
        state.snow_density[i, j] = (swe_old * ρ_old + snowfall * fresh_density) / swe_new
    else
        state.snow_density[i, j] = fresh_density
    end

    # Update depth
    state.snow_depth[i, j] = swe_new * 1000.0 / state.snow_density[i, j]

    # Reset age and albedo
    state.snow_age[i, j] = 0.0
    state.snow_albedo[i, j] = _fresh_snow_albedo(model)
end

_fresh_snow_density(::GenericSnowModel) = 100.0
_fresh_snow_density(::Snow17) = 100.0
_fresh_snow_density(model::UtahEnergyBalance) = model.fresh_snow_density

_fresh_snow_albedo(::GenericSnowModel) = 0.85
_fresh_snow_albedo(::Snow17) = 0.85
_fresh_snow_albedo(model::UtahEnergyBalance) = model.fresh_snow_albedo

"""
    update_snow_melt!(state, terrain, forcing, model, timestep_hours)

Compute snowmelt for each cell and route melt to surface water.
"""
function update_snow_melt!(
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain,
    air_temperature::Union{AbstractMatrix,Real},
    model::GenericSnowModel,
    timestep_hours::Real,
)
    nx, ny = size(terrain)

    Threads.@threads for j in 1:ny
        for i in 1:nx
            swe = state.snow_water_equivalent[i, j]
            if swe > 0
                T_air = air_temperature isa Real ? air_temperature : air_temperature[i, j]

                melt = snow_melt(model, swe, T_air, timestep_hours)

                # Remove melt from snowpack
                state.snow_water_equivalent[i, j] = max(0.0, swe - melt)

                # Update depth
                if state.snow_water_equivalent[i, j] > 0
                    state.snow_depth[i, j] = state.snow_water_equivalent[i, j] * 1000.0 / state.snow_density[i, j]
                else
                    state.snow_depth[i, j] = 0.0
                end

                # Route melt to surface water
                state.surface_water[i, j] += melt

                # Age the snow
                state.snow_age[i, j] += timestep_hours / 24.0
            end
        end
    end
end

"""
    update_snow_albedo!(state, model, timestep_hours)

Decay snow albedo with age.
"""
function update_snow_albedo!(
    state::SpatialMicroState,
    model::UtahEnergyBalance,
    timestep_hours::Real,
)
    nx, ny, _ = size(state)

    decay_rate = model.albedo_decay_rate
    α_old = model.old_snow_albedo
    dt_days = timestep_hours / 24.0

    Threads.@threads for j in 1:ny
        for i in 1:nx
            if state.snow_water_equivalent[i, j] > 0
                α = state.snow_albedo[i, j]
                state.snow_albedo[i, j] = α_old + (α - α_old) * exp(-decay_rate * dt_days)
            end
        end
    end
end

"""
    update_surface_temperature_for_snow!(state, terrain)

Modify surface soil temperature to account for snow insulation.

Snow acts as an insulating layer, decoupling soil from atmospheric forcing.
"""
function update_surface_temperature_for_snow!(
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain,
)
    nx, ny, _ = size(state)

    # Snow insulation: soil temperature relaxes toward 0°C under deep snow
    insulation_depth = 0.3  # m, characteristic depth for full insulation

    Threads.@threads for j in 1:ny
        for i in 1:nx
            snow_depth = state.snow_depth[i, j]
            if snow_depth > 0.01  # > 1cm
                T_soil = state.soil_temperature[i, j, 1]
                T_snow_base = 273.15  # K, snow base at 0°C when melting

                # Insulation fraction increases with snow depth
                insulation_fraction = 1 - exp(-snow_depth / insulation_depth)

                # Blend soil temperature toward snow base temperature
                state.soil_temperature[i, j, 1] = T_soil * (1 - insulation_fraction) + T_snow_base * insulation_fraction
            end
        end
    end
end

"""
    get_effective_surface_albedo(state, base_albedo, i, j)

Get effective surface albedo accounting for snow cover.
"""
function get_effective_surface_albedo(
    state::SpatialMicroState,
    base_albedo::Real,
    i::Int, j::Int,
)
    snow_fraction = snow_cover_fraction(state, i, j)
    snow_albedo = state.snow_albedo[i, j]

    return base_albedo * (1 - snow_fraction) + snow_albedo * snow_fraction
end

"""
    apply_cold_air_to_snow!(state, terrain, cold_air_temperature)

Cold air pooling reduces snow melt by lowering local air temperature.

This coupling means valleys with cold air pools retain snow longer.
"""
function apply_cold_air_to_snow!(
    state::SpatialMicroState,
    terrain::SpatialMicroTerrain,
    cold_air_temperature::Real = 268.15,  # K, ~-5°C
)
    nx, ny, _ = size(state)

    mixing_depth = 1.0  # m

    Threads.@threads for j in 1:ny
        for i in 1:nx
            cold_depth = state.cold_air_depth[i, j]
            if cold_depth > 0 && state.snow_water_equivalent[i, j] > 0
                # Cold air reduces snow temperature (slowing melt)
                mixing_fraction = 1 - exp(-cold_depth / mixing_depth)

                T_snow = state.snow_temperature[i, j] + 273.15  # to K
                T_cold = cold_air_temperature

                T_new = T_snow * (1 - mixing_fraction) + T_cold * mixing_fraction
                state.snow_temperature[i, j] = T_new - 273.15  # back to °C

                # Colder snow = increased cold content
                if state.snow_temperature[i, j] < 0
                    # Simplified: add to cold content
                    state.snow_cold_content[i, j] += abs(state.snow_temperature[i, j]) * 0.1
                end
            end
        end
    end
end
