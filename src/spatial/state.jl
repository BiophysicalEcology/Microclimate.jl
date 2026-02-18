"""
Spatial state for grid-based microclimate simulation.

State variables are stored as grids that can be updated in-place
during simulation. The `pool` from point-based soil_water_balance!
feeds into lateral surface_water_flow.
"""

"""
    SpatialMicroState

Grid-based state for spatial microclimate simulation.

# Fields
## 3D fields (nx, ny, nz) - soil layers
- `soil_temperature`: Soil temperature at each depth (K)
- `soil_moisture`: Volumetric soil moisture (m³/m³)
- `soil_water_potential`: Soil water potential (J/kg)
- `accumulated_latent_heat`: Accumulated latent heat for phase transitions (J)

## 2D fields (nx, ny) - surface
- `surface_water`: Pooled water depth from infiltration excess (m)
- `cold_air_depth`: Cold air pool depth from drainage (m)
- `soil_wetness`: Surface wetness fraction for evaporation (0-1)
"""
struct SpatialMicroState{T,A3<:AbstractArray{T,3},A2<:AbstractArray{T,2}}
    # 3D: soil layers (nx, ny, nz)
    soil_temperature::A3
    soil_moisture::A3
    soil_water_potential::A3
    accumulated_latent_heat::A3

    # 2D: surface (nx, ny)
    surface_water::A2
    cold_air_depth::A2
    soil_wetness::A2
end

function SpatialMicroState(
    nx::Int, ny::Int, nz::Int;
    T::Type = Float64,
    initial_soil_temperature = 280.0,  # K
    initial_soil_moisture = 0.3,       # m³/m³
    initial_water_potential = -10.0,   # J/kg
)
    SpatialMicroState(
        fill(T(initial_soil_temperature), nx, ny, nz),
        fill(T(initial_soil_moisture), nx, ny, nz),
        fill(T(initial_water_potential), nx, ny, nz),
        zeros(T, nx, ny, nz),  # accumulated_latent_heat
        zeros(T, nx, ny),      # surface_water
        zeros(T, nx, ny),      # cold_air_depth
        zeros(T, nx, ny),      # soil_wetness
    )
end

"""
    SpatialMicroState(dem; nz=10, kw...)

Create spatial state matching the size of a DEM.
"""
function SpatialMicroState(dem::AbstractMatrix; nz::Int = 10, kw...)
    nx, ny = size(dem)
    SpatialMicroState(nx, ny, nz; kw...)
end

Base.size(s::SpatialMicroState) = size(s.soil_temperature)

# Get surface temperature (top soil layer)
surface_temperature(s::SpatialMicroState) = @view s.soil_temperature[:, :, 1]

# Get surface moisture (top soil layer)
surface_moisture(s::SpatialMicroState) = @view s.soil_moisture[:, :, 1]
