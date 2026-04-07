"""
Spatial state for grid-based microclimate simulation.

State variables are stored as grids that can be updated in-place
during simulation. The `pool` from point-based soil_water_balance!
feeds into lateral surface_water_flow.
"""

using Unitful: Unitful, ustrip, @u_str

"""
    KLAM21State

State variables for KLAM21 cold air drainage model.

All fields are 2D arrays with appropriate units.

# Fields
- `heat_deficit`: Heat deficit E (J/m²) - primary state variable
- `depth`: Cold air layer depth H (m) - diagnosed from heat_deficit
- `temperature_disturbance`: Temperature disturbance ΔT (K) - diagnosed from depth, positive = colder than ambient
- `velocity_x`: East-west velocity component (m/s)
- `velocity_y`: North-south velocity component (m/s)
"""
struct KLAM21State{E,H,T,V}
    heat_deficit::E           # J/m²
    depth::H                  # m
    temperature_disturbance::T  # K
    velocity_x::V             # m/s
    velocity_y::V             # m/s
end

"""
    KLAM21State(nx, ny)

Create a KLAM21State with zero-initialized arrays of the given size.
All arrays are initialized with appropriate Unitful types.
"""
function KLAM21State(nx::Int, ny::Int)
    KLAM21State(
        fill(0.0u"J/m^2", nx, ny),    # heat_deficit
        fill(0.0u"m", nx, ny),         # depth
        fill(0.0u"K", nx, ny),         # temperature_disturbance
        fill(0.0u"m/s", nx, ny),       # velocity_x
        fill(0.0u"m/s", nx, ny),       # velocity_y
    )
end

"""
    KLAM21State(nx, ny, ArrayType)

Create a KLAM21State using the specified array constructor.

This allows creating state on GPU by passing a GPU array type:

```julia
using CUDA
state = KLAM21State(100, 100, CuArray)
```

The ArrayType should be a callable that creates an array from an existing array,
e.g., `CuArray`, `ROCArray`, or `MtlArray`.
"""
function KLAM21State(nx::Int, ny::Int, ArrayType)
    KLAM21State(
        ArrayType(fill(0.0u"J/m^2", nx, ny)),    # heat_deficit
        ArrayType(fill(0.0u"m", nx, ny)),         # depth
        ArrayType(fill(0.0u"K", nx, ny)),         # temperature_disturbance
        ArrayType(fill(0.0u"m/s", nx, ny)),       # velocity_x
        ArrayType(fill(0.0u"m/s", nx, ny)),       # velocity_y
    )
end

"""
    SpatialMicroState

Grid-based state for spatial microclimate simulation.

All fields use Unitful types for physical quantities.

# Fields
## 3D fields (nx, ny, nz) - soil layers
- `soil_temperature`: Soil temperature at each depth (K)
- `soil_moisture`: Volumetric soil moisture (m³/m³)
- `soil_water_potential`: Soil water potential (J/kg)
- `accumulated_latent_heat`: Accumulated latent heat for phase transitions (J)

## 2D fields (nx, ny) - surface
- `surface_water`: Pooled water depth from infiltration excess (m)
- `soil_wetness`: Surface wetness fraction for evaporation (0-1, dimensionless)

## Nested state
- `cold_air`: KLAM21State for cold air drainage
"""
struct SpatialMicroState{ST,SM,WP,LH,SW,W,K<:KLAM21State}
    # 3D: soil layers (nx, ny, nz)
    soil_temperature::ST      # K
    soil_moisture::SM         # m³/m³
    soil_water_potential::WP  # J/kg
    accumulated_latent_heat::LH  # J

    # 2D: surface (nx, ny)
    surface_water::SW         # m
    soil_wetness::W           # dimensionless (0-1)

    # Nested: cold air drainage state
    cold_air::K
end

"""
    SpatialMicroState(nx, ny, nz; initial_soil_temperature=280.0u"K", ...)

Create a SpatialMicroState with initialized arrays of the given size.
All arrays are initialized with appropriate Unitful types.
"""
function SpatialMicroState(
    nx::Int, ny::Int, nz::Int;
    initial_soil_temperature = 280.0u"K",
    initial_soil_moisture = 0.3u"m^3/m^3",
    initial_water_potential = -10.0u"J/kg",
)
    SpatialMicroState(
        fill(initial_soil_temperature, nx, ny, nz),
        fill(initial_soil_moisture, nx, ny, nz),
        fill(initial_water_potential, nx, ny, nz),
        fill(0.0u"J", nx, ny, nz),      # accumulated_latent_heat
        fill(0.0u"m", nx, ny),           # surface_water
        zeros(Float64, nx, ny),          # soil_wetness (dimensionless)
        KLAM21State(nx, ny),             # cold_air
    )
end

"""
    SpatialMicroState(elevation; nz=10, kw...)

Create spatial state matching the size of an elevation grid (DEM).
"""
function SpatialMicroState(elevation::AbstractMatrix; nz::Int = 10, kw...)
    nx, ny = size(elevation)
    SpatialMicroState(nx, ny, nz; kw...)
end

Base.size(s::SpatialMicroState) = size(s.soil_temperature)

# Get surface temperature (top soil layer)
surface_temperature(s::SpatialMicroState) = @view s.soil_temperature[:, :, 1]

# Get surface moisture (top soil layer)
surface_moisture(s::SpatialMicroState) = @view s.soil_moisture[:, :, 1]

# Convenience accessor for cold air depth (from nested KLAM21State)
cold_air_depth(s::SpatialMicroState) = s.cold_air.depth
