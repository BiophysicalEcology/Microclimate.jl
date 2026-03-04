"""
Spatial terrain for grid-based microclimate simulation.

Holds DEM and precomputed terrain metrics from Geomorphometry.jl.
These are static (read-only) during simulation.
"""

using Geomorphometry: slope, aspect, sky_view_factor, filldepressions, basin_depth

"""
    SpatialMicroTerrain

Static terrain data for spatial microclimate simulation.

# Fields
- `dem`: Digital elevation model (m)
- `slope`: Terrain slope (radians)
- `aspect`: Terrain aspect (radians from north)
- `sky_view_factor`: Fraction of sky visible (0-1)
- `filled_dem`: Depression-filled DEM for flow routing
- `basin_depth`: Depth below pour point (m) - identifies pooling zones
- `cellsize`: Grid cell dimensions (dx, dy) in meters

# Constructor
    SpatialMicroTerrain(dem; cellsize=(30.0, 30.0))

Computes all derived terrain metrics from the DEM.
"""
struct SpatialMicroTerrain{T,A<:AbstractArray{T,2},B<:AbstractArray{<:Real,2}}
    dem::A
    slope::B
    aspect::B
    sky_view_factor::B
    filled_dem::A
    basin_depth::A
    cellsize::Tuple{T,T}
end

function SpatialMicroTerrain(
    dem::AbstractMatrix{T};
    cellsize::Tuple{Real,Real} = (30.0, 30.0),
) where T
    cs = T.(cellsize)

    # Compute terrain metrics using Geomorphometry.jl
    s = slope(dem; cellsize=cs)
    a = aspect(dem; cellsize=cs)
    svf = sky_view_factor(dem; cellsize=cs)
    filled = filldepressions(dem)
    bd = basin_depth(dem; filled=filled)

    SpatialMicroTerrain(
        dem,
        s,
        a,
        svf,
        filled,
        bd,
        cs,
    )
end

Base.size(t::SpatialMicroTerrain) = size(t.dem)

"""
    elevation(terrain, i, j)

Get elevation at grid cell (i, j).
"""
elevation(t::SpatialMicroTerrain, i, j) = t.dem[i, j]
