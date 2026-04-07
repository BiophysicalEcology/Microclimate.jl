"""
Spatial terrain for grid-based microclimate simulation.

Holds elevation (DEM) and precomputed terrain metrics from Geomorphometry.jl.
These are static (read-only) during simulation.
"""

using Geomorphometry: Geomorphometry
using Geomorphometry: slope as geo_slope, aspect as geo_aspect, filldepressions, depression_depth
using Unitful: Unitful, ustrip, @u_str

"""
    SpatialMicroTerrain

Static terrain data for spatial microclimate simulation.

# Fields
- `dem`: Digital elevation model / elevation grid (m)
- `slope`: Terrain slope (radians)
- `aspect`: Terrain aspect (radians from north)
- `sky_view_factor`: Fraction of sky visible (0-1, dimensionless)
- `filled_dem`: Depression-filled elevation for flow routing (m)
- `basin_depth`: Depth below pour point (m) - identifies pooling zones
- `cellsize`: Grid cell dimensions (dx, dy) with length units

# Constructor
    SpatialMicroTerrain(elevation; cellsize=(30.0u"m", 30.0u"m"))

Computes all derived terrain metrics from the elevation grid.
"""
struct SpatialMicroTerrain{E,S,A,V,F,B,C}
    dem::E              # elevation (m)
    slope::S            # radians
    aspect::A           # radians
    sky_view_factor::V  # dimensionless (0-1)
    filled_dem::F       # elevation (m)
    basin_depth::B      # depth (m)
    cellsize::C         # (dx, dy) with units
end

"""
    SpatialMicroTerrain(elevation; cellsize=(30.0u"m", 30.0u"m"))

Create terrain from an elevation grid with units.
Computes slope, aspect, sky view factor, and depression metrics.
"""
function SpatialMicroTerrain(
    elevation::AbstractMatrix;
    cellsize = (30.0u"m", 30.0u"m"),
)
    # Strip units for Geomorphometry.jl (it expects Float64)
    elev_stripped = ustrip.(u"m", elevation)
    cs_stripped = (ustrip(u"m", cellsize[1]), ustrip(u"m", cellsize[2]))

    # Compute terrain metrics using Geomorphometry.jl
    s = geo_slope(elev_stripped; cellsize=cs_stripped)
    a = geo_aspect(elev_stripped; cellsize=cs_stripped)
    # Sky view factor: simplified approximation from slope
    # svf ≈ cos(slope/2)² for unobstructed terrain
    svf = cos.(s ./ 2).^2
    filled = filldepressions(elev_stripped)
    bd = depression_depth(elev_stripped; filled=filled)

    # Add units back where appropriate
    SpatialMicroTerrain(
        elevation,           # keep original with units
        s,                   # radians (dimensionless in Geomorphometry)
        a,                   # radians
        svf,                 # dimensionless
        filled .* u"m",      # add units back
        bd .* u"m",          # add units back
        cellsize,            # keep units
    )
end

Base.size(t::SpatialMicroTerrain) = size(t.dem)

"""
    elevation(terrain, i, j)

Get elevation at grid cell (i, j).
"""
elevation(t::SpatialMicroTerrain, i, j) = t.dem[i, j]
