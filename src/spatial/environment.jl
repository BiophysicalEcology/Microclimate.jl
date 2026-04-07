"""
Unified environment type for spatial microclimate simulation.

Consolidates all external inputs:
- Weather (temporal): temperature, wind, humidity, radiation, precipitation
- Terrain (spatial): DEM and derived metrics (in SpatialMicroTerrain)
- Land surface (spatial): roughness, albedo, friction coefficients
- Soil (spatial): hydraulic properties

All components are parametric to support Unitful types.
"""

# =============================================================================
# Weather (Temporal Forcing)
# =============================================================================

"""
    SpatialWeather

Weather forcing for spatial simulation.

Can be spatially uniform (scalars) or spatially varying (grids).
Temporal variation is handled at the solver level.

# Fields
- `air_temperature`: Air temperature (K or °C)
- `wind_speed`: Wind speed (m/s)
- `wind_direction`: Wind direction (radians from north, optional)
- `humidity`: Relative humidity (0-1) or specific humidity
- `pressure`: Atmospheric pressure (Pa)
- `global_radiation`: Incoming shortwave radiation (W/m²)
- `rainfall`: Precipitation rate (m/hour or kg/m²/hour)
- `cloud_cover`: Cloud fraction (0-1)
"""
struct SpatialWeather{Ta,W,Wd,H,P,R,Rf,C}
    air_temperature::Ta
    wind_speed::W
    wind_direction::Wd  # optional, for KLAM21 regional wind term
    humidity::H
    pressure::P
    global_radiation::R
    rainfall::Rf
    cloud_cover::C
end

function SpatialWeather(;
    air_temperature,
    wind_speed,
    wind_direction = nothing,
    humidity = 0.5,
    pressure = 101325.0,
    global_radiation = 0.0,
    rainfall = 0.0,
    cloud_cover = 0.0,
)
    SpatialWeather(
        air_temperature,
        wind_speed,
        wind_direction,
        humidity,
        pressure,
        global_radiation,
        rainfall,
        cloud_cover,
    )
end

# =============================================================================
# Land Surface Properties (Spatial)
# =============================================================================

"""
    LandSurface

Land surface properties that affect energy and water balance.

All fields can be scalar (uniform) or matrix (spatially varying).

# Fields
- `roughness_length`: Aerodynamic roughness length z₀ (m) - for KLAM21 friction
- `albedo`: Surface albedo (0-1)
- `emissivity`: Surface emissivity (0-1)
- `surface_water_friction`: Friction coefficient for surface water flow (0-1)
"""
struct LandSurface{Z,A,E,F}
    roughness_length::Z
    albedo::A
    emissivity::E
    surface_water_friction::F
end

function LandSurface(;
    roughness_length = 0.05,  # m, typical grass
    albedo = 0.2,
    emissivity = 0.96,
    surface_water_friction = 0.1,
)
    LandSurface(roughness_length, albedo, emissivity, surface_water_friction)
end

# =============================================================================
# Soil Properties (Spatial)
# =============================================================================

"""
    SpatialSoilProperties

Soil hydraulic and thermal properties for spatial simulation.

All fields can be scalar (uniform) or matrix (spatially varying).

# Fields
- `saturated_hydraulic_conductivity`: K_sat (m/hour)
- `saturation_moisture`: θ_sat (m³/m³)
- `infiltration_rate`: Base infiltration rate (m/hour)
"""
struct SpatialSoilProperties{K,Θ,I}
    saturated_hydraulic_conductivity::K
    saturation_moisture::Θ
    infiltration_rate::I
end

function SpatialSoilProperties(;
    saturated_hydraulic_conductivity = 0.001,  # m/hour
    saturation_moisture = 0.45,
    infiltration_rate = 0.001,
)
    SpatialSoilProperties(
        saturated_hydraulic_conductivity,
        saturation_moisture,
        infiltration_rate,
    )
end

# =============================================================================
# Unified Environment
# =============================================================================

"""
    SpatialEnvironment

Unified environment for spatial microclimate simulation.

Consolidates all external inputs into a single type with nested components.

# Fields
- `terrain`: SpatialMicroTerrain with DEM and derived metrics
- `weather`: SpatialWeather with temporal forcing
- `land_surface`: LandSurface with surface properties
- `soil`: SpatialSoilProperties with soil hydraulic properties
"""
struct SpatialEnvironment{T<:SpatialMicroTerrain,W<:SpatialWeather,L<:LandSurface,S<:SpatialSoilProperties}
    terrain::T
    weather::W
    land_surface::L
    soil::S
end

function SpatialEnvironment(
    terrain::SpatialMicroTerrain;
    weather = nothing,
    land_surface = LandSurface(),
    soil = SpatialSoilProperties(),
)
    if isnothing(weather)
        error("weather must be provided")
    end
    SpatialEnvironment(terrain, weather, land_surface, soil)
end

# =============================================================================
# Accessor functions for convenient nested access
# =============================================================================

# Terrain accessors (delegate to SpatialMicroTerrain)
dem(env::SpatialEnvironment) = env.terrain.dem
slope(env::SpatialEnvironment) = env.terrain.slope
aspect(env::SpatialEnvironment) = env.terrain.aspect
sky_view_factor(env::SpatialEnvironment) = env.terrain.sky_view_factor
cellsize(env::SpatialEnvironment) = env.terrain.cellsize

# Weather accessors
air_temperature(env::SpatialEnvironment) = env.weather.air_temperature
wind_speed(env::SpatialEnvironment) = env.weather.wind_speed
wind_direction(env::SpatialEnvironment) = env.weather.wind_direction
humidity(env::SpatialEnvironment) = env.weather.humidity
rainfall(env::SpatialEnvironment) = env.weather.rainfall

# Land surface accessors
roughness_length(env::SpatialEnvironment) = env.land_surface.roughness_length
albedo(env::SpatialEnvironment) = env.land_surface.albedo
emissivity(env::SpatialEnvironment) = env.land_surface.emissivity
surface_water_friction(env::SpatialEnvironment) = env.land_surface.surface_water_friction

# Soil accessors
infiltration_rate(env::SpatialEnvironment) = env.soil.infiltration_rate
saturation_moisture(env::SpatialEnvironment) = env.soil.saturation_moisture

# =============================================================================
# Grid/scalar value access helpers
# =============================================================================

"""
    cell_value(x, i, j)

Get value at grid cell (i, j), handling both scalars and matrices.

For scalars (including Unitful quantities), returns the value unchanged.
For matrices, returns `x[i, j]`.
"""
cell_value(x::AbstractMatrix, i, j) = x[i, j]
cell_value(::Nothing, i, j) = nothing
cell_value(x, i, j) = x  # scalar fallback (works for Real, Quantity, etc.)

# =============================================================================
# Derived quantities
# =============================================================================

"""
    vapor_pressure(temperature, relative_humidity)

Compute vapor pressure from temperature and relative humidity.

Uses the Magnus formula for saturation vapor pressure.

# Arguments
- `temperature`: Air temperature (K)
- `relative_humidity`: Relative humidity (0-1)

# Returns
Vapor pressure (Pa)
"""
function vapor_pressure(temperature, relative_humidity)
    # Convert to Celsius for Magnus formula
    T_C = ustrip(u"K", temperature) - 273.15
    # Magnus formula: e_sat = 611.2 * exp(17.67 * T / (T + 243.5))
    e_sat = 611.2 * exp(17.67 * T_C / (T_C + 243.5))  # Pa
    return relative_humidity * e_sat * u"Pa"
end

"""
    vapor_pressure(env::SpatialEnvironment)
    vapor_pressure(env::SpatialEnvironment, i, j)

Compute vapor pressure from environment's temperature and humidity.
"""
function vapor_pressure(env::SpatialEnvironment)
    T = env.weather.air_temperature
    rh = env.weather.humidity
    if T isa AbstractMatrix
        return vapor_pressure.(T, rh isa AbstractMatrix ? rh : Ref(rh))
    else
        return vapor_pressure(T, rh isa AbstractMatrix ? rh[1,1] : rh)
    end
end

function vapor_pressure(env::SpatialEnvironment, i, j)
    T = get_value(env.weather.air_temperature, i, j)
    rh = get_value(env.weather.humidity, i, j)
    return vapor_pressure(T, rh)
end

"""
    regional_wind(env::SpatialEnvironment)

Get regional wind as (u, v) velocity components from speed and direction.

Returns `nothing` if wind_direction is not specified.
"""
function regional_wind(env::SpatialEnvironment)
    ws = env.weather.wind_speed
    wd = env.weather.wind_direction
    isnothing(wd) && return nothing

    # Convert direction (from north, clockwise) to velocity components
    # u = east-west component, v = north-south component
    u = ws * sin(wd)
    v = ws * cos(wd)
    return (u, v)
end
