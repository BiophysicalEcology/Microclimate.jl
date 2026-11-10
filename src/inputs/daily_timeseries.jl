"""
    DailyTimeseries(; shade, soil_wetness, surface_emissivity, cloud_emissivity,
                      rainfall, deep_soil_temperature, leaf_area_index, albedo, roughness_height)

One value per simulated day of the slower-varying forcings that aren't
resolved hourly. Each field is a vector of length `length(days)`.

- `shade` — fractional shade cast by vegetation (0–1)
- `soil_wetness` — fractional surface wetness (0–1)
- `surface_emissivity`, `cloud_emissivity` — longwave emissivities (0–1)
- `rainfall` — daily rainfall total (kg/m²)
- `deep_soil_temperature` — boundary condition at the deepest soil node (°C)
- `leaf_area_index` — leaf area index (m²/m²)
- `albedo` — ground albedo (0–1). Snow cover overrides it.
- `roughness_height` — ground roughness height (m)
"""
@kwdef struct DailyTimeseries{Sh,SW,SE,CE,R,DST,LAI,AB,RH} <: AbstractEnvironment
    shade::Sh
    soil_wetness::SW
    surface_emissivity::SE
    cloud_emissivity::CE
    rainfall::R
    deep_soil_temperature::DST
    leaf_area_index::LAI
    albedo::AB
    roughness_height::RH
end

"""
    example_daily_environment(days=DEFAULT_DAYS; kwargs...)

Example [`DailyTimeseries`](@ref) for Madison, Wisconsin, USA — one representative
day per month by default. `kwargs` override any field.
"""
function example_daily_environment(days=DEFAULT_DAYS;
    shade = fill(0.0, length(days)), # fractional shade cast by vegetation
    soil_wetness = fill(0.0, length(days)), # fractional surface wetness
    surface_emissivity = fill(0.96, length(days)), # - surface emissivity
    cloud_emissivity = fill(0.96, length(days)), # - cloud emissivity
    rainfall = ([28, 28.2, 54.6, 79.7, 81.3, 100.1, 101.3, 102.5, 89.7, 62.4, 54.9, 41.2])u"kg/m^2",
    deep_soil_temperature = fill(7.741666u"°C", length(days)),
    leaf_area_index = fill(0.1, length(days)),
    albedo = fill(0.15, length(days)),
    roughness_height = fill(0.004u"m", length(days)),
)
    DailyTimeseries(;
        shade, soil_wetness, surface_emissivity, cloud_emissivity,
        rainfall, deep_soil_temperature, leaf_area_index, albedo, roughness_height,
    )
end
