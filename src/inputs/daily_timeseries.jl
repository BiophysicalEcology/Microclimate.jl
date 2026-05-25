@kwdef struct DailyTimeseries{Sh,SW,SE,CE,R,DST,LAI,AB} <: AbstractEnvironment
    shade::Sh
    soil_wetness::SW
    surface_emissivity::SE
    cloud_emissivity::CE
    rainfall::R
    deep_soil_temperature::DST
    leaf_area_index::LAI
    # Daily surface albedo. `nothing` (default) → fall back to the scalar
    # `Site.albedo`. A vector of length `length(days)` supplies time-varying
    # albedo (e.g. derived from a MODIS BRDF/Albedo product). Snow-cover
    # overrides still take precedence when snow is active.
    albedo::AB = nothing
end

function example_daily_environment(days=DEFAULT_DAYS;
    shade = fill(0.0, length(days)), # fractional shade cast by vegetation
    soil_wetness = fill(0.0, length(days)), # fractional surface wetness
    surface_emissivity = fill(0.96, length(days)), # - surface emissivity
    cloud_emissivity = fill(0.96, length(days)), # - cloud emissivity
    rainfall = ([28, 28.2, 54.6, 79.7, 81.3, 100.1, 101.3, 102.5, 89.7, 62.4, 54.9, 41.2])u"kg/m^2",
    deep_soil_temperature = fill(7.741666u"°C", length(days)),
    leaf_area_index = fill(0.1, length(days)),
    albedo = nothing,
)
    DailyTimeseries(;
        shade, soil_wetness, surface_emissivity, cloud_emissivity,
        rainfall, deep_soil_temperature, leaf_area_index, albedo,
    )
end
