@kwdef struct DailyTimeseries{Sh,SW,SE,CE,R,DST,LAI} <: AbstractEnvironment
    shade::Sh
    soil_wetness::SW
    surface_emissivity::SE
    cloud_emissivity::CE
    rainfall::R
    deep_soil_temperature::DST
    leaf_area_index::LAI
end
