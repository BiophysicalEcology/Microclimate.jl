"""
    LayeredRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2")

Rain falls at angle `atan(wind_speed / drop_speed)`, with drop speed `3.78·rainfall^0.067`,
and each layer catches it like a direct beam ([`ellipsoidal_extinction_coefficient`](@ref)).
Each layer stores up to `leaf_water_storage_capacity` times its plant area index; the rest
drips to the layer below.
"""
@kwdef struct LayeredRainInterception{C} <: AbstractStorageInterception
    leaf_water_storage_capacity::C = 0.1u"kg/m^2"
end

function allocate_interception(::LayeredRainInterception, canopy_height, plant_area_index, heights, n_layers, boundary_layer_model)
    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    return (;
        layer_plant_area_index,
        leaf_surface_water = zeros(typeof(0.0u"kg/m^2"), n_layers),
        leaf_surface_water_day_start = zeros(typeof(0.0u"kg/m^2"), n_layers),
        leaf_standing_dew = zeros(typeof(0.0u"kg/m^2"), n_layers),
        leaf_standing_dew_day_start = zeros(typeof(0.0u"kg/m^2"), n_layers),
        leaf_standing_frost = zeros(typeof(0.0u"kg/m^2"), n_layers),
        leaf_standing_frost_day_start = zeros(typeof(0.0u"kg/m^2"), n_layers),
        throughfall = zeros(typeof(0.0u"kg/m^2"), n_layers + 1),
    )
end

function canopy_interception!(buffers, model::LayeredRainInterception, canopy_projection_ratio;
    rainfall, wind_speed,
)
    (; layer_plant_area_index, leaf_surface_water, throughfall) = buffers
    n_layers = length(leaf_surface_water)
    throughfall[1] = rainfall

    if rainfall <= 0.0u"kg/m^2"
        for i in 1:n_layers
            throughfall[i + 1] = throughfall[i]
        end
        return (; ground_throughfall = throughfall[n_layers + 1])
    end

    raindrop_fall_velocity = 3.78 * ustrip(u"kg/m^2", rainfall)^0.067 * u"m/s"
    capacity = model.leaf_water_storage_capacity
    for i in 1:n_layers
        rain_zenith_angle = atan(ustrip(u"m/s", wind_speed[i]), ustrip(u"m/s", raindrop_fall_velocity)) * u"rad"
        _layer_rain_intercept!(leaf_surface_water, throughfall, i, rain_zenith_angle,
            canopy_projection_ratio, layer_plant_area_index[i], capacity)
    end

    return (; ground_throughfall = throughfall[n_layers + 1])
end
