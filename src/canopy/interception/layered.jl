"""
    LayeredRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2")

Per-layer rain interception: extinction at a rain-fall angle
(`atan(wind_speed / raindrop_fall_velocity)`, drop velocity from uncited eq.
`3.78 * rainfall^0.067`, matching micropoint's `rainintercept`
exactly) using the same leaf-angle machinery as direct-beam shortwave
([`ellipsoidal_extinction_coefficient`](@ref)); `wind_speed` is supplied by
the caller, not read from a model here. Storage capacity scales with each
layer's own PAI (mass per unit leaf area, 0.1 kg/m² ~ 0.1 mm film); drip
(whatever exceeds capacity) is passed to the next layer as ordinary
incident rain, not a distinct (larger, near-vertical) drop stream. Wetness
([`wet_canopy_fraction`](@ref)) is a linear fractional-saturation blend.

Rain-fall angle assumes drops track local wind instantaneously -- a
trajectory approximation, not a physical drop model.
"""
@kwdef struct LayeredRainInterception{C} <: AbstractCanopyInterceptionModel
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
        @inbounds for i in 1:n_layers
            throughfall[i + 1] = throughfall[i]
        end
        return (; ground_throughfall = throughfall[n_layers + 1])
    end

    raindrop_fall_velocity = 3.78 * ustrip(u"kg/m^2", rainfall)^0.067 * u"m/s"
    capacity = model.leaf_water_storage_capacity
    @inbounds for i in 1:n_layers
        rain_zenith_angle = atan(ustrip(u"m/s", wind_speed[i]), ustrip(u"m/s", raindrop_fall_velocity)) * u"rad"
        _layer_rain_intercept!(leaf_surface_water, throughfall, i, rain_zenith_angle,
            canopy_projection_ratio, layer_plant_area_index[i], capacity)
    end

    return (; ground_throughfall = throughfall[n_layers + 1])
end
