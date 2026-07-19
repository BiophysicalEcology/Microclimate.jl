"""
    LayeredRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2")

Per-layer rain interception: extinction at a rain-fall angle
(`atan(local_wind_speed / raindrop_fall_velocity)`, drop velocity from
Best 1950's `3.78 * rainfall^0.067`) using the same leaf-angle machinery as
direct-beam shortwave ([`ellipsoidal_extinction_coefficient`](@ref)), local
wind from [`CanopyWindAttenuation`](@ref). Storage capacity scales with each
layer's own PAI (mass per unit leaf area, 0.1 kg/m² ~ 0.1 mm film); drip is
whatever exceeds capacity each layer. Wetness ([`wet_canopy_fraction`](@ref))
is a linear fractional-saturation blend, not an empirical decay constant.

# References
- Best, A. C. (1950). The size distribution of raindrops. *Quarterly
  Journal of the Royal Meteorological Society*, 76(327), 16-36.
"""
@kwdef struct LayeredRainInterception{C} <: AbstractCanopyInterceptionModel
    leaf_water_storage_capacity::C = 0.1u"kg/m^2"
end

function allocate_interception(::LayeredRainInterception, canopy_height, plant_area_index, heights, n_layers, boundary_layer_model)
    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    return (;
        layer_plant_area_index,
        leaf_surface_water = zeros(typeof(0.0u"kg/m^2"), n_layers),
        throughfall = zeros(typeof(0.0u"kg/m^2"), n_layers + 1),
    )
end

function canopy_interception!(buffers, model::LayeredRainInterception, leaf_angle_distribution_parameter;
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
        rain_zenith_angle = atan(ustrip(u"m/s", wind_speed[i]) / ustrip(u"m/s", raindrop_fall_velocity)) * u"rad"
        extinction_coefficient = ellipsoidal_extinction_coefficient(rain_zenith_angle, leaf_angle_distribution_parameter)
        transmission = exp(-extinction_coefficient * layer_plant_area_index[i])

        intercepted = (1.0 - transmission) * throughfall[i]
        layer_capacity = capacity * layer_plant_area_index[i]
        new_storage = leaf_surface_water[i] + intercepted
        drip = max(new_storage - layer_capacity, 0.0u"kg/m^2")
        leaf_surface_water[i] = new_storage - drip
        throughfall[i + 1] = (throughfall[i] - intercepted) + drip
    end

    return (; ground_throughfall = throughfall[n_layers + 1])
end

function wet_canopy_fraction(model::LayeredRainInterception, buffers, layer)
    layer_capacity = model.leaf_water_storage_capacity * buffers.layer_plant_area_index[layer]
    iszero(layer_capacity) && return 0.0
    return clamp(buffers.leaf_surface_water[layer] / layer_capacity, 0.0, 1.0)
end

function deplete_canopy_water!(::LayeredRainInterception, buffers, layer, evaporated_mass)
    buffers.leaf_surface_water[layer] = max(buffers.leaf_surface_water[layer] - evaporated_mass, 0.0u"kg/m^2")
    return nothing
end

function reset_interception!(::LayeredRainInterception, buffers)
    fill!(buffers.leaf_surface_water, 0.0u"kg/m^2")
    return nothing
end
