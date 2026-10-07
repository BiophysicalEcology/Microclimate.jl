"""
    VerticalRainInterception(; leaf_water_storage_capacity=0.1u"kg/m^2",
                                rain_zenith_angle=0.0u"°")

Per-layer rain interception, same extinction/storage/drip bookkeeping as
[`LayeredRainInterception`](@ref), but at a fixed `rain_zenith_angle`
(default vertical) instead of a wind-driven one -- e.g. Shaw's `RAINC`,
which reuses the canopy's own direct-beam transmittance geometry at a
near-vertical, slope-adjusted angle rather than a terminal-velocity
relation. Set `rain_zenith_angle` to the site's `slope` for that
equivalence (rain falls vertically in the world frame, so its angle
relative to the canopy's local zenith equals terrain slope).
"""
@kwdef struct VerticalRainInterception{C,A} <: AbstractCanopyInterceptionModel
    leaf_water_storage_capacity::C = 0.1u"kg/m^2"
    rain_zenith_angle::A = 0.0u"°"
end

allocate_interception(::VerticalRainInterception, canopy_height, plant_area_index, heights, n_layers, boundary_layer_model) =
    allocate_interception(LayeredRainInterception(), canopy_height, plant_area_index, heights, n_layers, boundary_layer_model)

function canopy_interception!(buffers, model::VerticalRainInterception, canopy_projection_ratio;
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

    capacity = model.leaf_water_storage_capacity
    rain_zenith_angle = model.rain_zenith_angle
    @inbounds for i in 1:n_layers
        _layer_rain_intercept!(leaf_surface_water, throughfall, i, rain_zenith_angle,
            canopy_projection_ratio, layer_plant_area_index[i], capacity)
    end

    return (; ground_throughfall = throughfall[n_layers + 1])
end

# Storage bookkeeping shared with LayeredRainInterception -- identical for
# both, only the interception mechanism above differs.
const _StorageBasedInterception = Union{LayeredRainInterception,VerticalRainInterception}

function wet_canopy_fraction(model::_StorageBasedInterception, buffers, layer)
    layer_capacity = model.leaf_water_storage_capacity * buffers.layer_plant_area_index[layer]
    iszero(layer_capacity) && return 0.0
    return clamp(buffers.leaf_surface_water[layer] / layer_capacity, 0.0, 1.0)
end

available_canopy_water(::_StorageBasedInterception, buffers, layer) = buffers.leaf_surface_water[layer]

canopy_water_capacity(model::_StorageBasedInterception, buffers, layer) =
    model.leaf_water_storage_capacity * buffers.layer_plant_area_index[layer]

function deplete_canopy_water!(::_StorageBasedInterception, buffers, layer, evaporated_mass)
    buffers.leaf_surface_water[layer] = max(buffers.leaf_surface_water[layer] - evaporated_mass, 0.0u"kg/m^2")
    return nothing
end

leaf_standing_dew(::_StorageBasedInterception, buffers, layer) = buffers.leaf_standing_dew[layer]
leaf_standing_frost(::_StorageBasedInterception, buffers, layer) = buffers.leaf_standing_frost[layer]

function add_leaf_standing_dew!(::_StorageBasedInterception, buffers, layer, Δ)
    buffers.leaf_standing_dew[layer] = max(buffers.leaf_standing_dew[layer] + Δ, 0.0u"kg/m^2")
    return nothing
end

function add_leaf_standing_frost!(::_StorageBasedInterception, buffers, layer, Δ)
    buffers.leaf_standing_frost[layer] = max(buffers.leaf_standing_frost[layer] + Δ, 0.0u"kg/m^2")
    return nothing
end

function clamp_leaf_standing_dew!(::_StorageBasedInterception, buffers, layer, ceiling)
    buffers.leaf_standing_dew[layer] = min(buffers.leaf_standing_dew[layer], ceiling)
    return nothing
end

function snapshot_interception!(::_StorageBasedInterception, buffers)
    buffers.leaf_surface_water_day_start .= buffers.leaf_surface_water
    buffers.leaf_standing_dew_day_start .= buffers.leaf_standing_dew
    buffers.leaf_standing_frost_day_start .= buffers.leaf_standing_frost
    return nothing
end

function restore_interception!(::_StorageBasedInterception, buffers)
    buffers.leaf_surface_water .= buffers.leaf_surface_water_day_start
    buffers.leaf_standing_dew .= buffers.leaf_standing_dew_day_start
    buffers.leaf_standing_frost .= buffers.leaf_standing_frost_day_start
    return nothing
end

function reset_interception!(::_StorageBasedInterception, buffers)
    fill!(buffers.leaf_surface_water, 0.0u"kg/m^2")
    fill!(buffers.leaf_standing_dew, 0.0u"kg/m^2")
    fill!(buffers.leaf_standing_frost, 0.0u"kg/m^2")
    return nothing
end
