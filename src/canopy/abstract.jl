abstract type AbstractCanopyModel end

"""
    allocate_canopy(model, heights, boundary_layer_model)

Pre-allocate a canopy model's per-layer buffers, sized to the sub-canopy
portion of `heights`. Structural quantities depending only on `model`/
`boundary_layer_model` are computed once here, not per hour.
"""
function allocate_canopy end

"""
    reset_canopy_scratch!(model, buffers)

Zero day-boundary-reset canopy state (e.g. interception storage). No-op for
`NoCanopy`.
"""
function reset_canopy_scratch! end

"""
    allocate_canopy_inputs(model; site, environment_instant, boundary_layer_model)

Build a canopy model's mutable per-hour inputs (`nothing` for `NoCanopy`;
a `CanopyEnergyBalanceInputs` for `MultilayerCanopy`).
"""
function allocate_canopy_inputs end

"""
    initial_ground_overrides(model)

Properly-typed placeholder `(; ground_shortwave_transmission,
ground_incoming_longwave)` for `SoilEnergyInputs`'s initial construction —
`nothing` for `NoCanopy`, concrete for `MultilayerCanopy`. Fixes the field's
type up front so later per-hour updates don't hit a `Nothing`-locked field.
"""
function initial_ground_overrides end

"""
    canopy_leaf_area_index(model)

Leaf area index for Campbell soil-hydraulics' Beer's-law soil/canopy
potential-ET split. `nothing` for `NoCanopy` (falls back to
`environment_instant.leaf_area_index`, the weather-driven daily series);
`sum(plant_area_index) * (1 - woody_area_fraction)` for `MultilayerCanopy`,
narrowing total PAI down to living-leaf LAI.
"""
function canopy_leaf_area_index end

"""
    canopy_shortwave!(buffers, model; zenith_angle, direct_horizontal_irradiance,
                       diffuse_horizontal_irradiance, ground_reflectance)

Layer-resolved absorbed shortwave through the canopy for the current hour;
writes per-layer results into `buffers`, returns `(;
ground_absorbed_shortwave, landscape_absorbed_shortwave)`. `ground_reflectance`
is supplied by the caller, not stored on the model. Despite the name,
`landscape_absorbed_shortwave` is net radiation retained by the whole
canopy+ground system below the top boundary (`incoming -
top_of_canopy_reflected`), not canopy-only absorption — for that, sum the
per-layer `absorbed_shortwave` buffer instead.
"""
function canopy_shortwave! end

"""
    canopy_wind_profile!(buffers, model, boundary_layer_model; site, environment_instant,
                          canopy_source_temperature, vapour_pressure_equation)

Per-layer wind-speed profile through the canopy for the current hour;
writes into `buffers`, returns `(; canopy_top_wind_speed,
canopy_top_air_temperature, friction_velocity)`. Above canopy top this is
`atmospheric_surface_profile!` with the canopy's displacement height and
roughness length; below top, the wind-attenuation shape scales that
boundary value down through the layers.
"""
function canopy_wind_profile! end

"""
    canopy_longwave!(buffers, model; leaf_temperature, ground_temperature, ground_emissivity,
                      site, environment_instant, vapour_pressure_equation)

Layer-resolved longwave exchange (leaf emission/absorption plus sky/ground
boundary terms) for the current hour; writes per-layer results into
`buffers`, returns `(; ground_absorbed_longwave, landscape_absorbed_longwave)`.
`ground_temperature`/`ground_emissivity` are the lagged (previous-hour)
soil-surface boundary condition. As with [`canopy_shortwave!`](@ref),
`landscape_absorbed_longwave` is net longwave retained by the whole
canopy+ground system below the top boundary (`incoming_longwave -
top_of_canopy_upward`), not canopy-only absorption — sum the per-layer
`absorbed_longwave` buffer for that.
"""
function canopy_longwave! end

"""
    canopy_layer_geometry(plant_area_index, n_layers)

Per-layer PAI geometry shared by every canopy sub-model: `layer_plant_area_index`
(one value per layer), cumulative PAI at each layer boundary (`n_layers+1`
values, boundary 1 = canopy top), and cumulative PAI above each layer's
midpoint. `plant_area_index` is a scalar total (split into `n_layers`
equal-PAI slabs) or an `AbstractVector` of `n_layers` per-layer values,
top-to-bottom.
"""
function canopy_layer_geometry(plant_area_index::Real, n_layers)
    return canopy_layer_geometry(fill(plant_area_index / n_layers, n_layers), n_layers)
end
function canopy_layer_geometry(layer_plant_area_index::AbstractVector, n_layers)
    length(layer_plant_area_index) == n_layers ||
        throw(ArgumentError("plant_area_index vector must have one entry per canopy layer"))
    layer_plant_area_index_boundaries = similar(layer_plant_area_index, n_layers + 1)
    layer_plant_area_index_boundaries[1] = zero(eltype(layer_plant_area_index))
    @inbounds for i in 1:n_layers
        layer_plant_area_index_boundaries[i + 1] = layer_plant_area_index_boundaries[i] + layer_plant_area_index[i]
    end
    layer_plant_area_index_above = [layer_plant_area_index_boundaries[i] + layer_plant_area_index[i] / 2 for i in 1:n_layers]
    return (; layer_plant_area_index, layer_plant_area_index_boundaries, layer_plant_area_index_above)
end

"""
    canopy_layer_heights(heights, canopy_height, n_layers; min_spacing=1.0e-3u"m")

Canopy layer node heights (top-to-bottom) and each layer's height thickness
(distance to the next node, or to the ground for the bottom layer), floored
at `min_spacing` (a node can legitimately sit exactly at `canopy_height`).
"""
function canopy_layer_heights(heights, canopy_height, n_layers; min_spacing=1.0e-3u"m")
    layer_heights = sort(heights[heights .<= canopy_height]; rev=true)
    length(layer_heights) == n_layers ||
        throw(ArgumentError("number of `heights` at or below `canopy_height` must equal n_layers"))

    layer_thickness = similar(layer_heights)
    @inbounds for i in 1:(n_layers - 1)
        layer_thickness[i] = max(layer_heights[i] - layer_heights[i + 1], min_spacing)
    end
    layer_thickness[n_layers] = max(layer_heights[n_layers], min_spacing)  # bottom layer down to the ground (z=0)
    return (; layer_heights, layer_thickness)
end

"""
    plant_area_index_from_density(density, heights, canopy_height)

Convert a leaf/plant area density profile (`density`, PAI per unit height,
one value per canopy layer, top-to-bottom) into the absolute per-layer PAI
vector `MultilayerCanopy` expects, weighting each layer by its own height
thickness. Call once at model-construction time, not per solve.
"""
function plant_area_index_from_density(density, heights, canopy_height)
    n_layers = length(density)
    (; layer_thickness) = canopy_layer_heights(heights, canopy_height, n_layers)
    return density .* layer_thickness
end
