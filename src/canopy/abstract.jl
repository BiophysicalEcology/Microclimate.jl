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
    canopy_shortwave!(buffers, model; zenith_angle, direct_horizontal_irradiance,
                       diffuse_horizontal_irradiance, ground_reflectance)

Layer-resolved absorbed shortwave through the canopy for the current hour;
writes per-layer results into `buffers`, returns
`(; ground_absorbed_shortwave, canopy_absorbed_shortwave)`.
`ground_reflectance` is supplied by the caller, not stored on the model.
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
`buffers`, returns `(; ground_absorbed_longwave, canopy_absorbed_longwave)`.
`ground_temperature`/`ground_emissivity` are the lagged (previous-hour)
soil-surface boundary condition.
"""
function canopy_longwave! end

"""
    canopy_layer_geometry(plant_area_index, n_layers)

Per-layer plant-area-index geometry shared by every canopy sub-model that
discretises the canopy into `n_layers` equal-thickness slabs:
`layer_thickness` (PAI per layer), `layer_plant_area_index` (one value per
layer), cumulative PAI at each layer boundary (`n_layers+1` values, boundary
1 = canopy top), and cumulative PAI above each layer's midpoint.
"""
function canopy_layer_geometry(plant_area_index, n_layers)
    layer_thickness = plant_area_index / n_layers
    layer_plant_area_index = fill(layer_thickness, n_layers)
    layer_plant_area_index_boundaries = [(i - 1) * layer_thickness for i in 1:(n_layers + 1)]
    layer_plant_area_index_above = [(i - 0.5) * layer_thickness for i in 1:n_layers]
    return (; layer_thickness, layer_plant_area_index, layer_plant_area_index_boundaries, layer_plant_area_index_above)
end
