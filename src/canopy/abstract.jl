abstract type AbstractCanopyModel end

"""
    allocate_canopy(model, heights, boundary_layer_model)

Pre-allocate the per-layer buffers a canopy model needs to evaluate its
physics each hour, sized to the sub-canopy portion of `heights` (the same
`MicroModel.heights` vector the boundary-layer model already uses for its
above-canopy profile). Structural quantities that depend only on `model`
and `boundary_layer_model` (radiative-transfer optics, wind-attenuation
shape, and similar canopy-structure constants) are computed once here, not
per hour.
"""
function allocate_canopy end

"""
    canopy_shortwave!(buffers, model; zenith_angle, direct_horizontal_irradiance,
                       diffuse_horizontal_irradiance, ground_reflectance)

Compute layer-resolved absorbed shortwave radiation through the canopy for
the current hour, writing per-layer results into `buffers` and returning
`(; ground_absorbed_shortwave, canopy_absorbed_shortwave)`. `ground_reflectance`
is supplied by the caller (e.g. from `Site.albedo` or the current
`environment_instant`) rather than stored on the canopy model, matching how
other radiation models here read albedo/emissivity from their caller.
"""
function canopy_shortwave! end

"""
    canopy_wind_profile!(buffers, model, boundary_layer_model; site, environment_instant,
                          canopy_source_temperature, vapour_pressure_equation)

Compute the per-layer wind-speed profile through the canopy for the current
hour, writing into `buffers` and returning `(; canopy_top_wind_speed,
canopy_top_air_temperature, friction_velocity)`. Above and at the canopy
top, this is `atmospheric_surface_profile!` evaluated with the canopy's
displacement height and roughness length; below canopy top, the structural
wind-attenuation shape scales that boundary value down through the layers.
"""
function canopy_wind_profile! end

"""
    canopy_longwave!(buffers, model; leaf_temperature, ground_temperature, ground_emissivity,
                      site, environment_instant, vapour_pressure_equation)

Compute layer-resolved longwave exchange (leaf emission/absorption plus sky
and ground boundary terms) for the current hour, writing per-layer results
into `buffers` and returning `(; ground_absorbed_longwave,
canopy_absorbed_longwave)`. `leaf_temperature` is the current Picard-pass
guess (or previous hour's converged value); `ground_temperature`/
`ground_emissivity` are the *lagged* soil-surface boundary condition (last
hour's converged soil state), matching how [`canopy_radiation!`](@ref)-style
models read boundary data from the caller rather than storing it.
"""
function canopy_longwave! end

"""
    canopy_layer_geometry(plant_area_index, n_layers)

Per-layer plant-area-index geometry shared by every canopy sub-model that
discretises the canopy into `n_layers` equal-thickness slabs (shortwave,
longwave, wind attenuation): `layer_thickness` (PAI per layer),
`layer_plant_area_index` (one value per layer, all equal to
`layer_thickness`), the cumulative PAI at each layer *boundary*
(`n_layers + 1` values, boundary 1 = canopy top), and the cumulative PAI
*above the midpoint* of each layer (`n_layers` values).
"""
function canopy_layer_geometry(plant_area_index, n_layers)
    layer_thickness = plant_area_index / n_layers
    layer_plant_area_index = fill(layer_thickness, n_layers)
    layer_plant_area_index_boundaries = [(i - 1) * layer_thickness for i in 1:(n_layers + 1)]
    layer_plant_area_index_above = [(i - 0.5) * layer_thickness for i in 1:n_layers]
    return (; layer_thickness, layer_plant_area_index, layer_plant_area_index_boundaries, layer_plant_area_index_above)
end
