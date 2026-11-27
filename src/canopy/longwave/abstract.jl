"""
    AbstractCanopyLongwaveModel

Longwave exchange through canopy layers: [`LayeredLongwaveExchange`](@ref) or
[`LayeredRadiosityExchange`](@ref).
"""
abstract type AbstractCanopyLongwaveModel end

"""
    allocate_longwave(longwave_model, plant_area_index, n_layers, canopy_projection_ratio)

Allocate a model's buffers. Models that depend on leaf angle use `canopy_projection_ratio`.
"""
function allocate_longwave end

canopy_incoming_longwave(model::AbstractCanopyLongwaveModel; site, environment_instant, vapour_pressure_equation=GoffGratch()) =
    precompute_longwave_sky(model.atmospheric_radiation_model;
        site, environment_instant, vapour_pressure_equation, shade=0.0).incoming_longwave
