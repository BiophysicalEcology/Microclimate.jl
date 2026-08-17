"""
    ExponentialCanopyWindAttenuation(; thermal_roughness_model=ScalarRoughnessRatio(; ratio=0.5),
                                       max_attenuation_coefficient=2.879)

In-canopy wind decays exponentially with cumulative plant area index from
canopy top: `wind_speed(cumulative_plant_area_index) = canopy_top_wind_speed *
exp(-attenuation_coefficient * cumulative_plant_area_index / total_plant_area_index)`.
`attenuation_coefficient = max_attenuation_coefficient * (1 - exp(-total_plant_area_index))`
saturates at `max_attenuation_coefficient` for a dense canopy and vanishes
for a sparse one -- unlike [`MixingLengthCanopyWindAttenuation`](@ref)'s shelter-factor
form, which over-attenuates wind at low total plant area index (see its own
docstring).

# References
- Nikolov, N., & Zeller, K. (2003). Modeling coupled interactions of carbon,
  water, and ozone exchange between terrestrial ecosystems and the
  atmosphere. I: Model description. *Environmental Pollution*, 124(2), 231-246.
"""
@kwdef struct ExponentialCanopyWindAttenuation{TRM,MAC} <: AbstractCanopyWindModel
    thermal_roughness_model::TRM = ScalarRoughnessRatio(; ratio=0.5)
    max_attenuation_coefficient::MAC = 2.879
end

# Wind fraction per layer (top-to-bottom, index 1 = canopy top): exponential
# decay with cumulative plant area index, bottom tenth replaced by a local
# log-law (see _apply_near_ground_loglaw!).
function exponential_wind_attenuation_profile(layer_plant_area_index, canopy_height, karman_constant, max_attenuation_coefficient)
    n = length(layer_plant_area_index)
    n >= 10 || throw(ArgumentError("exponential_wind_attenuation_profile requires at least 10 canopy layers"))
    total_plant_area_index = sum(layer_plant_area_index)
    iszero(total_plant_area_index) && return ones(n)

    attenuation_coefficient = max_attenuation_coefficient * (1.0 - exp(-total_plant_area_index))
    cumulative_plant_area_index_above = 0.0
    relative_wind_top_to_bottom = map(layer_plant_area_index) do plant_area_index_layer
        cumulative_plant_area_index_at_midpoint = cumulative_plant_area_index_above + plant_area_index_layer / 2
        cumulative_plant_area_index_above += plant_area_index_layer
        exp(-attenuation_coefficient * cumulative_plant_area_index_at_midpoint / total_plant_area_index)
    end

    relative_wind = reverse(relative_wind_top_to_bottom)  # bottom-to-top, matching _apply_near_ground_loglaw!'s convention
    _apply_near_ground_loglaw!(relative_wind, canopy_height, karman_constant, n ÷ 10)
    return reverse(relative_wind)
end

canopy_wind_profile!(buffers, ::ExponentialCanopyWindAttenuation, boundary_layer_model; kw...) =
    _shape_based_canopy_wind_profile!(buffers, boundary_layer_model; kw...)

function allocate_wind(model::ExponentialCanopyWindAttenuation, canopy_height, plant_area_index, boundary_layer_model, heights, n_layers)
    canopy_height = max(canopy_height, 1.0e-3u"m")  # avoid 0/0 in the roughness formula below
    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    total_plant_area_index = sum(plant_area_index)

    displacement_height = zero_plane_displacement(canopy_height, total_plant_area_index)
    roughness_length = canopy_roughness_length(
        canopy_height, total_plant_area_index, displacement_height, boundary_layer_model.karman_constant,
    )
    wind_attenuation = exponential_wind_attenuation_profile(
        layer_plant_area_index, canopy_height, boundary_layer_model.karman_constant, model.max_attenuation_coefficient,
    )

    (; layer_heights) = canopy_layer_heights(heights, canopy_height, n_layers)
    ground_layer_height = layer_heights[n_layers]

    return (;
        canopy_height, displacement_height, roughness_length, wind_attenuation, ground_layer_height,
        above_canopy_profile = allocate_profile(boundary_layer_model, [canopy_height, last(heights)]),
        wind_speed = zeros(typeof(0.0u"m/s"), n_layers),
        thermal_roughness_model = model.thermal_roughness_model,
    )
end
