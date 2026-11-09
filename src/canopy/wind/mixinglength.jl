"""
    MixingLengthCanopyWindAttenuation(; thermal_roughness_model=ScalarRoughnessRatio(),
                            shelter_floor=0.003, shelter_plant_area_coefficient=0.1,
                            mixing_length_coefficient=2.0, mixing_length_plant_area_coefficient=0.25)

Below the canopy top, wind decays by a mixing-length shape; at and above it, the
boundary-layer model runs with the canopy's displacement height and roughness length.
`thermal_roughness_model` replaces the boundary-layer model's for the canopy top only.
The other keywords are empirical shelter and mixing-length constants. May
over-attenuate wind at low plant area index.
"""
@kwdef struct MixingLengthCanopyWindAttenuation{TRM,SF,SPC,MLC,MLPC} <: AbstractCanopyWindModel
    thermal_roughness_model::TRM = ScalarRoughnessRatio()
    shelter_floor::SF = 0.003
    shelter_plant_area_coefficient::SPC = 0.1
    mixing_length_coefficient::MLC = 2.0
    mixing_length_plant_area_coefficient::MLPC = 0.25
end

"""
    zero_plane_displacement(canopy_height, plant_area_index)

Zero-plane displacement height from canopy height and total plant area
index (structural, computed once).

# References
- Raupach, M. R. (1994). Simplified expressions for vegetation roughness
  length and zero-plane displacement as functions of canopy height and
  area index. *Boundary-Layer Meteorology*, 71(1-2), 211-216.
"""
function zero_plane_displacement(canopy_height, plant_area_index)
    plant_area_index = max(plant_area_index, 0.001)  # numerical floor, not a physical parameter
    x = sqrt(7.5 * plant_area_index)
    return (1.0 - (1.0 - exp(-x)) / x) * canopy_height
end

"""
    canopy_roughness_length(canopy_height, plant_area_index, displacement_height, karman_constant;
                             shelter_floor=0.003, shelter_plant_area_coefficient=0.1)

Neutral roughness length from canopy structure.
"""
function canopy_roughness_length(canopy_height, plant_area_index, displacement_height, karman_constant;
    shelter_floor=0.003, shelter_plant_area_coefficient=0.1,
)
    canopy_element_shelter = sqrt(shelter_floor + shelter_plant_area_coefficient * plant_area_index)
    canopy_depth = canopy_height - displacement_height
    roughness = canopy_depth * exp(-karman_constant / canopy_element_shelter)
    return clamp(roughness, 0.0005u"m", 0.9 * canopy_depth)  # free/tunable floor and ceiling, uncited
end

# Attenuation coefficient for one plant-area-index element (structural
# mixing-length argument): larger local plant area density -> more
# attenuation per unit height.
function _element_attenuation(element_plant_area_index, canopy_height,
    shelter_floor, shelter_plant_area_coefficient, mixing_length_coefficient, mixing_length_plant_area_coefficient,
)
    element_shelter = sqrt(shelter_floor + shelter_plant_area_coefficient * element_plant_area_index)
    area_density = element_plant_area_index / canopy_height
    mixing_length = mixing_length_coefficient * element_shelter^3 / (mixing_length_plant_area_coefficient * area_density)
    return element_shelter * canopy_height / mixing_length
end

"""
    wind_attenuation_profile(layer_plant_area_index, canopy_height, plant_area_index_total, karman_constant;
                              shelter_floor=0.003, shelter_plant_area_coefficient=0.1,
                              mixing_length_coefficient=2.0, mixing_length_plant_area_coefficient=0.25)

Each layer's wind speed as a fraction of the canopy-top wind, top-to-bottom. Needs at
least 10 layers; the bottom tenth uses a log law so wind doesn't reach zero.
"""
function wind_attenuation_profile(layer_plant_area_index, canopy_height, plant_area_index_total, karman_constant;
    shelter_floor=0.003, shelter_plant_area_coefficient=0.1, mixing_length_coefficient=2.0, mixing_length_plant_area_coefficient=0.25,
)
    n = length(layer_plant_area_index)
    n >= 10 || throw(ArgumentError("wind_attenuation_profile requires at least 10 canopy layers"))

    # internally bottom-to-top (index n = canopy top = unattenuated)
    whole_canopy_attenuation = _element_attenuation(plant_area_index_total, canopy_height,
        shelter_floor, shelter_plant_area_coefficient, mixing_length_coefficient, mixing_length_plant_area_coefficient)
    layer_attenuation = [_element_attenuation(element_plant_area_index, canopy_height,
        shelter_floor, shelter_plant_area_coefficient, mixing_length_coefficient, mixing_length_plant_area_coefficient)
        for element_plant_area_index in Iterators.reverse(layer_plant_area_index)]
    layer_attenuation_sum = sum(layer_attenuation)
    # zero total PAI (leafless canopy): every element is already 0, no attenuation to rescale
    if !iszero(layer_attenuation_sum)
        layer_attenuation .*= whole_canopy_attenuation / layer_attenuation_sum
    end

    n_ground_layers = n ÷ 10
    relative_wind = ones(n)
    for i in n:-1:(n_ground_layers + 1)
        relative_wind[i - 1] = relative_wind[i] * (1.0 - layer_attenuation[i - 1])
    end
    _apply_near_ground_loglaw!(relative_wind, canopy_height, karman_constant, n_ground_layers)
    return reverse(relative_wind)
end

# Replaces the bottom n_ground_layers of a bottom-to-top relative-wind array
# (index n = canopy top) with a local log-law, calibrated to match whatever
# shape produced relative_wind[n_ground_layers] -- avoids wind staying
# unrealistically high right at the ground (a plain decay-with-height/PAI
# shape need not reach zero there, especially for a sparse canopy).
function _apply_near_ground_loglaw!(relative_wind, canopy_height, karman_constant, n_ground_layers)
    # 10.0 here is the same "bottom tenth" ratio as n_ground_layers = n÷10;
    # 20.0 sets ground_roughness ~ half the bottom layer's own thickness.
    κ = karman_constant
    ground_roughness = canopy_height / (20.0 * n_ground_layers)
    friction_scale = κ * relative_wind[n_ground_layers] / log(canopy_height / (10.0 * ground_roughness))
    for i in 1:n_ground_layers
        z = i * canopy_height / (10.0 * n_ground_layers)
        relative_wind[i] = (friction_scale / κ) * log(z / ground_roughness)
    end
    return relative_wind
end

function allocate_wind(model::MixingLengthCanopyWindAttenuation, canopy_height, plant_area_index, boundary_layer_model, heights, n_layers)
    canopy_height = max(canopy_height, 1.0e-3u"m")  # avoid 0/0 in the per-element/roughness formulas below
    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    total_plant_area_index = sum(plant_area_index)

    (; shelter_floor, shelter_plant_area_coefficient, mixing_length_coefficient, mixing_length_plant_area_coefficient) = model

    displacement_height = zero_plane_displacement(canopy_height, total_plant_area_index)
    roughness_length = canopy_roughness_length(
        canopy_height, total_plant_area_index, displacement_height, boundary_layer_model.karman_constant;
        shelter_floor, shelter_plant_area_coefficient,
    )
    wind_attenuation = wind_attenuation_profile(
        layer_plant_area_index, canopy_height, total_plant_area_index, boundary_layer_model.karman_constant;
        shelter_floor, shelter_plant_area_coefficient, mixing_length_coefficient, mixing_length_plant_area_coefficient,
    )
    # Ground-most layer's height above the soil surface — not every
    # air_profile_model buffer stores this (KTheoryAirProfile doesn't), so
    # compute it once here where it's cheap (allocation-time only).
    (; layer_heights) = canopy_layer_heights(heights, canopy_height, n_layers)
    ground_layer_height = layer_heights[n_layers]

    return (;
        canopy_height, displacement_height, roughness_length, wind_attenuation, ground_layer_height,
        above_canopy_profile = allocate_profile(boundary_layer_model, [canopy_height, last(heights)]),
        wind_speed = zeros(typeof(0.0u"m/s"), n_layers),
        thermal_roughness_model = model.thermal_roughness_model,
    )
end

canopy_wind_profile!(buffers, ::MixingLengthCanopyWindAttenuation, boundary_layer_model; kw...) =
    _shape_based_canopy_wind_profile!(buffers, boundary_layer_model; kw...)

# Shared by every AbstractCanopyWindModel whose allocate_wind produces a
# precomputed wind_attenuation shape array (MixingLengthCanopyWindAttenuation,
# ExponentialCanopyWindAttenuation) -- only the shape differs between them.
function _shape_based_canopy_wind_profile!(buffers, boundary_layer_model;
    site, environment_instant, canopy_source_temperature, vapour_pressure_equation=GoffGratch(),
    inverse_obukhov_length=nothing,
)
    profile = atmospheric_surface_profile!(boundary_layer_model, buffers.above_canopy_profile;
        site, environment_instant, surface_temperature=canopy_source_temperature, vapour_pressure_equation,
        displacement_height=buffers.displacement_height, roughness_length=buffers.roughness_length,
        thermal_roughness_model=buffers.thermal_roughness_model,
        temperature_anchor_height=buffers.canopy_height,
        stability_variable=InverseObukhovLength(),
        stability=inverse_obukhov_length,
    )
    canopy_top_wind_speed = first(profile.wind_speed)
    canopy_top_air_temperature = first(profile.air_temperature)
    canopy_top_relative_humidity = first(profile.relative_humidity)

    wind_speed = buffers.wind_speed
    for i in eachindex(wind_speed)
        wind_speed[i] = buffers.wind_attenuation[i] * canopy_top_wind_speed
    end

    return (; canopy_top_wind_speed, canopy_top_air_temperature, canopy_top_relative_humidity,
        friction_velocity=profile.friction_velocity,
        inverse_obukhov_length=profile.inverse_obukhov_length)
end

"""
    canopy_top_flux_boundary(boundary_layer_model, canopy_height, displacement_height, reference_height,
                              friction_velocity, inverse_obukhov_length, sensible_heat_flux, vapor_flux,
                              air_volumetric_heat_capacity, reference_temperature, reference_vapor_density,
                              previous_top_temperature, previous_top_vapor_density; max_resistance=2000.0u"s/m")

The air temperature and vapor density at the top of the canopy.

`sensible_heat_flux` and `vapor_flux` are what the leaves and the ground send up through the
canopy top. The canopy-top value is the one that carries them on to the reference height.

`max_resistance` caps the resistance in calm air. The result moves from the previous value
towards the balance only as far as the air could in an hour.

Returns `(; top_temperature, top_vapor_density, sensible_flux)`, where `sensible_flux` is
the upward heat flux from the canopy top to the reference height.
"""
function canopy_top_flux_boundary(boundary_layer_model, canopy_height, displacement_height, reference_height,
    friction_velocity, inverse_obukhov_length, sensible_heat_flux, vapor_flux, air_volumetric_heat_capacity,
    reference_temperature, reference_vapor_density, previous_top_temperature, previous_top_vapor_density;
    max_resistance=2000.0u"s/m",
)
    (; karman_constant, stability_functions) = boundary_layer_model
    top_depth = max(canopy_height - displacement_height, 1.0e-3u"m")
    reference_depth = max(reference_height - displacement_height, 1.0e-3u"m")
    heat_correction(depth) = heat_stability_correction(stability_functions,
        stability_parameter_at_height(InverseObukhovLength(), depth, inverse_obukhov_length))

    κ = karman_constant
    u★ = friction_velocity
    ρcₚ = air_volumetric_heat_capacity
    R_top = min((log(reference_depth / top_depth) - heat_correction(reference_depth) + heat_correction(top_depth)) /
        (κ * u★), max_resistance)

    balanced_temperature = reference_temperature + sensible_heat_flux * R_top / ρcₚ
    balanced_vapor_density = reference_vapor_density + vapor_flux * R_top

    residence_time = top_depth * R_top
    fraction_reached = 1.0 - exp(-CANOPY_TIMESTEP / residence_time)
    top_temperature = previous_top_temperature + (balanced_temperature - previous_top_temperature) * fraction_reached
    top_vapor_density = previous_top_vapor_density + (balanced_vapor_density - previous_top_vapor_density) * fraction_reached
    sensible_flux = ρcₚ * (top_temperature - reference_temperature) / R_top
    return (; top_temperature, top_vapor_density, sensible_flux)
end
