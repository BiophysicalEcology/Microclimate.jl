"""
    CanopyWindAttenuation()

Wind speed below canopy top decays from the canopy-top value via a
structural attenuation-shape profile (Massman/Katul-style mixing-length
argument); above and at canopy top, wind is the shared
`AbstractBoundaryLayerModel` (e.g. `MoninObukhov`) evaluated with the
canopy's displacement height and roughness length. No free parameters of
its own.
"""
struct CanopyWindAttenuation <: AbstractCanopyWindModel end

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
    pai = max(plant_area_index, 0.001)  # numerical floor, not a physical parameter
    return (1.0 - (1.0 - exp(-sqrt(7.5 * pai))) / sqrt(7.5 * pai)) * canopy_height
end

"""
    canopy_roughness_length(canopy_height, plant_area_index, displacement_height, karman_constant)

Neutral roughness length from canopy structure. Atmospheric stability is
applied separately by `atmospheric_surface_profile!`'s own Monin-Obukhov
machinery, not folded in here.
"""
function canopy_roughness_length(canopy_height, plant_area_index, displacement_height, karman_constant)
    canopy_element_shelter = sqrt(0.003 + 0.1 * plant_area_index)  # shelter-factor form, empirical, source unconfirmed
    canopy_depth = canopy_height - displacement_height
    roughness = canopy_depth * exp(-karman_constant / canopy_element_shelter)
    return clamp(roughness, 0.0005u"m", 0.9 * canopy_depth)  # free/tunable floor and ceiling, uncited
end

# Attenuation coefficient for one plant-area-index element (structural
# mixing-length argument): larger local plant area density -> more
# attenuation per unit height. 0.003/0.1 shelter-factor and 2.0/0.25
# mixing-length constants are empirical, source unconfirmed.
function _element_attenuation(element_plant_area_index, canopy_height)
    element_shelter = sqrt(0.003 + 0.1 * element_plant_area_index)
    area_density = element_plant_area_index / canopy_height
    mixing_length = 2.0 * element_shelter^3 / (0.25 * area_density)
    return element_shelter * canopy_height / mixing_length
end

"""
    wind_attenuation_profile(layer_plant_area_index, canopy_height, plant_area_index_total, karman_constant)

Shape function giving each layer's wind speed as a fraction of the
canopy-top wind speed (structural, computed once). `layer_plant_area_index`
is ordered top-to-bottom (index 1 = canopy top), matching the radiation
layer convention. Requires at least 10 layers; the bottom tenth uses a
local log-law rather than the shelter recurrence, to avoid the wind speed
collapsing to zero.
"""
function wind_attenuation_profile(layer_plant_area_index, canopy_height, plant_area_index_total, karman_constant)
    n = length(layer_plant_area_index)
    n >= 10 || throw(ArgumentError("wind_attenuation_profile requires at least 10 canopy layers"))

    # internally bottom-to-top (index n = canopy top = unattenuated)
    whole_canopy_attenuation = _element_attenuation(plant_area_index_total, canopy_height)
    layer_attenuation = [_element_attenuation(pai, canopy_height) for pai in Iterators.reverse(layer_plant_area_index)]
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

    # 10.0 here is the same "bottom tenth" ratio as n_ground_layers = n÷10;
    # 20.0 sets ground_roughness ~ half the bottom layer's own thickness.
    ground_roughness = canopy_height / (20.0 * n_ground_layers)
    friction_scale = karman_constant * relative_wind[n_ground_layers] / log(canopy_height / (10.0 * ground_roughness))
    for i in 1:n_ground_layers
        z = i * canopy_height / (10.0 * n_ground_layers)
        relative_wind[i] = (friction_scale / karman_constant) * log(z / ground_roughness)
    end
    return reverse(relative_wind)
end

function allocate_wind(::CanopyWindAttenuation, canopy_height, plant_area_index, boundary_layer_model, heights, n_layers)
    canopy_height = max(canopy_height, 1.0e-3u"m")  # avoid 0/0 in the per-element/roughness formulas below
    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    total_plant_area_index = sum(plant_area_index)

    displacement_height = zero_plane_displacement(canopy_height, total_plant_area_index)
    roughness_length = canopy_roughness_length(
        canopy_height, total_plant_area_index, displacement_height, boundary_layer_model.karman_constant,
    )
    wind_attenuation = wind_attenuation_profile(
        layer_plant_area_index, canopy_height, total_plant_area_index, boundary_layer_model.karman_constant,
    )
    # Ground-most layer's height above the soil surface — not every
    # air_profile_model buffer stores this (KTheoryAirProfile doesn't), so
    # compute it once here where it's cheap (allocation-time only).
    (; layer_heights) = canopy_layer_heights(heights, canopy_height, n_layers)
    ground_layer_height = layer_heights[n_layers]

    return (;
        displacement_height, roughness_length, wind_attenuation, ground_layer_height,
        above_canopy_profile = allocate_profile(boundary_layer_model, [canopy_height, last(heights)]),
        wind_speed = zeros(typeof(0.0u"m/s"), n_layers),
    )
end

function canopy_wind_profile!(buffers, ::CanopyWindAttenuation, boundary_layer_model;
    site, environment_instant, canopy_source_temperature, vapour_pressure_equation=GoffGratch(),
)
    profile = atmospheric_surface_profile!(boundary_layer_model, buffers.above_canopy_profile;
        site, environment_instant, surface_temperature=canopy_source_temperature, vapour_pressure_equation,
        displacement_height=buffers.displacement_height, roughness_length=buffers.roughness_length,
    )
    canopy_top_wind_speed = first(profile.wind_speed)
    canopy_top_air_temperature = first(profile.air_temperature)
    canopy_top_relative_humidity = first(profile.relative_humidity)

    wind_speed = buffers.wind_speed
    @inbounds for i in eachindex(wind_speed)
        wind_speed[i] = buffers.wind_attenuation[i] * canopy_top_wind_speed
    end

    return (; canopy_top_wind_speed, canopy_top_air_temperature, canopy_top_relative_humidity,
        friction_velocity=profile.friction_velocity, obukhov_length=profile.obukhov_length)
end

"""
    canopy_top_flux_boundary(boundary_layer_model, canopy_height, displacement_height, roughness_length,
                              reference_height, friction_velocity, obukhov_length,
                              total_sensible_flux, total_latent_flux, ground_temperature, ground_vapour_density,
                              reference_temperature, reference_vapour_density; max_iter=10, tol=0.01u"K")

Canopy-top boundary temperature/vapour density: resistance-weighted from the
whole canopy's own aggregated source flux (`total_sensible_flux`/
`total_latent_flux`, per unit ground area, summed over every layer) plus a
ground flux, referenced to `reference_temperature`/`reference_vapour_density`
via the canopy-top-to-reference aerodynamic resistance. Ground-to-canopy-top
resistance uses a bulk top-of-canopy eddy diffusivity (same form as
`_canopy_direct_conductances`'s `eddy_diffusivity_top`), not a specific
`air_profile_model`'s own layered resistance -- this boundary is shared by
every `air_profile_model`. Small fixed point (the ground flux depends on the
boundary temperature); converges in a handful of passes.
"""
function canopy_top_flux_boundary(boundary_layer_model, canopy_height, displacement_height, roughness_length,
    reference_height, friction_velocity, obukhov_length,
    total_sensible_flux, total_latent_flux, ground_temperature, ground_vapour_density,
    reference_temperature, reference_vapour_density; max_iter=10, tol=0.01u"K", relax=0.5,
)
    κ = boundary_layer_model.karman_constant
    γ = boundary_layer_model.dyer_constant
    z_top = max(canopy_height - displacement_height, 1.0e-3u"m")
    z_ref = max(reference_height - displacement_height, 1.0e-3u"m")

    ground_resistance = canopy_height / (κ * friction_velocity * z_top)

    log_ratio = log(z_ref / roughness_length)
    stable = !isfinite(obukhov_length) || obukhov_length >= 0.0u"m"
    ψ_h = stable ? 0.0 : calc_ψ_h(calc_φ_m(z_ref, γ, obukhov_length))
    top_resistance = max(log_ratio - ψ_h, 0.1 * log_ratio, 1.0e-6) / (κ * friction_velocity)

    # Self-referential (ground_flux depends on top_temperature) fixed point --
    # top_resistance/ground_resistance can both grow large under near-zero
    # wind, so an unrelaxed step can overshoot badly; damp every step, not
    # just check tolerance on the raw one.
    top_temperature = reference_temperature
    top_vapour_density = reference_vapour_density
    for _ in 1:max_iter
        ρ_cp = calc_ρ_cp(top_temperature)
        ground_flux = (ρ_cp / ground_resistance) * (ground_temperature - top_temperature)
        ground_vapour_flux = (ground_vapour_density - top_vapour_density) / ground_resistance
        raw_temperature = reference_temperature + (total_sensible_flux + ground_flux) * top_resistance / ρ_cp
        raw_vapour_density = reference_vapour_density + (total_latent_flux + ground_vapour_flux) * top_resistance
        new_temperature = top_temperature + relax * (raw_temperature - top_temperature)
        new_vapour_density = top_vapour_density + relax * (raw_vapour_density - top_vapour_density)
        converged = abs(new_temperature - top_temperature) < tol
        top_temperature = new_temperature
        top_vapour_density = new_vapour_density
        converged && break
    end
    # Hard backstop against a runaway fixed point, same bracket-clamp spirit
    # as leaf_temperature's own clamp in _canopy_picard_pass!.
    bound_lo = min(ground_temperature, reference_temperature) - 40.0u"K"
    bound_hi = max(ground_temperature, reference_temperature) + 40.0u"K"
    top_temperature = clamp(top_temperature, bound_lo, bound_hi)
    top_vapour_density = max(top_vapour_density, zero(top_vapour_density))
    return (; top_temperature, top_vapour_density)
end
