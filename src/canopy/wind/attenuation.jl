"""
    CanopyWindAttenuation(; thermal_roughness_model=ScalarRoughnessRatio())

Wind speed below canopy top decays from the canopy-top value via a
structural attenuation-shape profile (Massman/Katul-style mixing-length
argument); above and at canopy top, wind is the shared
`AbstractBoundaryLayerModel` (e.g. `MoninObukhov`) evaluated with the
canopy's displacement height and roughness length.

`thermal_roughness_model` overrides the shared `boundary_layer_model`'s own
`thermal_roughness_model` for this canopy-top evaluation specifically (bare-
ground/soil surface fluxes elsewhere keep using the `boundary_layer_model`'s
own setting). Defaults to `ScalarRoughnessRatio()`.
"""
@kwdef struct CanopyWindAttenuation{TRM} <: AbstractCanopyWindModel
    thermal_roughness_model::TRM = ScalarRoughnessRatio()
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

function allocate_wind(model::CanopyWindAttenuation, canopy_height, plant_area_index, boundary_layer_model, heights, n_layers)
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
        canopy_height, displacement_height, roughness_length, wind_attenuation, ground_layer_height,
        above_canopy_profile = allocate_profile(boundary_layer_model, [canopy_height, last(heights)]),
        wind_speed = zeros(typeof(0.0u"m/s"), n_layers),
        thermal_roughness_model = model.thermal_roughness_model,
    )
end

function canopy_wind_profile!(buffers, model::CanopyWindAttenuation, boundary_layer_model;
    site, environment_instant, canopy_source_temperature, vapour_pressure_equation=GoffGratch(),
)
    profile = atmospheric_surface_profile!(boundary_layer_model, buffers.above_canopy_profile;
        site, environment_instant, surface_temperature=canopy_source_temperature, vapour_pressure_equation,
        displacement_height=buffers.displacement_height, roughness_length=buffers.roughness_length,
        thermal_roughness_model=model.thermal_roughness_model,
        temperature_anchor_height=buffers.canopy_height,
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
                              reference_temperature, reference_vapour_density; max_iter=10, tol=0.01u"K",
                              relax=0.5, max_resistance=2000.0u"s/m")

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

`top_resistance`/`ground_resistance` grow as `1/friction_velocity`, unbounded
under calm/stable conditions; `top_resistance` directly multiplies
`total_sensible_flux`, so that alone can push `top_temperature` to the `±40K`
clamp. `max_resistance` caps both at a literature-typical ceiling (Thom 1975;
Monteith & Unsworth 2013).
"""
function canopy_top_flux_boundary(boundary_layer_model, canopy_height, displacement_height, roughness_length,
    reference_height, friction_velocity, obukhov_length,
    total_sensible_flux, total_latent_flux, ground_temperature, ground_vapour_density,
    reference_temperature, reference_vapour_density; max_iter=10, tol=0.01u"K", relax=0.5, max_resistance=2000.0u"s/m",
)
    κ = boundary_layer_model.karman_constant
    γ = boundary_layer_model.dyer_constant
    z_top = max(canopy_height - displacement_height, 1.0e-3u"m")
    z_ref = max(reference_height - displacement_height, 1.0e-3u"m")

    ground_resistance = min(canopy_height / (κ * friction_velocity * z_top), max_resistance)

    log_ratio = log(z_ref / roughness_length)
    ψ_h = isfinite(obukhov_length) ?
        calc_ψ_h(z_ref, γ, obukhov_length, boundary_layer_model.stable_beta, boundary_layer_model.turbulent_prandtl_number) : 0.0
    top_resistance = min(max(log_ratio - ψ_h, 0.1 * log_ratio, 1.0e-6) / (κ * friction_velocity), max_resistance)

    # ρ_cp fixed from ambient data, not the iterate: calc_ρ_cp ∝ 1/T, so
    # recomputing it from top_temperature each pass lets one oversized step
    # near T=0 blow it up or flip its sign, diverging regardless of ω.
    ρ_cp = calc_ρ_cp((reference_temperature + ground_temperature) / 2)
    top_temperature = reference_temperature
    top_vapour_density = reference_vapour_density
    # Aitken Δ² acceleration: ω adapts each pass from the residual history,
    # damping the self-referential ground_flux/top_temperature loop under
    # low-wind (high-resistance) overshoot.
    ω_temp, ω_vapour = relax, relax
    residual_temp_prev = zero(top_temperature)
    residual_vapour_prev = zero(top_vapour_density)
    step_temp = zero(top_temperature)
    for iter in 1:max_iter
        ground_flux = (ρ_cp / ground_resistance) * (ground_temperature - top_temperature)
        ground_vapour_flux = (ground_vapour_density - top_vapour_density) / ground_resistance
        residual_temp = reference_temperature + (total_sensible_flux + ground_flux) * top_resistance / ρ_cp - top_temperature
        residual_vapour = reference_vapour_density + (total_latent_flux + ground_vapour_flux) * top_resistance - top_vapour_density
        if iter > 1
            Δr_temp = residual_temp - residual_temp_prev
            Δr_vapour = residual_vapour - residual_vapour_prev
            iszero(Δr_temp) || (ω_temp = clamp(-ω_temp * (residual_temp_prev / Δr_temp), 0.02, 0.9))
            iszero(Δr_vapour) || (ω_vapour = clamp(-ω_vapour * (residual_vapour_prev / Δr_vapour), 0.02, 0.9))
        end
        step_temp = ω_temp * residual_temp
        top_temperature += step_temp
        top_vapour_density += ω_vapour * residual_vapour
        residual_temp_prev = residual_temp
        residual_vapour_prev = residual_vapour
        abs(step_temp) < tol && break
    end
    # Hard backstop against a runaway fixed point, same bracket-clamp spirit
    # as leaf_temperature's own clamp in _canopy_picard_pass!.
    # Open issue (kept as a warning, not fixed): under low friction_velocity,
    # :bulk-mode RaupachLTheoryAirProfile pins layer 1 (canopy top) exactly
    # to T_top but jumps 5-6K to layer 2 instead of transitioning smoothly
    # (unlike R's LangrangianOne, whose CfTh-CnTh+CnT near-field kernel
    # guarantees continuity there) -- that artificial layer-1/2 ΔT drives an
    # extreme total_sensible_flux, which feeds back into this clamp.
    bound_lo = min(ground_temperature, reference_temperature) - 40.0u"K"
    bound_hi = max(ground_temperature, reference_temperature) + 40.0u"K"
    if top_temperature < bound_lo
        @warn "canopy_top_flux_boundary: hit the ±40K clamp (cold)" top_temperature friction_velocity top_resistance ground_resistance total_sensible_flux total_latent_flux ground_temperature reference_temperature maxlog=40 _id=:canopy_top_clamp_cold
    elseif top_temperature > bound_hi
        @warn "canopy_top_flux_boundary: hit the ±40K clamp (hot)" top_temperature friction_velocity top_resistance ground_resistance total_sensible_flux total_latent_flux ground_temperature reference_temperature maxlog=40 _id=:canopy_top_clamp_hot
    end
    top_temperature = clamp(top_temperature, bound_lo, bound_hi)
    top_vapour_density = max(top_vapour_density, zero(top_vapour_density))
    return (; top_temperature, top_vapour_density, bound_lo, bound_hi)
end
