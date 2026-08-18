"""
    MixingLengthCanopyWindAttenuation(; thermal_roughness_model=ScalarRoughnessRatio(),
                            shelter_floor=0.003, shelter_pai_coefficient=0.1,
                            mixing_length_coefficient=2.0, mixing_length_pai_coefficient=0.25)

Wind speed below canopy top decays from the canopy-top value via a
structural attenuation-shape profile (Massman/Katul-style mixing-length
argument); above and at canopy top, wind is the shared
`AbstractBoundaryLayerModel` (e.g. `MoninObukhov`) evaluated with the
canopy's displacement height and roughness length.

`thermal_roughness_model` overrides the shared `boundary_layer_model`'s own
`thermal_roughness_model` for this canopy-top evaluation specifically (bare-
ground/soil surface fluxes elsewhere keep using the `boundary_layer_model`'s
own setting). Defaults to `ScalarRoughnessRatio()`.

`shelter_floor`/`shelter_pai_coefficient`/`mixing_length_coefficient`/
`mixing_length_pai_coefficient` are the empirical shelter-factor/mixing-
length constants in [`_element_attenuation`](@ref)/[`canopy_roughness_length`](@ref).
Can over-attenuate wind at low total PAI -- tune per site if so.
"""
@kwdef struct MixingLengthCanopyWindAttenuation{TRM,SF,SPC,MLC,MLPC} <: AbstractCanopyWindModel
    thermal_roughness_model::TRM = ScalarRoughnessRatio()
    shelter_floor::SF = 0.003
    shelter_pai_coefficient::SPC = 0.1
    mixing_length_coefficient::MLC = 2.0
    mixing_length_pai_coefficient::MLPC = 0.25
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
    canopy_roughness_length(canopy_height, plant_area_index, displacement_height, karman_constant;
                             shelter_floor=0.003, shelter_pai_coefficient=0.1)

Neutral roughness length from canopy structure. Atmospheric stability is
applied separately by `atmospheric_surface_profile!`'s own Monin-Obukhov
machinery, not folded in here.
"""
function canopy_roughness_length(canopy_height, plant_area_index, displacement_height, karman_constant;
    shelter_floor=0.003, shelter_pai_coefficient=0.1,
)
    canopy_element_shelter = sqrt(shelter_floor + shelter_pai_coefficient * plant_area_index)
    canopy_depth = canopy_height - displacement_height
    roughness = canopy_depth * exp(-karman_constant / canopy_element_shelter)
    return clamp(roughness, 0.0005u"m", 0.9 * canopy_depth)  # free/tunable floor and ceiling, uncited
end

# Attenuation coefficient for one plant-area-index element (structural
# mixing-length argument): larger local plant area density -> more
# attenuation per unit height.
function _element_attenuation(element_plant_area_index, canopy_height,
    shelter_floor, shelter_pai_coefficient, mixing_length_coefficient, mixing_length_pai_coefficient,
)
    element_shelter = sqrt(shelter_floor + shelter_pai_coefficient * element_plant_area_index)
    area_density = element_plant_area_index / canopy_height
    mixing_length = mixing_length_coefficient * element_shelter^3 / (mixing_length_pai_coefficient * area_density)
    return element_shelter * canopy_height / mixing_length
end

"""
    wind_attenuation_profile(layer_plant_area_index, canopy_height, plant_area_index_total, karman_constant;
                              shelter_floor=0.003, shelter_pai_coefficient=0.1,
                              mixing_length_coefficient=2.0, mixing_length_pai_coefficient=0.25)

Shape function giving each layer's wind speed as a fraction of the
canopy-top wind speed (structural, computed once). `layer_plant_area_index`
is ordered top-to-bottom (index 1 = canopy top), matching the radiation
layer convention. Requires at least 10 layers; the bottom tenth uses a
local log-law rather than the shelter recurrence, to avoid the wind speed
collapsing to zero.
"""
function wind_attenuation_profile(layer_plant_area_index, canopy_height, plant_area_index_total, karman_constant;
    shelter_floor=0.003, shelter_pai_coefficient=0.1, mixing_length_coefficient=2.0, mixing_length_pai_coefficient=0.25,
)
    n = length(layer_plant_area_index)
    n >= 10 || throw(ArgumentError("wind_attenuation_profile requires at least 10 canopy layers"))

    # internally bottom-to-top (index n = canopy top = unattenuated)
    whole_canopy_attenuation = _element_attenuation(plant_area_index_total, canopy_height,
        shelter_floor, shelter_pai_coefficient, mixing_length_coefficient, mixing_length_pai_coefficient)
    layer_attenuation = [_element_attenuation(pai, canopy_height,
        shelter_floor, shelter_pai_coefficient, mixing_length_coefficient, mixing_length_pai_coefficient)
        for pai in Iterators.reverse(layer_plant_area_index)]
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
    ground_roughness = canopy_height / (20.0 * n_ground_layers)
    friction_scale = karman_constant * relative_wind[n_ground_layers] / log(canopy_height / (10.0 * ground_roughness))
    for i in 1:n_ground_layers
        z = i * canopy_height / (10.0 * n_ground_layers)
        relative_wind[i] = (friction_scale / karman_constant) * log(z / ground_roughness)
    end
    return relative_wind
end

function allocate_wind(model::MixingLengthCanopyWindAttenuation, canopy_height, plant_area_index, boundary_layer_model, heights, n_layers)
    canopy_height = max(canopy_height, 1.0e-3u"m")  # avoid 0/0 in the per-element/roughness formulas below
    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    total_plant_area_index = sum(plant_area_index)

    (; shelter_floor, shelter_pai_coefficient, mixing_length_coefficient, mixing_length_pai_coefficient) = model

    displacement_height = zero_plane_displacement(canopy_height, total_plant_area_index)
    roughness_length = canopy_roughness_length(
        canopy_height, total_plant_area_index, displacement_height, boundary_layer_model.karman_constant;
        shelter_floor, shelter_pai_coefficient,
    )
    wind_attenuation = wind_attenuation_profile(
        layer_plant_area_index, canopy_height, total_plant_area_index, boundary_layer_model.karman_constant;
        shelter_floor, shelter_pai_coefficient, mixing_length_coefficient, mixing_length_pai_coefficient,
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
)
    profile = atmospheric_surface_profile!(boundary_layer_model, buffers.above_canopy_profile;
        site, environment_instant, surface_temperature=canopy_source_temperature, vapour_pressure_equation,
        displacement_height=buffers.displacement_height, roughness_length=buffers.roughness_length,
        thermal_roughness_model=buffers.thermal_roughness_model,
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
        friction_velocity=profile.friction_velocity, friction_velocity_raw=profile.friction_velocity_raw,
        obukhov_length=profile.obukhov_length)
end

"""
    canopy_top_flux_boundary(boundary_layer_model, canopy_height, displacement_height, roughness_length,
                              reference_height, friction_velocity, obukhov_length,
                              total_sensible_flux, total_latent_flux, ground_temperature, ground_vapour_density,
                              reference_temperature, reference_vapour_density, previous_top_temperature,
                              previous_top_vapour_density; max_resistance=2000.0u"s/m")

Canopy-top boundary temperature/vapour density: resistance-weighted from the
whole canopy's own aggregated source flux (`total_sensible_flux`/
`total_latent_flux`, per unit ground area, summed over every layer) plus a
ground flux, referenced to `reference_temperature`/`reference_vapour_density`
via the canopy-top-to-reference aerodynamic resistance. Ground-to-canopy-top
resistance uses a bulk top-of-canopy eddy diffusivity (same form as
`_canopy_direct_conductances`'s `eddy_diffusivity_top`), not a specific
`air_profile_model`'s own layered resistance -- this boundary is shared by
every `air_profile_model`. With `ρ_cp` and both resistances fixed (not
themselves iterated), the balance is linear in `top_temperature`/
`top_vapour_density`, so it's solved directly rather than iteratively.

`top_resistance`/`ground_resistance` grow as `1/friction_velocity`, unbounded
under calm/stable conditions; `top_resistance` directly multiplies
`total_sensible_flux`, so that alone can push `top_temperature` to the `±40K`
clamp. `max_resistance` caps both at a literature-typical ceiling (Thom 1975;
Monteith & Unsworth 2013).

Even at that cap, `top_resistance` can imply a residence time
(`z_top·top_resistance`) longer than `CANOPY_TIMESTEP` -- under calm, dense
canopy conditions the algebraic steady state can then land well within the
±40K bound but still be physically unreachable in an hour. Damped
exponentially toward `previous_top_temperature`/`previous_top_vapour_density`
by that residence time before the bound is applied.
"""
function canopy_top_flux_boundary(boundary_layer_model, canopy_height, displacement_height, roughness_length,
    reference_height, friction_velocity, obukhov_length,
    total_sensible_flux, total_latent_flux, ground_temperature, ground_vapour_density,
    reference_temperature, reference_vapour_density, previous_top_temperature, previous_top_vapour_density;
    max_resistance=2000.0u"s/m",
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

    # ρ_cp fixed from ambient data, not the solved top_temperature: calc_ρ_cp
    # ∝ 1/T, so tying it to the unknown would make the balance nonlinear again.
    ρ_cp = calc_ρ_cp((reference_temperature + ground_temperature) / 2)

    # Bracket, same spirit as leaf_temperature's own clamp.
    # Open issue: :bulk-mode RaupachLTheoryAirProfile can jump 5-6K between
    # layer 1 and layer 2 under low friction_velocity, driving an extreme
    # total_sensible_flux that feeds back into this clamp.
    bound_lo = min(ground_temperature, reference_temperature) - 40.0u"K"
    bound_hi = max(ground_temperature, reference_temperature) + 40.0u"K"

    # Closed-form solution of the linear balance:
    #   Tt = Tr + (H + (ρcp/Rg)(Tg-Tt)) Rt/ρcp   =>   Tt(1+Rt/Rg) = Tr + H·Rt/ρcp + (Rt/Rg)Tg
    # and the analogous vapour-density equation.
    Rtg = top_resistance / ground_resistance
    undamped_top_temperature = (reference_temperature + total_sensible_flux * top_resistance / ρ_cp + Rtg * ground_temperature) / (1.0 + Rtg)
    undamped_top_vapour_density = (reference_vapour_density + total_latent_flux * top_resistance + Rtg * ground_vapour_density) / (1.0 + Rtg)
    top_temperature = undamped_top_temperature
    top_vapour_density = undamped_top_vapour_density

    # Residence time of a z_top-deep column via top_resistance; can exceed an
    # hour under calm conditions (see docstring).
    τ = z_top * top_resistance
    physical_weight = 1.0 - exp(-ustrip(u"s", CANOPY_TIMESTEP) / ustrip(u"s", τ))
    top_temperature = previous_top_temperature + (top_temperature - previous_top_temperature) * physical_weight
    top_vapour_density = previous_top_vapour_density + (top_vapour_density - previous_top_vapour_density) * physical_weight

    # Final backstop: the damping step above can walk back outside the
    # bracket if previous_top_temperature itself sits near/at the edge.
    clamped = top_temperature < bound_lo ? :cold : top_temperature > bound_hi ? :hot : :none
    top_temperature = clamp(top_temperature, bound_lo, bound_hi)
    top_vapour_density = max(top_vapour_density, zero(top_vapour_density))
    return (; top_temperature, top_vapour_density, bound_lo, bound_hi, clamped,
        undamped_top_temperature, previous_top_temperature, physical_weight, τ, top_resistance, ground_resistance)
end
