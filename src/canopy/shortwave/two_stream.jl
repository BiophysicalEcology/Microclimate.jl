"""
    TwoStreamRadiation(; leaf_reflectance=0.05, leaf_transmittance=0.05)

Dickinson/Sellers two-stream canopy radiative transfer (the same closed-form
solution used in Dickinson (1983), Sellers (1985), and ClimaLand's
`TwoStreamModel`).

- `leaf_reflectance`, `leaf_transmittance` — leaf shortwave optical properties

Leaf angle distribution (`x`) is not stored here — it's a geometric leaf
trait shared with the rain-interception model, so it lives on
[`LeafParameters`](@ref). Ground reflectance is likewise supplied by the
caller to [`canopy_shortwave!`](@ref) (e.g. `Site.albedo`), not stored here.

# References
- Dickinson, R. E. (1983). Land surface processes and climate—surface
  albedos and energy balance. *Advances in Geophysics*, 25, 305–353.
- Sellers, P. J. (1985). Canopy reflectance, photosynthesis and
  transpiration. *International Journal of Remote Sensing*, 6(8), 1335–1372.
- Campbell, G. S. (1990). Derivation of an angle density function for
  canopies with ellipsoidal leaf angle distributions. *Agricultural and
  Forest Meteorology*, 49(3), 173–176.
"""
@kwdef struct TwoStreamRadiation{LR,LT} <: AbstractCanopyShortwaveModel
    leaf_reflectance::LR = 0.05
    leaf_transmittance::LT = 0.05
end

"""
    ellipsoidal_extinction_coefficient(zenith_angle, leaf_angle_distribution_parameter)

Campbell's (1990) ellipsoidal-leaf-angle-distribution extinction coefficient
for a direct beam at `zenith_angle`, for a canopy with leaf angle distribution
parameter `x`. `x == 1` (spherical), `x == 0` (horizontal leaves), and
`x == Inf` (vertical leaves) are closed-form special cases.
"""
function ellipsoidal_extinction_coefficient(zenith_angle, leaf_angle_distribution_parameter)
    x = leaf_angle_distribution_parameter
    if isinf(x)
        return 1.0
    elseif x == 1.0
        return 1.0 / (2.0 * cos(zenith_angle))
    elseif x == 0.0
        return tan(zenith_angle)
    else
        return sqrt(x^2 + tan(zenith_angle)^2) / (x + 1.774 * (x + 1.182)^(-0.733))
    end
end

# Mean leaf-orientation factor J (Campbell 1990); x==1 is closed-form (J=1/3)
function mean_leaf_orientation_factor(leaf_angle_distribution_parameter)
    x = leaf_angle_distribution_parameter
    if x == 1.0
        return 1.0 / 3.0
    else
        mean_leaf_inclination_angle = min(9.65 * (3.0 + x)^(-1.65), π / 2)
        return cos(mean_leaf_inclination_angle)^2
    end
end

"""
    leaf_optics(plant_area_index, leaf_angle_distribution_parameter,
                leaf_reflectance, leaf_transmittance)

Structural two-stream quantities that depend only on canopy leaf optical
properties and total plant area index — not on solar position or ground
reflectance (which can change hour to hour with soil wetness or snow) — so
this is computed once (in `allocate_radiation`) and reused every hour.
"""
function leaf_optics(plant_area_index, leaf_angle_distribution_parameter, leaf_reflectance, leaf_transmittance)
    leaf_scattering_coefficient = leaf_reflectance + leaf_transmittance
    leaf_absorptance = 1.0 - leaf_scattering_coefficient
    leaf_scattering_asymmetry = leaf_reflectance - leaf_transmittance
    leaf_orientation_factor = mean_leaf_orientation_factor(leaf_angle_distribution_parameter)
    diffuse_backscatter_coefficient = 0.5 * (leaf_scattering_coefficient + leaf_orientation_factor * leaf_scattering_asymmetry)
    diffuse_attenuation_rate = sqrt(leaf_absorptance^2 + 2.0 * leaf_absorptance * diffuse_backscatter_coefficient)
    diffuse_transmission_factor = exp(-diffuse_attenuation_rate * plant_area_index)
    return (;
        leaf_scattering_coefficient, leaf_absorptance, leaf_scattering_asymmetry, leaf_orientation_factor,
        diffuse_backscatter_coefficient, diffuse_attenuation_rate, diffuse_transmission_factor,
    )
end

"""
    diffuse_two_stream_optics(leaf_optics, ground_reflectance)

Two-stream boundary-condition coefficients for *diffuse* radiation. Unlike
`leaf_optics`, these depend on ground reflectance, so they're recomputed
every hour rather than cached — cheap scalar ops, no allocation.
"""
function diffuse_two_stream_optics(leaf_optics, ground_reflectance)
    (; leaf_absorptance, diffuse_backscatter_coefficient, diffuse_attenuation_rate, diffuse_transmission_factor) = leaf_optics

    upward_diffuse_boundary_term = leaf_absorptance + diffuse_backscatter_coefficient * (1.0 - 1.0 / ground_reflectance)
    downward_diffuse_boundary_term = leaf_absorptance + diffuse_backscatter_coefficient * (1.0 - ground_reflectance)

    upward_diffuse_normalization =
        (leaf_absorptance + diffuse_backscatter_coefficient + diffuse_attenuation_rate) * (upward_diffuse_boundary_term - diffuse_attenuation_rate) / diffuse_transmission_factor -
        (leaf_absorptance + diffuse_backscatter_coefficient - diffuse_attenuation_rate) * (upward_diffuse_boundary_term + diffuse_attenuation_rate) * diffuse_transmission_factor
    downward_diffuse_normalization =
        (downward_diffuse_boundary_term + diffuse_attenuation_rate) / diffuse_transmission_factor -
        (downward_diffuse_boundary_term - diffuse_attenuation_rate) * diffuse_transmission_factor

    diffuse_reflected_decay_coefficient = (diffuse_backscatter_coefficient / (upward_diffuse_normalization * diffuse_transmission_factor)) * (upward_diffuse_boundary_term - diffuse_attenuation_rate)
    diffuse_reflected_growth_coefficient = (-diffuse_backscatter_coefficient * diffuse_transmission_factor / upward_diffuse_normalization) * (upward_diffuse_boundary_term + diffuse_attenuation_rate)
    diffuse_transmitted_decay_coefficient = (1.0 / (downward_diffuse_normalization * diffuse_transmission_factor)) * (downward_diffuse_boundary_term + diffuse_attenuation_rate)
    diffuse_transmitted_growth_coefficient = (-diffuse_transmission_factor / downward_diffuse_normalization) * (downward_diffuse_boundary_term - diffuse_attenuation_rate)

    return (;
        upward_diffuse_boundary_term, downward_diffuse_boundary_term,
        upward_diffuse_normalization, downward_diffuse_normalization,
        diffuse_reflected_decay_coefficient, diffuse_reflected_growth_coefficient,
        diffuse_transmitted_decay_coefficient, diffuse_transmitted_growth_coefficient,
    )
end

"""
    direct_two_stream_optics(plant_area_index, direct_beam_extinction_coefficient,
                              ground_reflectance, leaf_optics, diffuse_optics)

Two-stream coefficients for the *direct* beam. Depends on solar zenith angle
(via `direct_beam_extinction_coefficient`), so recomputed every hour like
`diffuse_optics`.

Note: the classic Dickinson/Sellers direct-beam particular solution doesn't
exactly satisfy the ground boundary condition (reflected upward flux =
`ground_reflectance` × downward flux at canopy base) — max relative error
in whole-canopy energy conservation is ~6e-4 across a zenith/reflectance
sweep (the diffuse-only solution closes to floating-point precision); see
`test/canopy_radiation_test.jl`.
"""
function direct_two_stream_optics(plant_area_index, direct_beam_extinction_coefficient, ground_reflectance, leaf_optics, diffuse_optics)
    (; leaf_scattering_coefficient, leaf_absorptance, leaf_scattering_asymmetry, leaf_orientation_factor,
       diffuse_backscatter_coefficient, diffuse_attenuation_rate, diffuse_transmission_factor) = leaf_optics
    (; upward_diffuse_boundary_term, downward_diffuse_boundary_term,
       upward_diffuse_normalization, downward_diffuse_normalization) = diffuse_optics
    k = direct_beam_extinction_coefficient

    direct_beam_solution_denominator = k^2 + diffuse_backscatter_coefficient^2 - (leaf_absorptance + diffuse_backscatter_coefficient)^2
    direct_upward_source_term = 0.5 * (leaf_scattering_coefficient + leaf_orientation_factor * leaf_scattering_asymmetry / k) * k
    direct_downward_source_term = leaf_scattering_coefficient * k - direct_upward_source_term
    direct_beam_transmission_factor = exp(-k * plant_area_index)

    direct_reflected_beam_coefficient = -direct_upward_source_term * (leaf_absorptance + diffuse_backscatter_coefficient - k) -
        diffuse_backscatter_coefficient * direct_downward_source_term
    v1 = direct_upward_source_term - direct_reflected_beam_coefficient * (leaf_absorptance + diffuse_backscatter_coefficient + k) / direct_beam_solution_denominator
    v2 = direct_upward_source_term - diffuse_backscatter_coefficient -
        (direct_reflected_beam_coefficient / direct_beam_solution_denominator) * (upward_diffuse_boundary_term + k)
    direct_reflected_decay_coefficient = (1.0 / upward_diffuse_normalization) *
        ((v1 / diffuse_transmission_factor) * (upward_diffuse_boundary_term - diffuse_attenuation_rate) -
         (leaf_absorptance + diffuse_backscatter_coefficient - diffuse_attenuation_rate) * direct_beam_transmission_factor * v2)
    direct_reflected_growth_coefficient = (-1.0 / upward_diffuse_normalization) *
        ((v1 * diffuse_transmission_factor) * (upward_diffuse_boundary_term + diffuse_attenuation_rate) -
         (leaf_absorptance + diffuse_backscatter_coefficient + diffuse_attenuation_rate) * direct_beam_transmission_factor * v2)

    direct_transmitted_beam_coefficient = direct_downward_source_term * (leaf_absorptance + diffuse_backscatter_coefficient + k) -
        diffuse_backscatter_coefficient * direct_upward_source_term
    v3 = (direct_downward_source_term + diffuse_backscatter_coefficient * ground_reflectance -
          (direct_transmitted_beam_coefficient / -direct_beam_solution_denominator) * (downward_diffuse_boundary_term - k)) * direct_beam_transmission_factor
    direct_transmitted_decay_coefficient = (-1.0 / downward_diffuse_normalization) *
        ((direct_transmitted_beam_coefficient / (-direct_beam_solution_denominator * diffuse_transmission_factor)) * (downward_diffuse_boundary_term + diffuse_attenuation_rate) + v3)
    direct_transmitted_growth_coefficient = (1.0 / downward_diffuse_normalization) *
        ((direct_transmitted_beam_coefficient * diffuse_transmission_factor / -direct_beam_solution_denominator) * (downward_diffuse_boundary_term - diffuse_attenuation_rate) + v3)

    return (;
        direct_beam_solution_denominator,
        direct_reflected_beam_coefficient, direct_reflected_decay_coefficient, direct_reflected_growth_coefficient,
        direct_transmitted_beam_coefficient, direct_transmitted_decay_coefficient, direct_transmitted_growth_coefficient,
    )
end

"""
    two_stream_flux_fractions(leaf_optics, diffuse_optics, direct_optics,
                               direct_beam_extinction_coefficient, maximum_layer_reflectance,
                               has_direct_beam, plant_area_index_above)

Diffuse and direct up/down flux fractions (of incoming diffuse/direct
horizontal irradiance respectively) at a given cumulative plant area index
`plant_area_index_above`. Shared between the per-layer profile, the layer
*boundaries* (for the flux-divergence absorption below), and the whole-canopy
ground-level totals, so the clamp/exp algebra is written once.
"""
function two_stream_flux_fractions(leaf_optics, diffuse_optics, direct_optics,
    direct_beam_extinction_coefficient, maximum_layer_reflectance, has_direct_beam, plant_area_index_above,
)
    (; diffuse_attenuation_rate) = leaf_optics
    (; diffuse_reflected_decay_coefficient, diffuse_reflected_growth_coefficient,
       diffuse_transmitted_decay_coefficient, diffuse_transmitted_growth_coefficient) = diffuse_optics
    diffuse_decay = exp(-diffuse_attenuation_rate * plant_area_index_above)
    diffuse_growth = exp(diffuse_attenuation_rate * plant_area_index_above)
    diffuse_up = clamp(diffuse_reflected_decay_coefficient * diffuse_decay + diffuse_reflected_growth_coefficient * diffuse_growth, 0.0, 1.0)
    diffuse_down = clamp(diffuse_transmitted_decay_coefficient * diffuse_decay + diffuse_transmitted_growth_coefficient * diffuse_growth, 0.0, 1.0)

    if has_direct_beam
        (; direct_beam_solution_denominator, direct_reflected_beam_coefficient, direct_reflected_decay_coefficient,
           direct_reflected_growth_coefficient, direct_transmitted_beam_coefficient, direct_transmitted_decay_coefficient,
           direct_transmitted_growth_coefficient) = direct_optics
        direct_transmission = exp(-direct_beam_extinction_coefficient * plant_area_index_above)
        direct_up = clamp(
            direct_reflected_beam_coefficient / direct_beam_solution_denominator * direct_transmission +
            direct_reflected_decay_coefficient * diffuse_decay + direct_reflected_growth_coefficient * diffuse_growth,
            0.0, maximum_layer_reflectance)
        direct_down = clamp(
            direct_transmitted_beam_coefficient / -direct_beam_solution_denominator * direct_transmission +
            direct_transmitted_decay_coefficient * diffuse_decay + direct_transmitted_growth_coefficient * diffuse_growth,
            0.0, maximum_layer_reflectance)
    else
        direct_transmission = 0.0
        direct_up = 0.0
        direct_down = 0.0
    end
    return (; diffuse_up, diffuse_down, direct_up, direct_down, direct_transmission)
end

# Downward (diffuse + direct) and upward irradiance through a horizontal
# plane at a plant-area-index level, given precomputed flux fractions there —
# used at layer boundaries for flux-divergence absorption, and at the canopy
# base/top for ground/sky fluxes. The uncollided direct beam's contribution
# uses `direct_horizontal_irradiance` (the horizontal-plane projection), not
# the beam-normal irradiance — the latter is only meaningful as the
# `downward_direct_shortwave` diagnostic, not for a ground-area energy budget.
function two_stream_irradiances(fractions, direct_horizontal_irradiance, diffuse_horizontal_irradiance)
    downward = fractions.diffuse_down * diffuse_horizontal_irradiance + fractions.direct_down * direct_horizontal_irradiance +
        direct_horizontal_irradiance * fractions.direct_transmission
    upward = fractions.diffuse_up * diffuse_horizontal_irradiance + fractions.direct_up * direct_horizontal_irradiance
    return (; downward, upward)
end

function allocate_shortwave(radiation_model::TwoStreamRadiation, canopy_height, plant_area_index, n_layers, leaf_angle_distribution_parameter)
    # boundary i is the PAI above the top of layer i (layer i spans boundary i to i+1)
    (; layer_plant_area_index_boundaries, layer_plant_area_index_above) = canopy_layer_geometry(plant_area_index, n_layers)

    leaf_optics_precomputed = leaf_optics(
        sum(plant_area_index), leaf_angle_distribution_parameter,
        radiation_model.leaf_reflectance, radiation_model.leaf_transmittance,
    )

    return (;
        layer_plant_area_index_boundaries, layer_plant_area_index_above,
        leaf_angle_distribution_parameter,
        leaf_optics = leaf_optics_precomputed,
        boundary_downward_shortwave = zeros(typeof(0.0u"W/m^2"), n_layers + 1),
        boundary_upward_shortwave = zeros(typeof(0.0u"W/m^2"), n_layers + 1),
        downward_direct_shortwave = zeros(typeof(0.0u"W/m^2"), n_layers),
        downward_diffuse_shortwave = zeros(typeof(0.0u"W/m^2"), n_layers),
        upward_diffuse_shortwave = zeros(typeof(0.0u"W/m^2"), n_layers),
        sunlit_fraction = zeros(n_layers),
        absorbed_shortwave = zeros(typeof(0.0u"W/m^2"), n_layers),
    )
end

function canopy_shortwave!(buffers, radiation_model::TwoStreamRadiation, plant_area_index;
    zenith_angle, direct_horizontal_irradiance, diffuse_horizontal_irradiance, ground_reflectance,
)
    (; layer_plant_area_index_boundaries, layer_plant_area_index_above, boundary_downward_shortwave,
       boundary_upward_shortwave, downward_direct_shortwave, downward_diffuse_shortwave,
       upward_diffuse_shortwave, sunlit_fraction, absorbed_shortwave) = buffers
    n_layers = length(layer_plant_area_index_above)

    fill!(downward_direct_shortwave, 0.0u"W/m^2")
    fill!(downward_diffuse_shortwave, 0.0u"W/m^2")
    fill!(upward_diffuse_shortwave, 0.0u"W/m^2")
    fill!(sunlit_fraction, 0.0)
    fill!(absorbed_shortwave, 0.0u"W/m^2")

    global_horizontal_irradiance = direct_horizontal_irradiance + diffuse_horizontal_irradiance
    if global_horizontal_irradiance <= 0.0u"W/m^2" || zenith_angle >= 90.0u"°"
        return (; ground_absorbed_shortwave = 0.0u"W/m^2", canopy_absorbed_shortwave = 0.0u"W/m^2")
    end

    leaf = buffers.leaf_optics
    diffuse = diffuse_two_stream_optics(leaf, ground_reflectance)
    direct_beam_extinction_coefficient = ellipsoidal_extinction_coefficient(zenith_angle, buffers.leaf_angle_distribution_parameter)
    has_direct_beam = direct_horizontal_irradiance > 0.0u"W/m^2"
    direct_beam_irradiance = has_direct_beam ? direct_horizontal_irradiance / cos(zenith_angle) : 0.0u"W/m^2"
    maximum_layer_reflectance = max(ground_reflectance, radiation_model.leaf_reflectance)
    direct = direct_two_stream_optics(sum(plant_area_index), direct_beam_extinction_coefficient, ground_reflectance, leaf, diffuse)

    # boundary fluxes, reused across adjacent layers (boundary i is shared by layers i-1 and i)
    @inbounds for i in eachindex(layer_plant_area_index_boundaries)
        fractions = two_stream_flux_fractions(leaf, diffuse, direct, direct_beam_extinction_coefficient,
            maximum_layer_reflectance, has_direct_beam, layer_plant_area_index_boundaries[i])
        irradiances = two_stream_irradiances(fractions, direct_horizontal_irradiance, diffuse_horizontal_irradiance)
        boundary_downward_shortwave[i] = irradiances.downward
        boundary_upward_shortwave[i] = irradiances.upward
    end

    @inbounds for i in 1:n_layers
        absorbed_shortwave[i] = (boundary_downward_shortwave[i] - boundary_downward_shortwave[i + 1]) +
                                 (boundary_upward_shortwave[i + 1] - boundary_upward_shortwave[i])

        mid_fractions = two_stream_flux_fractions(leaf, diffuse, direct, direct_beam_extinction_coefficient,
            maximum_layer_reflectance, has_direct_beam, layer_plant_area_index_above[i])
        downward_diffuse_shortwave[i] = mid_fractions.diffuse_down * diffuse_horizontal_irradiance + mid_fractions.direct_down * direct_horizontal_irradiance
        upward_diffuse_shortwave[i] = mid_fractions.diffuse_up * diffuse_horizontal_irradiance + mid_fractions.direct_up * direct_horizontal_irradiance
        downward_direct_shortwave[i] = direct_beam_irradiance * mid_fractions.direct_transmission
        sunlit_fraction[i] = mid_fractions.direct_transmission
    end

    ground_absorbed_shortwave = (1.0 - ground_reflectance) * boundary_downward_shortwave[n_layers + 1]
    canopy_absorbed_shortwave = global_horizontal_irradiance - boundary_upward_shortwave[1]

    return (; ground_absorbed_shortwave, canopy_absorbed_shortwave)
end
