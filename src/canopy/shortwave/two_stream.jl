"""
    TwoStreamRadiation(; leaf_reflectance=0.25, leaf_transmittance=0.25)

Dickinson/Sellers two-stream canopy radiative transfer: how direct and diffuse
sunlight is split between leaf layers and the ground.

- `leaf_reflectance`, `leaf_transmittance` — leaf shortwave optical properties.

Leaf angle comes from [`LeafParameters`](@ref), ground reflectance from the caller.

# References
- Dickinson, R. E. (1983). Land surface processes and climate—surface
  albedos and energy balance. *Advances in Geophysics*, 25, 305–353.
- Sellers, P. J. (1985). Canopy reflectance, photosynthesis and
  transpiration. *International Journal of Remote Sensing*, 6(8), 1335–1372.
"""
@kwdef struct TwoStreamRadiation{LR,LT} <: AbstractCanopyShortwaveModel
    leaf_reflectance::LR = 0.25
    leaf_transmittance::LT = 0.25
end

"""
    ellipsoidal_extinction_coefficient(zenith_angle, canopy_projection_ratio)

Campbell's (1990) ellipsoidal-leaf-angle-distribution extinction coefficient
for a direct beam at `zenith_angle`, for a canopy where the ratio of horizontal
to vertical projections of canopy elements is `canopy_projection_ratio`.
A ratio of 1 (spherical) and `Inf` (horizontal leaves) are closed-form special
cases.

# References
- Campbell, G. S. (1990). Derivation of an angle density function for
  canopies with ellipsoidal leaf angle distributions. *Agricultural and
  Forest Meteorology*, 49(3), 173–176.
"""
function ellipsoidal_extinction_coefficient(zenith_angle, canopy_projection_ratio)
    χ = canopy_projection_ratio
    if isinf(χ)
        return 1.0
    elseif χ == 1.0
        return 1.0 / (2.0 * cos(zenith_angle))
    else
        Λ = χ + 1.774 * (χ + 1.182)^(-0.733) # Campbell eq. 14
        return sqrt(χ^2 + tan(zenith_angle)^2) / Λ
    end
end

"""
    leaf_orientation_factor(canopy_projection_ratio)

Squared cosine of the mean leaf angle (Campbell 1990, eq. 16), or 1/3 for spherical
leaves (`canopy_projection_ratio == 1`).

# References
- Campbell, G. S. (1990). Derivation of an angle density function for
  canopies with ellipsoidal leaf angle distributions. *Agricultural and
  Forest Meteorology*, 49(3), 173–176.
"""
function leaf_orientation_factor(canopy_projection_ratio)
    χ = canopy_projection_ratio
    if χ == 1.0
        return 1.0 / 3.0
    else
        mean_leaf_angle = min(9.65 * (3.0 + χ)^(-1.65), π / 2)
        return cos(mean_leaf_angle)^2
    end
end

"""
    leaf_optics(plant_area_index, canopy_projection_ratio,
                leaf_reflectance, leaf_transmittance)

Structural two-stream quantities that depend only on canopy leaf optical
properties and total plant area index, not on solar position or ground
reflectance. Compute once.
"""
function leaf_optics(plant_area_index, canopy_projection_ratio, leaf_reflectance, leaf_transmittance)
    PAI = plant_area_index
    ρₗ = leaf_reflectance
    τₗ = leaf_transmittance

    ω = ρₗ + τₗ # scattered fraction
    α = 1.0 - ω # absorbed fraction
    ζ = ρₗ - τₗ # scattering asymmetry

    g = leaf_orientation_factor(canopy_projection_ratio)
    β = 0.5 * (ω + g * ζ) # effective diffuse backscatter coefficient
    λ = sqrt(α^2 + 2.0 * α * β) # diffuse attenuation rate

    τ = exp(-λ * PAI) # whole-canopy diffuse transmission

    return (;
        leaf_scattering_coefficient = ω,
        leaf_absorptance = α,
        leaf_scattering_asymmetry = ζ,
        leaf_orientation_factor = g,
        diffuse_backscatter_coefficient = β,
        diffuse_attenuation_rate = λ,
        diffuse_transmission_factor = τ,
    )
end

# Floors |x| away from zero, sign-preserving (0 -> +ε) -- a plain max() floor
# would flip the sign of a legitimately-negative value instead.
_safe_nonzero(x, ε=1.0e-6) = abs(x) < ε ? (x < 0 ? -ε : ε) : x

"""
    diffuse_two_stream_optics(leaf_optics, ground_reflectance)

Two-stream boundary-condition coefficients for *diffuse* radiation, which
depend on ground reflectance.
"""
function diffuse_two_stream_optics(leaf_optics, ground_reflectance)
    (; leaf_absorptance, diffuse_backscatter_coefficient, diffuse_attenuation_rate, diffuse_transmission_factor) = leaf_optics

    α = leaf_absorptance
    β = diffuse_backscatter_coefficient
    λ = diffuse_attenuation_rate
    τ = diffuse_transmission_factor
    ρg = ground_reflectance

    Bᵤ = α + β * (1.0 - 1.0 / ρg)
    B_d = α + β * (1.0 - ρg)

    # divided into below and in direct_two_stream_optics -- can hit zero
    Nᵤ = _safe_nonzero((α + β + λ) * (Bᵤ - λ) / τ - (α + β - λ) * (Bᵤ + λ) * τ)
    N_d = _safe_nonzero((B_d + λ) / τ - (B_d - λ) * τ)

    diffuse_reflected_decay_coefficient = β * (Bᵤ - λ) / (Nᵤ * τ)
    diffuse_reflected_growth_coefficient = -β * τ * (Bᵤ + λ) / Nᵤ

    diffuse_transmitted_decay_coefficient = (B_d + λ) / (N_d * τ)
    diffuse_transmitted_growth_coefficient = -τ * (B_d - λ) / N_d

    return (;
        upward_diffuse_boundary_term = Bᵤ,
        downward_diffuse_boundary_term = B_d,
        upward_diffuse_normalization = Nᵤ,
        downward_diffuse_normalization = N_d,
        diffuse_reflected_decay_coefficient,
        diffuse_reflected_growth_coefficient,
        diffuse_transmitted_decay_coefficient,
        diffuse_transmitted_growth_coefficient,
    )
end

"""
    direct_two_stream_optics(plant_area_index, direct_beam_extinction_coefficient,
                              ground_reflectance, leaf_optics, diffuse_optics)

Two-stream coefficients for the *direct* beam as a function of solar zenith angle.
"""
function direct_two_stream_optics(plant_area_index, direct_beam_extinction_coefficient, ground_reflectance, leaf_optics, diffuse_optics)
    (; leaf_scattering_coefficient, leaf_absorptance, leaf_scattering_asymmetry, leaf_orientation_factor,
       diffuse_backscatter_coefficient, diffuse_attenuation_rate, diffuse_transmission_factor) = leaf_optics
    (; upward_diffuse_boundary_term, downward_diffuse_boundary_term,
       upward_diffuse_normalization, downward_diffuse_normalization) = diffuse_optics

    # Optical properties and structural terms.
    PAI = plant_area_index
    kᵦ = direct_beam_extinction_coefficient
    ρg = ground_reflectance

    ω = leaf_scattering_coefficient
    α = leaf_absorptance
    ζ = leaf_scattering_asymmetry
    g = leaf_orientation_factor

    β = diffuse_backscatter_coefficient
    λ = diffuse_attenuation_rate
    τ = diffuse_transmission_factor

    # Diffuse-solution boundary and normalization terms.
    Bᵤ = upward_diffuse_boundary_term
    B_d = downward_diffuse_boundary_term
    Nᵤ = upward_diffuse_normalization
    N_d = downward_diffuse_normalization

    # Particular-solution denominator. This has a known Dickinson/Sellers
    # singularity at particular values of kᵦ encountered as zenith angle varies.
    Δ = _safe_nonzero(kᵦ^2 + β^2 - (α + β)^2)

    # Direct-beam source terms for upward and downward diffuse radiation.
    Sᵤ = 0.5 * (ω + g * ζ / kᵦ) * kᵦ
    S_d = ω * kᵦ - Sᵤ

    # Direct-beam transmission through the complete canopy.
    Tᵦ = exp(-kᵦ * PAI)

    # Reflected direct-beam particular-solution coefficient.
    Cᵣ = -Sᵤ * (α + β - kᵦ) - β * S_d
    V₁ = Sᵤ - Cᵣ * (α + β + kᵦ) / Δ
    V₂ = Sᵤ - β - (Cᵣ / Δ) * (Bᵤ + kᵦ)

    # Coefficients of the decaying and growing reflected solutions.
    Rᵦ⁻ = ((V₁ / τ) * (Bᵤ - λ) - (α + β - λ) * Tᵦ * V₂) / Nᵤ
    Rᵦ⁺ = -((V₁ * τ) * (Bᵤ + λ) - (α + β + λ) * Tᵦ * V₂) / Nᵤ

    # Transmitted direct-beam particular-solution coefficient.
    Cₜ = S_d * (α + β + kᵦ) + β * Sᵤ

    V₃ = (S_d + β * ρg - (Cₜ / -Δ) * (B_d - kᵦ)) * Tᵦ

    # Coefficients of the decaying and growing transmitted solutions.
    Tᵦ⁻ = -((Cₜ / (-Δ * τ)) * (B_d + λ) + V₃) / N_d
    Tᵦ⁺ = ((Cₜ * τ / -Δ) * (B_d - λ) + V₃) / N_d

    return (;
        direct_beam_solution_denominator = Δ,
        direct_reflected_beam_coefficient = Cᵣ,
        direct_reflected_decay_coefficient = Rᵦ⁻,
        direct_reflected_growth_coefficient = Rᵦ⁺,
        direct_transmitted_beam_coefficient = Cₜ,
        direct_transmitted_decay_coefficient = Tᵦ⁻,
        direct_transmitted_growth_coefficient = Tᵦ⁺,
    )
end

"""
    two_stream_flux_fractions(leaf_optics, diffuse_optics, direct_optics,
                               direct_beam_extinction_coefficient, maximum_layer_reflectance,
                               has_direct_beam, plant_area_index_above)

Calculate upward and downward two-stream flux fractions at cumulative plant
area index `plant_area_index_above`.

Diffuse fluxes are expressed as fractions of incoming diffuse horizontal
irradiance. Direct-beam-generated fluxes and direct transmission are expressed
as fractions of incoming direct horizontal irradiance.
"""
function two_stream_flux_fractions(leaf_optics, diffuse_optics, direct_optics,
    direct_beam_extinction_coefficient, maximum_layer_reflectance, has_direct_beam, plant_area_index_above,
)
    (; diffuse_attenuation_rate) = leaf_optics
    (; diffuse_reflected_decay_coefficient, diffuse_reflected_growth_coefficient,
       diffuse_transmitted_decay_coefficient, diffuse_transmitted_growth_coefficient) = diffuse_optics

    L = plant_area_index_above
    λ = diffuse_attenuation_rate
    kᵦ = direct_beam_extinction_coefficient
    ρmax = maximum_layer_reflectance
    R⁻ = diffuse_reflected_decay_coefficient
    R⁺ = diffuse_reflected_growth_coefficient
    T⁻ = diffuse_transmitted_decay_coefficient
    T⁺ = diffuse_transmitted_growth_coefficient

    E⁻ = exp(-λ * L) # diffuse decay
    E⁺ = exp(λ * L) # diffuse growth

    diffuse_up = clamp(R⁻ * E⁻ + R⁺ * E⁺, 0.0, 1.0)
    diffuse_down = clamp(T⁻ * E⁻ + T⁺ * E⁺, 0.0, 1.0)

    if has_direct_beam
        (; direct_beam_solution_denominator, direct_reflected_beam_coefficient, direct_reflected_decay_coefficient,
           direct_reflected_growth_coefficient, direct_transmitted_beam_coefficient, direct_transmitted_decay_coefficient,
           direct_transmitted_growth_coefficient) = direct_optics

        Δ = direct_beam_solution_denominator
        Cᵣ = direct_reflected_beam_coefficient
        Cₜ = direct_transmitted_beam_coefficient
        Rᵦ⁻ = direct_reflected_decay_coefficient
        Rᵦ⁺ = direct_reflected_growth_coefficient
        Tᵦ⁻ = direct_transmitted_decay_coefficient
        Tᵦ⁺ = direct_transmitted_growth_coefficient

        Eᵦ = exp(-kᵦ * L)
        direct_up = clamp((Cᵣ / Δ) * Eᵦ + Rᵦ⁻ * E⁻ + Rᵦ⁺ * E⁺, 0.0, ρmax)
        direct_down = clamp(-(Cₜ / Δ) * Eᵦ + Tᵦ⁻ * E⁻ + Tᵦ⁺ * E⁺, 0.0, ρmax)
    else
        Eᵦ = 0.0
        direct_up = 0.0
        direct_down = 0.0
    end
    return (; diffuse_up, diffuse_down, direct_up, direct_down, direct_transmission = Eᵦ)
end

# Downward (diffuse + direct) and upward irradiance through a horizontal
# plane at a plant-area-index level, given precomputed flux fractions there —
# used at layer boundaries for flux-divergence absorption, and at the canopy
# base/top for ground/sky fluxes.
function two_stream_irradiances(fractions, direct_horizontal_irradiance, diffuse_horizontal_irradiance)
    downward = fractions.diffuse_down * diffuse_horizontal_irradiance + fractions.direct_down * direct_horizontal_irradiance +
        direct_horizontal_irradiance * fractions.direct_transmission
    upward = fractions.diffuse_up * diffuse_horizontal_irradiance + fractions.direct_up * direct_horizontal_irradiance
    return (; downward, upward)
end

function allocate_shortwave(radiation_model::TwoStreamRadiation, canopy_height, plant_area_index, n_layers, canopy_projection_ratio)
    # boundary i is the PAI above the top of layer i (layer i spans boundary i to i+1)
    (; layer_plant_area_index_boundaries, layer_plant_area_index_above) = canopy_layer_geometry(plant_area_index, n_layers)

    leaf_optics_precomputed = leaf_optics(
        sum(plant_area_index), canopy_projection_ratio,
        radiation_model.leaf_reflectance, radiation_model.leaf_transmittance,
    )

    return (;
        layer_plant_area_index_boundaries, layer_plant_area_index_above,
        canopy_projection_ratio,
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

    PAI = sum(plant_area_index)
    θ = zenith_angle
    ρg = max(ground_reflectance, 1.0e-4) # the two-stream solution divides by ground reflectance
    Idir = direct_horizontal_irradiance
    Idiff = diffuse_horizontal_irradiance
    Iglobal = Idir + Idiff
    F_d = boundary_downward_shortwave
    F_u = boundary_upward_shortwave

    n_layers = length(layer_plant_area_index_above)

    fill!(downward_direct_shortwave, 0.0u"W/m^2")
    fill!(downward_diffuse_shortwave, 0.0u"W/m^2")
    fill!(upward_diffuse_shortwave, 0.0u"W/m^2")
    fill!(sunlit_fraction, 0.0)
    fill!(absorbed_shortwave, 0.0u"W/m^2")
    fill!(F_d, 0.0u"W/m^2")
    fill!(F_u, 0.0u"W/m^2")

    # No shortwave radiation reaches the canopy
    if Iglobal <= 0.0u"W/m^2" || θ >= 90.0u"°"
        return (;
            ground_absorbed_shortwave = 0.0u"W/m^2",
            net_absorbed_shortwave = 0.0u"W/m^2",
        )
    end

    # two-stream optical coefficients
    leaf = buffers.leaf_optics
    diffuse = diffuse_two_stream_optics(leaf, ρg)

    kᵦ = ellipsoidal_extinction_coefficient(θ, buffers.canopy_projection_ratio)
    has_direct_beam = Idir > 0.0u"W/m^2"

    # convert direct horizontal irradiance to irradiance normal to the beam
    beam_normal_irradiance = has_direct_beam ? Idir / cos(θ) : 0.0u"W/m^2"

    ρmax = max(ρg, radiation_model.leaf_reflectance)
    direct = direct_two_stream_optics(PAI, kᵦ, ρg, leaf, diffuse)

    # total downward and upward shortwave fluxes at every layer
    for i in eachindex(layer_plant_area_index_boundaries)
        fractions = two_stream_flux_fractions(leaf, diffuse, direct, kᵦ, ρmax, has_direct_beam, layer_plant_area_index_boundaries[i])
        irradiances = two_stream_irradiances(fractions, Idir, Idiff)
        F_d[i] = irradiances.downward
        F_u[i] = irradiances.upward
    end

    for i in 1:n_layers
        # shortwave absorption is the flux convergence across layer i
        absorbed_shortwave[i] = (F_d[i] - F_d[i + 1]) + (F_u[i + 1] - F_u[i])
        fractions = two_stream_flux_fractions(leaf, diffuse, direct, kᵦ, ρmax, has_direct_beam, layer_plant_area_index_above[i])
        downward_diffuse_shortwave[i] = fractions.diffuse_down * Idiff + fractions.direct_down * Idir
        upward_diffuse_shortwave[i] = fractions.diffuse_up * Idiff + fractions.direct_up * Idir
        # unscattered direct-beam irradiance normal to the solar beam
        downward_direct_shortwave[i] = beam_normal_irradiance * fractions.direct_transmission
        # direct-beam transmission is also the modelled sunlit fraction at this canopy depth.
        sunlit_fraction[i] = fractions.direct_transmission
    end

    ground_absorbed_shortwave = (1.0 - ρg) * F_d[n_layers + 1]
    net_absorbed_shortwave = Iglobal - F_u[1]

    return (; ground_absorbed_shortwave, net_absorbed_shortwave)
end
