"""
    RaupachLTheoryAirProfile(; ground_velocity_std_factor=0.25, canopy_top_velocity_std_factor=1.25,
                              min_ground_resistance=2.0u"s/m", relaxation=0.5, near_field_subdivisions=20)

Raupach's (1989) localized near-field (LNF) Lagrangian in-canopy scalar
transport model: each layer's temperature and humidity is the sum of a
far-field contribution (a resistance-weighted integral of all sources, as if
fully mixed) and a near-field contribution (a direct near-source correction,
since K-theory's gradient-diffusion assumption breaks down close to a
source in a canopy). Unlike `KTheoryAirProfile`, transport is not local —
every layer's concentration depends on every other layer's source directly.

`ground_velocity_std_factor`/`canopy_top_velocity_std_factor` (Raupach's own
`a0`/`a1`) set the vertical velocity standard deviation profile: `σ_w(z) =
(a1+a0)/2·u* + (a1-a0)/2·u*·cos(π(1 - z/h))`, minimum at the ground
(`a0·u*`), maximum at canopy top (`a1·u*`). `min_ground_resistance` floors the
ground-to-layer aerodynamic resistance (avoids a runaway ground flux for a
very open canopy).

The near-field kernel is singular at zero distance, so a layer's contribution
to a point within (or coincident with) itself can't be evaluated as a single
point source the way every other (i,j) pair is. Instead that self term is
resolved by subdividing the source layer's own thickness into
`near_field_subdivisions` sub-points and averaging the kernel over them — a
plain midpoint quadrature of the (integrable) singularity, not a fitted
correction. Raising `near_field_subdivisions` trades runtime for accuracy;
the O(n²) near-field double sum only needs it for the O(n) diagonal terms; the
off-diagonal terms use a single point-source evaluation each.

`relaxation` under-relaxes each Picard pass against the previous pass's own
solution: `x_new = relaxation·x_solved + (1-relaxation)·x_prev`. Needed
because, unlike K-theory's direct tridiagonal solve, this Lagrangian solve is
not guaranteed monotonically convergent pass-to-pass.

Humidity is transported as vapor density (kg/m³), not vapor pressure: a
density gradient is directly a mass flux (Fick's law), the same scalar
`KTheoryAirProfile`'s own vapor solve uses, with no scale factor analogous
to heat's `ρ_cp` needed. A mixing-ratio-style concentration (`e·ρ_air·L/P`)
would need a standard 0.622 factor that isn't easily verified against a
closed form; vapor density sidesteps that conversion entirely. Boundary
conditions and the relative-humidity readout go through
`wet_air_properties`'s vapor-pressure/vapor-density primitives, the same
ones `monin_obukhov.jl` uses.

The per-layer cumulative source sum and far-field integral are each an
O(n) prefix/suffix sum, not recomputed per layer. The near-field double sum
stays O(n²) — its kernel is nonlinear in both indices — but heat and vapor
share one kernel evaluation per (i,j) pair. Out-of-range results are
asserted finite rather than silently clipped.

# References
- Raupach, M. R. (1989). A practical Lagrangian method for relating scalar
  concentrations to source distributions in vegetation canopies. *Quarterly
  Journal of the Royal Meteorological Society*, 115(487), 609-632.
"""
@kwdef struct RaupachLTheoryAirProfile{GVSF,CVSF,GR,RL,NS} <: AbstractCanopyAirProfileModel
    ground_velocity_std_factor::GVSF = 0.25
    canopy_top_velocity_std_factor::CVSF = 1.25
    min_ground_resistance::GR = 2.0u"s/m"
    relaxation::RL = 0.5
    near_field_subdivisions::NS = 20
end

function allocate_air_profile(::RaupachLTheoryAirProfile, canopy_height, plant_area_index, heights, n_layers, boundary_layer_model)
    layer_heights = sort(heights[heights .<= canopy_height]; rev=true)  # top-to-bottom, matching K-theory/radiation/wind convention
    length(layer_heights) == n_layers ||
        throw(ArgumentError("number of `heights` at or below `canopy_height` must equal n_layers"))

    min_spacing = 1.0e-3u"m"
    layer_thickness = similar(layer_heights)
    @inbounds for i in 1:(n_layers - 1)
        layer_thickness[i] = max(layer_heights[i] - layer_heights[i + 1], min_spacing)
    end
    layer_thickness[n_layers] = max(layer_heights[n_layers], min_spacing)  # bottom layer down to the ground (z=0)

    return (;
        layer_heights, layer_thickness,
        vertical_velocity_std = zeros(typeof(0.0u"m/s"), n_layers),
        inv_near_field_length = zeros(typeof(1.0 / (1.0u"m/s" * 1.0u"s")), n_layers),
        eddy_diffusivity = zeros(typeof(0.0u"m^2/s"), n_layers),
        layer_resistance = zeros(typeof(0.0u"s/m"), n_layers),
        resistance_from_top = zeros(typeof(0.0u"s/m"), n_layers),
        resistance_to_ground = zeros(typeof(0.0u"s/m"), n_layers),
        cumulative_sensible_below = zeros(typeof(0.0u"W/m^2"), n_layers),
        sensible_near_field_weight = zeros(typeof(0.0u"W/m^2" / 1.0u"m/s"), n_layers),
        cumulative_latent_below = zeros(typeof(0.0u"kg/m^2/s"), n_layers),
        latent_near_field_weight = zeros(typeof(1.0u"kg/m^2/s" / 1.0u"m/s"), n_layers),
        air_temperature = zeros(typeof(0.0u"K"), n_layers),
        air_temperature_prev = zeros(typeof(0.0u"K"), n_layers),
        vapour_density = zeros(typeof(0.0u"kg/m^3"), n_layers),
        vapour_density_prev = zeros(typeof(0.0u"kg/m^3"), n_layers),
        relative_humidity = zeros(n_layers),
    )
end

# Raupach's (1989) near-field kernel, eq. (34): kn(ζ) ≈ c1·ln(1-e^-ζ) + c2·e^-ζ,
# ζ>0 the distance normalized by σ_w·T_L. c1/c2 are the closed forms he gives
# so that ∫₀^∞ kn(ζ)dζ = 0.5 (his eq. 35), not independent fitted constants.
const _RAUPACH_KERNEL_C1 = -1 / sqrt(2π)
const _RAUPACH_KERNEL_C2 = 0.5 - π^2 / (6 * sqrt(2π))

@inline function _raupach_kernel(ζ)
    e = exp(-ζ)
    return _RAUPACH_KERNEL_C1 * log(1.0 - e) + _RAUPACH_KERNEL_C2 * e
end

# Kernel weight for a source layer's contribution to a point at `eval_height`,
# including the reflected (ground image-source) term. `source_idx` indexes into
# `layer_heights`/`layer_thickness`/`inv_near_field_length`.
@inline function _raupach_kernel_weight(eval_height, source_idx, layer_heights, inv_near_field_length)
    inv_len = inv_near_field_length[source_idx]
    signed_direct = (eval_height - layer_heights[source_idx]) * inv_len
    reflected = (eval_height + layer_heights[source_idx]) * inv_len
    return _raupach_kernel(abs(signed_direct)) * (signed_direct + reflected)
end

# Self-term: the kernel is singular at zero distance, so a source layer's
# contribution to a point within its own thickness can't use a single
# point-source evaluation. Subdivides the layer into `subdivisions` sub-points
# and averages — a plain midpoint quadrature of the singularity.
@inline function _raupach_self_kernel_weight(eval_height, source_idx, layer_heights, layer_thickness, inv_near_field_length, subdivisions)
    Δz = layer_thickness[source_idx]
    z0 = layer_heights[source_idx]
    inv_len = inv_near_field_length[source_idx]
    total = 0.0
    @inbounds for k in 1:subdivisions
        z_sub = z0 - Δz / 2 + Δz * (k - 0.5) / subdivisions
        signed_direct = (eval_height - z_sub) * inv_len
        ζ = max(abs(signed_direct), 1.0e-9)
        reflected = (eval_height + z_sub) * inv_len
        total += _raupach_kernel(ζ) * (signed_direct + reflected)
    end
    return total / subdivisions
end

function canopy_air_profile!(buffers, model::RaupachLTheoryAirProfile, boundary_layer_model;
    canopy_height, displacement_height, friction_velocity,
    canopy_top_air_temperature, canopy_top_relative_humidity, ground_temperature, ground_relative_humidity,
    sensible_heat_source, evaporation_mass_flow, obukhov_length, atmospheric_pressure, vapour_pressure_equation=GoffGratch(),
)
    (; layer_heights, layer_thickness, vertical_velocity_std, inv_near_field_length, eddy_diffusivity,
       layer_resistance, resistance_from_top, resistance_to_ground, cumulative_sensible_below,
       sensible_near_field_weight, cumulative_latent_below, latent_near_field_weight,
       air_temperature, air_temperature_prev, vapour_density, vapour_density_prev, relative_humidity) = buffers
    n = length(air_temperature)
    length(sensible_heat_source) == n || throw(ArgumentError("sensible_heat_source must have one entry per canopy layer"))
    length(evaporation_mass_flow) == n || throw(ArgumentError("evaporation_mass_flow must have one entry per canopy layer"))

    # Snapshot before overwriting: ground-flux terms below read last pass's
    # state, not the value being solved this pass (avoids self-reference).
    # Also the relaxation reference for this pass's under-relaxed update.
    air_temperature_prev .= air_temperature
    # Unlike air_temperature/relative_humidity, energy_balance! doesn't seed
    # vapour_density (model-specific state) — seed it from the boundary
    # condition on first use (still exactly zero, its allocate default).
    if all(iszero, vapour_density)
        fill!(vapour_density, wet_air_properties(canopy_top_air_temperature, canopy_top_relative_humidity,
            atmospheric_pressure; vapour_pressure_equation).vapour_density)
    end
    vapour_density_prev .= vapour_density

    (; ground_velocity_std_factor, canopy_top_velocity_std_factor, min_ground_resistance, relaxation, near_field_subdivisions) = model
    a0, a1 = ground_velocity_std_factor, canopy_top_velocity_std_factor  # Raupach's own notation, aliased for the formulas below
    γ = boundary_layer_model.dyer_constant

    # T_L is held constant through the canopy (a2·h/u*), not a function of z:
    # a2 is set by requiring eddy diffusivity to be continuous at the canopy
    # top, matching the below-canopy L-theory value σ_w(h)²·T_L against the
    # above-canopy Monin-Obukhov value κu*(h-d)/Φ_h — using σ_w(h) = a1·u*
    # (Ogée et al. 2003).
    z_eval = max(canopy_height - displacement_height, 1.0e-3u"m")
    Φ_h = calc_Φ_h(z_eval, γ, obukhov_length)
    a2 = Φ_h * boundary_layer_model.karman_constant * (1.0 - displacement_height / canopy_height) / a1^2
    TL = a2 * canopy_height / friction_velocity

    σ_w_mean = (a1 + a0) * 0.5 * friction_velocity
    σ_w_amplitude = (a1 - a0) * 0.5 * friction_velocity
    @inbounds for i in 1:n
        σ_w = σ_w_mean + σ_w_amplitude * cos(π * (1.0 - layer_heights[i] / canopy_height))
        vertical_velocity_std[i] = σ_w
        inv_near_field_length[i] = 1.0 / (σ_w * TL)
        eddy_diffusivity[i] = TL * σ_w^2
        layer_resistance[i] = layer_thickness[i] / eddy_diffusivity[i]
        sensible_near_field_weight[i] = sensible_heat_source[i] / σ_w
        latent_near_field_weight[i] = (evaporation_mass_flow[i] / 1.0u"m^2") / σ_w
    end

    # Prefix (canopy top -> i) and suffix (i -> ground) cumulative resistance.
    resistance_from_top[1] = layer_resistance[1]
    @inbounds for i in 2:n
        resistance_from_top[i] = resistance_from_top[i - 1] + layer_resistance[i]
    end
    resistance_to_ground[n] = layer_resistance[n]
    @inbounds for i in (n - 1):-1:1
        resistance_to_ground[i] = resistance_to_ground[i + 1] + layer_resistance[i]
    end

    cumulative_sensible_below[n] = sensible_heat_source[n]
    cumulative_latent_below[n] = evaporation_mass_flow[n] / 1.0u"m^2"
    @inbounds for i in (n - 1):-1:1
        cumulative_sensible_below[i] = cumulative_sensible_below[i + 1] + sensible_heat_source[i]
        cumulative_latent_below[i] = cumulative_latent_below[i + 1] + evaporation_mass_flow[i] / 1.0u"m^2"
    end

    ρ_cp_top = calc_ρ_cp(canopy_top_air_temperature)
    concentration_top = ρ_cp_top * canopy_top_air_temperature
    concentration_top_latent = wet_air_properties(canopy_top_air_temperature, canopy_top_relative_humidity,
        atmospheric_pressure; vapour_pressure_equation).vapour_density

    ground_vapour_density = wet_air_properties(ground_temperature, ground_relative_humidity,
        atmospheric_pressure; vapour_pressure_equation).vapour_density

    # Near-field concentration at the canopy top, subtracted from every
    # layer's own near-field term below. Layers whose height is (numerically)
    # coincident with canopy_height use the subdivided self kernel; every
    # other layer uses a single point-source evaluation. Heat and vapor share
    # one kernel evaluation per layer (same turbulence field).
    near_field_top = zero(concentration_top)
    near_field_top_latent = zero(concentration_top_latent)
    @inbounds for i in 1:n
        signed_direct = (canopy_height - layer_heights[i]) * inv_near_field_length[i]
        kernel_weight = if abs(signed_direct) > 1.0e-9
            _raupach_kernel_weight(canopy_height, i, layer_heights, inv_near_field_length)
        else
            _raupach_self_kernel_weight(canopy_height, i, layer_heights, layer_thickness, inv_near_field_length, near_field_subdivisions)
        end
        near_field_top += sensible_near_field_weight[i] * kernel_weight
        near_field_top_latent += latent_near_field_weight[i] * kernel_weight
    end

    @inbounds for i in 1:n
        # Heat and vapor share the ground-to-layer resistance (K_H = K_V assumed).
        ground_resistance = max(resistance_to_ground[i], min_ground_resistance)
        ground_flux = (calc_ρ_cp(air_temperature_prev[i]) / ground_resistance) * (ground_temperature - air_temperature_prev[i])
        far_field_source = cumulative_sensible_below[i] + ground_flux
        far_field = far_field_source * resistance_from_top[i]

        ground_vapour_flux = (ground_vapour_density - vapour_density_prev[i]) / ground_resistance
        far_field_source_latent = cumulative_latent_below[i] + ground_vapour_flux
        far_field_latent = far_field_source_latent * resistance_from_top[i]

        near_field = zero(concentration_top)
        near_field_latent = zero(concentration_top_latent)
        for j in 1:n
            kernel_weight = if j == i
                _raupach_self_kernel_weight(layer_heights[i], j, layer_heights, layer_thickness, inv_near_field_length, near_field_subdivisions)
            else
                _raupach_kernel_weight(layer_heights[i], j, layer_heights, inv_near_field_length)
            end
            near_field += sensible_near_field_weight[j] * kernel_weight
            near_field_latent += latent_near_field_weight[j] * kernel_weight
        end

        concentration = concentration_top - near_field_top + far_field + near_field
        new_air_temperature = concentration / calc_ρ_cp(air_temperature_prev[i])
        air_temperature[i] = relaxation * new_air_temperature + (1.0 - relaxation) * air_temperature_prev[i]

        concentration_latent = concentration_top_latent - near_field_top_latent + far_field_latent + near_field_latent
        vapour_density[i] = relaxation * concentration_latent + (1.0 - relaxation) * vapour_density_prev[i]
    end

    all(isfinite, ustrip.(u"K", air_temperature)) ||
        throw(DomainError(air_temperature, "RaupachLTheoryAirProfile produced a non-finite air temperature"))
    all(isfinite, ustrip.(u"kg/m^3", vapour_density)) ||
        throw(DomainError(vapour_density, "RaupachLTheoryAirProfile produced a non-finite vapor density"))

    @inbounds for i in 1:n
        ρv_sat = wet_air_properties(air_temperature[i], 1.0, atmospheric_pressure; vapour_pressure_equation).vapour_density
        relative_humidity[i] = clamp(ustrip(u"kg/m^3", vapour_density[i]) / ustrip(u"kg/m^3", ρv_sat), 0.0, 1.0)
    end

    return (; air_temperature, relative_humidity)
end
