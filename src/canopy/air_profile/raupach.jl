"""
    RaupachLTheoryAirProfile(; ground_velocity_std_factor=0.25, canopy_top_velocity_std_factor=1.25,
                              min_ground_resistance=2.0u"s/m", relaxation=0.5,
                              aitken_omega_min=0.02, aitken_omega_max=0.9,
                              aitken_weight_bottom=0.05, aitken_weight_top=0.8,
                              aitken_bottom_emphasis=10.0, near_field_subdivisions=20)

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
ground-most layer's own aerodynamic resistance, used once to evaluate the
ground flux Fg (avoids a runaway Fg for a very open canopy).

The near-field kernel is singular at zero distance, so a layer's contribution
to a point within (or coincident with) itself can't be evaluated as a single
point source the way every other (i,j) pair is. Instead that self term is
resolved by subdividing the source layer's own thickness into
`near_field_subdivisions` sub-points and averaging the kernel over them — a
plain midpoint quadrature of the (integrable) singularity, not a fitted
correction. Raising `near_field_subdivisions` trades runtime for accuracy;
the O(n²) near-field double sum only needs it for the O(n) diagonal terms; the
off-diagonal terms use a single point-source evaluation each.

Fg (ground flux) depends on the same air_temperature it helps solve,
coupled through up to the whole canopy's resistance — gain through that
loop can exceed 1 for a tall/dense canopy, which no fixed relaxation
factor stabilizes. Adaptive Aitken Δ²-acceleration (ports R micropoint's
`aitkin_weightdif`) handles this: `omega` is relearned each pass from the
residual history, clamped to `[aitken_omega_min, aitken_omega_max]`, and
weighted to under-relax hardest at the ground (`aitken_weight_bottom`)
and lightest at canopy top (`aitken_weight_top`); `aitken_bottom_emphasis`
weights ground-proximate residuals more heavily in the omega fit.
`relaxation` seeds the first pass, before there's a residual history.

All keyword defaults above are free/tunable, not literature-derived
values — no citation for the specific magnitudes.

Humidity is transported as vapor density (kg/m³), not vapor pressure: a
density gradient is directly a mass flux (Fick's law), the same scalar
`KTheoryAirProfile`'s own vapor solve uses, with no scale factor analogous
to heat's `ρ_cp` needed. A mixing-ratio-style concentration (`e·ρ_air·L/P`)
would need a standard 0.622 factor that isn't easily verified against a
closed form; vapor density sidesteps that conversion entirely. Boundary
conditions and the relative-humidity readout go through
`wet_air_properties`'s vapor-pressure/vapor-density primitives, the same
ones `monin_obukhov.jl` uses.

`far_field_mode` selects the far-field/ground-flux formula. `Val(:exact)`
(default) is eq. 19a's genuine integral: each layer's own local flux over
its own local resistance, summed (O(n) prefix sum), with a single ground
flux evaluated once from the ground-most layer. `Val(:bulk)` matches R
micropoint's `LangrangianOne` (and its C++ port): the flux *at* the query
layer times the total resistance from that layer to the top (or to the
ground, for the ground flux), recomputed per layer — structurally the same
flux-times-total-resistance shortcut :exact avoids, not the per-layer
integral — plus a near-field small-sample-size correction (`mu`, a function
of layer count alone) that :exact omits. Neither is strictly "more correct":
:exact is the literature-faithful form, :bulk trades that for matching R's
actual reference output, useful when comparing against data or trading
speed for realism (:bulk is a bit cheaper — no prefix *and* suffix sum of
resistance, just the two combined). The near-field double sum stays O(n²)
regardless of mode — its kernel is nonlinear in both indices — but heat and
vapor share one kernel evaluation per (i,j) pair. Out-of-range results are
asserted finite rather than silently clipped.

# References
- Raupach, M. R. (1989). A practical Lagrangian method for relating scalar
  concentrations to source distributions in vegetation canopies. *Quarterly
  Journal of the Royal Meteorological Society*, 115(487), 609-632.
"""
@kwdef struct RaupachLTheoryAirProfile{GVSF,CVSF,GR,RL,OMIN,OMAX,WBOT,WTOP,BETA,NS,FFM} <: AbstractCanopyAirProfileModel
    ground_velocity_std_factor::GVSF = 0.25
    canopy_top_velocity_std_factor::CVSF = 1.25
    min_ground_resistance::GR = 2.0u"s/m"
    relaxation::RL = 0.5
    aitken_omega_min::OMIN = 0.02
    aitken_omega_max::OMAX = 0.9
    aitken_weight_bottom::WBOT = 0.05
    aitken_weight_top::WTOP = 0.8
    aitken_bottom_emphasis::BETA = 10.0
    near_field_subdivisions::NS = 20
    far_field_mode::FFM = Val(:exact)
end

function allocate_air_profile(model::RaupachLTheoryAirProfile, canopy_height, plant_area_index, heights, n_layers, boundary_layer_model)
    (; layer_heights, layer_thickness) = canopy_layer_heights(heights, canopy_height, n_layers)
    initial_omega = clamp(Float64(model.relaxation), model.aitken_omega_min, model.aitken_omega_max)

    return (;
        layer_heights, layer_thickness,
        vertical_velocity_std = zeros(typeof(0.0u"m/s"), n_layers),
        inv_near_field_length = zeros(typeof(1.0 / (1.0u"m/s" * 1.0u"s")), n_layers),
        eddy_diffusivity = zeros(typeof(0.0u"m^2/s"), n_layers),
        layer_resistance = zeros(typeof(0.0u"s/m"), n_layers),
        far_field_accum = zeros(typeof(0.0u"J/m^3"), n_layers),
        far_field_accum_latent = zeros(typeof(0.0u"kg/m^3"), n_layers),
        # Only used by far_field_mode==:bulk; allocated unconditionally to keep
        # the buffer NamedTuple's shape (and canopy_air_profile!'s type) fixed
        # across mode choices.
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
        # Aitken state persists across Picard passes and hours (mirrors R
        # micropoint's WAitkenState) -- omega warm-starts hour to hour.
        aitken_omega = Ref(initial_omega),
        aitken_omega_latent = Ref(initial_omega),
        aitken_have_prev = Ref(false),
        aitken_have_prev_latent = Ref(false),
        aitken_residual_prev = zeros(typeof(0.0u"K"), n_layers),
        aitken_residual_prev_latent = zeros(typeof(0.0u"kg/m^3"), n_layers),
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

# kn is odd (kn(-ζ) = -kn(ζ); eq. 34 gives it for ζ>0 only).
@inline _raupach_kernel_signed(ζ) = ζ >= zero(ζ) ? _raupach_kernel(ζ) : -_raupach_kernel(-ζ)

# Eq. 37: Cn = ∫ [S(z0)/σw(z0)] · [kn{(z-z0)/λ} + kn{(z+z0)/λ}] dz0 — direct-path
# and reflected (ground image-source) terms each get their own kernel evaluation.
# `source_idx` indexes into `layer_heights`/`layer_thickness`/`inv_near_field_length`.
@inline function _raupach_kernel_weight(eval_height, source_idx, layer_heights, inv_near_field_length)
    inv_len = inv_near_field_length[source_idx]
    signed_direct = (eval_height - layer_heights[source_idx]) * inv_len
    reflected = (eval_height + layer_heights[source_idx]) * inv_len
    return _raupach_kernel_signed(signed_direct) + _raupach_kernel_signed(reflected)
end

# Self-term: the kernel is singular at zero distance, so a source layer's
# contribution to a point within its own thickness can't use a single
# point-source evaluation. Subdivides the layer into `subdivisions` sub-points
# and averages — a plain midpoint quadrature of the singularity. Only the
# direct-path term is ever near-singular; the reflected term isn't.
@inline function _raupach_self_kernel_weight(eval_height, source_idx, layer_heights, layer_thickness, inv_near_field_length, subdivisions)
    Δz = layer_thickness[source_idx]
    z0 = layer_heights[source_idx]
    inv_len = inv_near_field_length[source_idx]
    total = 0.0
    @inbounds for k in 1:subdivisions
        z_sub = z0 - Δz / 2 + Δz * (k - 0.5) / subdivisions
        signed_direct = (eval_height - z_sub) * inv_len
        direct_sign = signed_direct >= zero(signed_direct) ? 1.0 : -1.0
        ζ = max(abs(signed_direct), 1.0e-9)  # numerical floor, not a physical parameter
        reflected = (eval_height + z_sub) * inv_len
        total += direct_sign * _raupach_kernel(ζ) + _raupach_kernel_signed(reflected)
    end
    return total / subdivisions
end

# Weighted Aitken Δ²-acceleration (ports R micropoint's aitkin_weightdif).
# `newv` holds the raw unrelaxed solve on entry, the blended result on exit.
# `unit` strips units for the (dimensionless) omega fit only.
function _aitken_relax!(newv, oldv, layer_heights, canopy_height, omega_ref, have_prev_ref, residual_prev,
    omega_min, omega_max, weight_bottom, weight_top, bottom_emphasis, unit)
    n = length(newv)
    Δw = weight_top - weight_bottom
    if !have_prev_ref[]
        ω = clamp(omega_ref[], omega_min, omega_max)
        @inbounds for i in 1:n
            r = newv[i] - oldv[i]
            residual_prev[i] = r
            wz = weight_bottom + Δw * (layer_heights[i] / canopy_height)^2
            newv[i] = oldv[i] + (ω * wz) * r
        end
        omega_ref[] = ω
        have_prev_ref[] = true
        return nothing
    end
    num = 0.0
    den = 0.0
    @inbounds for i in 1:n
        r = newv[i] - oldv[i]
        dr = r - residual_prev[i]
        t = 1.0 - layer_heights[i] / canopy_height
        g = 1.0 + bottom_emphasis * t^2
        num += g * ustrip(unit, residual_prev[i]) * ustrip(unit, dr)
        den += g * ustrip(unit, dr)^2
    end
    ω = den > 0.0 ? -omega_ref[] * (num / den) : omega_ref[]
    ω = clamp(ω, omega_min, omega_max)
    @inbounds for i in 1:n
        r = newv[i] - oldv[i]
        residual_prev[i] = r
        wz = weight_bottom + Δw * (layer_heights[i] / canopy_height)^2
        newv[i] = oldv[i] + (ω * wz) * r
    end
    omega_ref[] = ω
    return nothing
end

# Eq. 19a: Cf(z)-Cf(zR) = ∫[z,zR] F(z')/Kf(z') dz' -- each layer's own local
# flux over its own local resistance, summed; not one flux value times the
# total resistance across the range. Fg is a single ground flux, evaluated
# once from the ground-most layer.
function _raupach_far_field!(::Val{:exact}, far_field_accum, far_field_accum_latent, _resistance_from_top, _resistance_to_ground,
    cumulative_sensible_below, cumulative_latent_below, layer_resistance,
    ground_temperature, ground_vapour_density, air_temperature_prev, vapour_density_prev, min_ground_resistance, n)
    ground_resistance = max(layer_resistance[n], min_ground_resistance)
    ground_flux = (calc_ρ_cp(air_temperature_prev[n]) / ground_resistance) * (ground_temperature - air_temperature_prev[n])
    ground_vapour_flux = (ground_vapour_density - vapour_density_prev[n]) / ground_resistance
    far_field_accum[1] = (cumulative_sensible_below[1] + ground_flux) * layer_resistance[1]
    far_field_accum_latent[1] = (cumulative_latent_below[1] + ground_vapour_flux) * layer_resistance[1]
    @inbounds for i in 2:n
        far_field_accum[i] = far_field_accum[i - 1] + (cumulative_sensible_below[i] + ground_flux) * layer_resistance[i]
        far_field_accum_latent[i] = far_field_accum_latent[i - 1] + (cumulative_latent_below[i] + ground_vapour_flux) * layer_resistance[i]
    end
    return nothing
end

# Matches R micropoint's LangrangianOne (and its C++ port): both the ground
# flux and the far-field term use the flux value at the query layer times
# the total resistance from that layer to the top/ground, recomputed per
# layer -- not eq. 19a's per-layer-weighted integral.
function _raupach_far_field!(::Val{:bulk}, far_field_accum, far_field_accum_latent, resistance_from_top, resistance_to_ground,
    cumulative_sensible_below, cumulative_latent_below, layer_resistance,
    ground_temperature, ground_vapour_density, air_temperature_prev, vapour_density_prev, min_ground_resistance, n)
    resistance_from_top[1] = layer_resistance[1]
    @inbounds for i in 2:n
        resistance_from_top[i] = resistance_from_top[i - 1] + layer_resistance[i]
    end
    resistance_to_ground[n] = layer_resistance[n]
    @inbounds for i in (n - 1):-1:1
        resistance_to_ground[i] = resistance_to_ground[i + 1] + layer_resistance[i]
    end
    @inbounds for i in 1:n
        ground_resistance = max(resistance_to_ground[i], min_ground_resistance)
        ground_flux = (calc_ρ_cp(air_temperature_prev[i]) / ground_resistance) * (ground_temperature - air_temperature_prev[i])
        ground_vapour_flux = (ground_vapour_density - vapour_density_prev[i]) / ground_resistance
        far_field_accum[i] = (cumulative_sensible_below[i] + ground_flux) * resistance_from_top[i]
        far_field_accum_latent[i] = (cumulative_latent_below[i] + ground_vapour_flux) * resistance_from_top[i]
    end
    return nothing
end

# R micropoint's near-field small-sample-size correction, a function of layer
# count alone; :exact stays eq.-19a/37-faithful and omits it.
@inline _raupach_near_field_scale(::Val{:exact}, n) = 1.0
@inline _raupach_near_field_scale(::Val{:bulk}, n) = 1.0 + 0.894 * exp(-0.01386 * n) + 9.82 * exp(-0.15 * n)

function canopy_air_profile!(buffers, model::RaupachLTheoryAirProfile, boundary_layer_model;
    canopy_height, displacement_height, friction_velocity,
    canopy_top_air_temperature, canopy_top_relative_humidity, ground_temperature, ground_relative_humidity,
    sensible_heat_source, evaporation_mass_flow, obukhov_length, atmospheric_pressure, vapour_pressure_equation=GoffGratch(),
)
    (; layer_heights, layer_thickness, vertical_velocity_std, inv_near_field_length, eddy_diffusivity,
       layer_resistance, far_field_accum, far_field_accum_latent, resistance_from_top, resistance_to_ground,
       cumulative_sensible_below, sensible_near_field_weight, cumulative_latent_below, latent_near_field_weight,
       air_temperature, air_temperature_prev, vapour_density, vapour_density_prev, relative_humidity,
       aitken_omega, aitken_omega_latent, aitken_have_prev, aitken_have_prev_latent,
       aitken_residual_prev, aitken_residual_prev_latent) = buffers
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

    (; ground_velocity_std_factor, canopy_top_velocity_std_factor, min_ground_resistance,
       aitken_omega_min, aitken_omega_max, aitken_weight_bottom, aitken_weight_top, aitken_bottom_emphasis,
       near_field_subdivisions) = model
    a0, a1 = ground_velocity_std_factor, canopy_top_velocity_std_factor  # Raupach's own notation, aliased for the formulas below
    γ = boundary_layer_model.dyer_constant

    # T_L is held constant through the canopy (a2·h/u*), not a function of z:
    # a2 is set by requiring eddy diffusivity to be continuous at the canopy
    # top, matching the below-canopy L-theory value σ_w(h)²·T_L against the
    # above-canopy value, using σ_w(h) = a1·u*. Φ_h multiplies the neutral
    # value (Ogée et al. 2003 eq. 8), not divides.
    z_eval = max(canopy_height - displacement_height, 1.0e-3u"m")  # numerical floor, not a physical parameter
    # obukhov_length is Inf on the stable/neutral branch; calc_Φ_h assumes unstable-only input.
    Φ_h = isfinite(obukhov_length) ? calc_Φ_h(z_eval, γ, obukhov_length) : 1.0
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

    cumulative_sensible_below[n] = sensible_heat_source[n]
    cumulative_latent_below[n] = evaporation_mass_flow[n] / 1.0u"m^2"
    @inbounds for i in (n - 1):-1:1
        cumulative_sensible_below[i] = cumulative_sensible_below[i + 1] + sensible_heat_source[i]
        cumulative_latent_below[i] = cumulative_latent_below[i + 1] + evaporation_mass_flow[i] / 1.0u"m^2"
    end

    ground_vapour_density = wet_air_properties(ground_temperature, ground_relative_humidity,
        atmospheric_pressure; vapour_pressure_equation).vapour_density

    _raupach_far_field!(model.far_field_mode, far_field_accum, far_field_accum_latent, resistance_from_top, resistance_to_ground,
        cumulative_sensible_below, cumulative_latent_below, layer_resistance,
        ground_temperature, ground_vapour_density, air_temperature_prev, vapour_density_prev, min_ground_resistance, n)
    near_field_scale = _raupach_near_field_scale(model.far_field_mode, n)

    ρ_cp_top = calc_ρ_cp(canopy_top_air_temperature)
    concentration_top = ρ_cp_top * canopy_top_air_temperature
    concentration_top_latent = wet_air_properties(canopy_top_air_temperature, canopy_top_relative_humidity,
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
    near_field_top *= near_field_scale
    near_field_top_latent *= near_field_scale

    @inbounds for i in 1:n
        far_field = far_field_accum[i]
        far_field_latent = far_field_accum_latent[i]

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
        near_field *= near_field_scale
        near_field_latent *= near_field_scale

        concentration = concentration_top - near_field_top + far_field + near_field
        air_temperature[i] = concentration / calc_ρ_cp(air_temperature_prev[i])  # raw; Aitken-blended below

        concentration_latent = concentration_top_latent - near_field_top_latent + far_field_latent + near_field_latent
        vapour_density[i] = concentration_latent  # raw; Aitken-blended below
    end

    # Aitken needs the whole residual vector to fit omega, so this is a
    # separate pass over the raw solve above, not folded into that loop.
    _aitken_relax!(air_temperature, air_temperature_prev, layer_heights, canopy_height,
        aitken_omega, aitken_have_prev, aitken_residual_prev,
        aitken_omega_min, aitken_omega_max, aitken_weight_bottom, aitken_weight_top, aitken_bottom_emphasis, u"K")
    _aitken_relax!(vapour_density, vapour_density_prev, layer_heights, canopy_height,
        aitken_omega_latent, aitken_have_prev_latent, aitken_residual_prev_latent,
        aitken_omega_min, aitken_omega_max, aitken_weight_bottom, aitken_weight_top, aitken_bottom_emphasis, u"kg/m^3")

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
