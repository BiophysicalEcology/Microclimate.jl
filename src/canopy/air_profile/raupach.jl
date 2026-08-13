"""
    RaupachLTheoryAirProfile(; ground_velocity_std_factor=0.25, canopy_top_velocity_std_factor=1.25,
                              min_ground_resistance=2.0u"s/m", relaxation=0.5,
                              aitken_omega_min=0.02, aitken_omega_max=0.9,
                              aitken_weight_bottom=0.05, aitken_weight_top=0.8,
                              aitken_bottom_emphasis=10.0, near_field_subdivisions=20,
                              max_lagrangian_timescale=200.0u"s", max_air_temperature_deviation=40.0u"K",
                              bulk_temperature_margin=2.0u"K")

Raupach's (1989) localized near-field (LNF) Lagrangian in-canopy scalar
transport model: each layer's temperature and humidity is the sum of a
far-field contribution (a resistance-weighted integral of all sources, as if
fully mixed) and a near-field contribution (a direct near-source correction,
since K-theory's gradient-diffusion assumption breaks down close to a
source in a canopy). Unlike `KTheoryAirProfile`, transport is not local —
every layer's concentration depends on every other layer's source directly.

`ground_velocity_std_factor`/`canopy_top_velocity_std_factor` (Raupach's
`a0`/`a1`) set the vertical velocity standard deviation profile: `σ_w(z) =
(a1+a0)/2·u* + (a1-a0)/2·u*·cos(π(1 - z/h))`, minimum at the ground
(`a0·u*`), maximum at canopy top (`a1·u*`). `min_ground_resistance` floors the
ground-most layer's own aerodynamic resistance, used once to evaluate the
ground flux Fg (avoids a runaway Fg for a very open canopy).

The near-field kernel is singular at zero distanceso its self term is
resolved by subdividing the source layer's own thickness into
`near_field_subdivisions` sub-points and averaging the kernel over them.
Raising `near_field_subdivisions` trades runtime for accuracy.

The recursive dependence of Fg (ground flux) on the air_temperature it helps 
solve and the potential runaway that causes for a tall/dense canopy is handled
with adaptive Aitken Δ²-acceleration: `omega` s relearned each pass from
the residual history, clamped to `[aitken_omega_min, aitken_omega_max]`, and 
weighted to under-relax hardest at the ground (`aitken_weight_bottom`) and 
lightest at canopy top (`aitken_weight_top`); `aitken_bottom_emphasis` weights 
ground-proximate residuals more heavily in the omega fit. `relaxation` seeds 
the first pass, before there's a residual history.

`T_L` (see `canopy_air_profile!`'s own derivation) is `∝ 1/friction_velocity`
and grows further under instability via `calc_Φ_h`, which can push `T_L`
(and the near-field kernel's `ζ = distance/(σ_w·T_L)` argument) toward
physically implausible eddy-turnover-time values under calm/unstable
conditions. `max_lagrangian_timescale` caps it; 200s is generous
relative to the seconds-to-tens-of-seconds range typical conditions produce.

The keyword defaults are free/tunable, not literature-derived values
— other than ground_velocity_std_factor, canopy_top_velocity_std_factor
which come from Raupach (1989).

`far_field_mode` selects between two self-consistent formulations of both
the far-field and near-field terms. `Val(:exact)` (default) follows Raupach 
(1989) exactly, using eq. 19a's integral, and near-field is eq. 37's 
separately-signed kernel evaluation with a coincident self term resolved
by subdivision quadrature. `Val(:bulk)` is an alternate formulation used in 
the micropoint model: far-field is the flux *at* the query layer times the total
resistance from that layer to the top (or to the ground, for the ground
flux), recomputed per layer, not the per-layer integral; near-field
evaluates the kernel once (at the direct distance) and multiplies by the
raw signed sum of direct+reflected normalized distances, and drops any self
term outright rather than resolving it; both terms also pick up a
near-field small-sample-size correction (`mu`, a function of layer count
alone) that :exact omits.

# References
- Raupach, M. R. (1989). A practical Lagrangian method for relating scalar
  concentrations to source distributions in vegetation canopies. *Quarterly
  Journal of the Royal Meteorological Society*, 115(487), 609-632.
"""
@kwdef struct RaupachLTheoryAirProfile{GVSF,CVSF,GR,RL,OMIN,OMAX,WBOT,WTOP,BETA,NS,FFM,MTL,MAD,BTM} <: AbstractCanopyAirProfileModel
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
    max_lagrangian_timescale::MTL = 200.0u"s"
    max_air_temperature_deviation::MAD = 40.0u"K"
    bulk_temperature_margin::BTM = 2.0u"K"
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
        vapor_density = zeros(typeof(0.0u"kg/m^3"), n_layers),
        vapor_density_prev = zeros(typeof(0.0u"kg/m^3"), n_layers),
        relative_humidity = zeros(n_layers),
        # Aitken state persists across Picard passes and hours -- omega
        # warm-starts hour to hour.
        aitken_omega = Ref(initial_omega),
        aitken_omega_latent = Ref(initial_omega),
        aitken_have_prev = Ref(false),
        aitken_have_prev_latent = Ref(false),
        aitken_residual_prev = zeros(typeof(0.0u"K"), n_layers),
        aitken_residual_prev_latent = zeros(typeof(0.0u"kg/m^3"), n_layers),
        ground_heat_conductance = Ref(0.0u"W/m^2/K"),
        ground_vapor_conductance = Ref(0.0u"m/s"),
    )
end

# Raupach's (1989) near-field kernel, eq. (34): kn(ζ) ≈ c1·ln(1-e^-ζ) + c2·e^-ζ,
# ζ>0 the distance normalized by σ_w·T_L. c1/c2 are the closed forms he gives
# so that ∫₀^∞ kn(ζ)dζ = 0.5 (his eq. 35), not independent fitted constants.
const _RAUPACH_KERNEL_C1 = -1 / sqrt(2π)
const _RAUPACH_KERNEL_C2 = 0.5 - π^2 / (6 * sqrt(2π))

@inline function _raupach_kernel(ζ)
    e = exp(-ζ)
    return _RAUPACH_KERNEL_C1 * log(1.0 - e) + _RAUPACH_KERNEL_C2 * e # eq. 34 Raupach 1989
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

# :bulk's near-field kernel evaluates kn once (at the direct distance only)
# and multiplies by the raw sum of direct+reflected normalized distances,
# not eq. 37's separately-signed evaluation -- the same shortcut as its
# far-field. It also drops any self term outright rather than resolving it
# by subdivision quadrature.
@inline function _raupach_bulk_kernel_weight(eval_height, source_idx, layer_heights, inv_near_field_length)
    inv_len = inv_near_field_length[source_idx]
    direct = (eval_height - layer_heights[source_idx]) * inv_len
    reflected = (eval_height + layer_heights[source_idx]) * inv_len
    return _raupach_kernel(abs(direct)) * (direct + reflected)
end

# Top-boundary near-field weight for source layer i, dispatched on
# far_field_mode. Layer 1 (closest to canopy_height) is always the self term
# -- matching _raupach_pair_kernel_weight's own j==i criterion, not a
# height-coincidence check, since layer_heights[1] need not sit exactly at
# canopy_height and the two criteria disagreeing double-counts/drops the
# singular kernel inconsistently between near_field_top and near_field[1].
@inline function _raupach_top_kernel_weight(::Val{:exact}, i, canopy_height, layer_heights, layer_thickness, inv_near_field_length, subdivisions)
    return i == 1 ?
        _raupach_self_kernel_weight(canopy_height, i, layer_heights, layer_thickness, inv_near_field_length, subdivisions) :
        _raupach_kernel_weight(canopy_height, i, layer_heights, inv_near_field_length)
end
@inline function _raupach_top_kernel_weight(::Val{:bulk}, i, canopy_height, layer_heights, layer_thickness, inv_near_field_length, subdivisions)
    return i == 1 ? 0.0 : _raupach_bulk_kernel_weight(canopy_height, i, layer_heights, inv_near_field_length)
end

# Near-field weight for one (eval layer i, source layer j) pair.
@inline function _raupach_pair_kernel_weight(::Val{:exact}, i, j, layer_heights, layer_thickness, inv_near_field_length, subdivisions)
    return j == i ?
        _raupach_self_kernel_weight(layer_heights[i], j, layer_heights, layer_thickness, inv_near_field_length, subdivisions) :
        _raupach_kernel_weight(layer_heights[i], j, layer_heights, inv_near_field_length)
end
@inline function _raupach_pair_kernel_weight(::Val{:bulk}, i, j, layer_heights, layer_thickness, inv_near_field_length, subdivisions)
    return j == i ? 0.0 : _raupach_bulk_kernel_weight(layer_heights[i], j, layer_heights, inv_near_field_length)
end

# Weighted Aitken Δ²-acceleration.
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
    cumulative_sensible_below, cumulative_latent_below, layer_resistance, _layer_thickness,
    ground_temperature, ground_vapor_density, air_temperature_prev, vapor_density_prev, ground_resistance, ρ_cp, n)
    # Same value as ground_heat_conductance/ground_vapor_conductance -- keeps soil/canopy flux consistent.
    ground_flux = (ρ_cp / ground_resistance) * (ground_temperature - air_temperature_prev[n])
    ground_vapor_flux = (ground_vapor_density - vapor_density_prev[n]) / ground_resistance
    far_field_accum[1] = (cumulative_sensible_below[1] + ground_flux) * layer_resistance[1]
    far_field_accum_latent[1] = (cumulative_latent_below[1] + ground_vapor_flux) * layer_resistance[1]
    @inbounds for i in 2:n
        far_field_accum[i] = far_field_accum[i - 1] + (cumulative_sensible_below[i] + ground_flux) * layer_resistance[i]
        far_field_accum_latent[i] = far_field_accum_latent[i - 1] + (cumulative_latent_below[i] + ground_vapor_flux) * layer_resistance[i]
    end
    return nothing
end

# :bulk's far-field: both the ground flux and the far-field term use the
# flux value at the query layer times the total resistance from that layer
# to the top/ground, recomputed per layer -- not eq. 19a's per-layer-weighted
# integral. Ground flux's ρ_cp is recomputed per layer from air_temperature_prev
# (unlike :exact's fixed ρ_cp), matching micropoint's LangrangianOne (ph/cp
# from tair[i]) exactly -- tried the fixed-ρ_cp form here too (removing the
# ρ_cp∝1/T self-reference), but it strengthens ground-temperature coupling
# into the outer canopy-top energy-balance loop and destabilizes it in
# practice, so this stays matched to R's own (empirically stable) form.
function _raupach_far_field!(::Val{:bulk}, far_field_accum, far_field_accum_latent, resistance_from_top, resistance_to_ground,
    cumulative_sensible_below, cumulative_latent_below, layer_resistance, _layer_thickness,
    ground_temperature, ground_vapor_density, air_temperature_prev, vapor_density_prev, ground_resistance, _ρ_cp, n)
    resistance_from_top[1] = layer_resistance[1]
    @inbounds for i in 2:n
        resistance_from_top[i] = resistance_from_top[i - 1] + layer_resistance[i]
    end
    # resistance_to_ground[i<n] adds on top of this, so already ≥ ground_resistance -- no extra floor needed.
    resistance_to_ground[n] = ground_resistance
    @inbounds for i in (n - 1):-1:1
        resistance_to_ground[i] = resistance_to_ground[i + 1] + layer_resistance[i]
    end
    @inbounds for i in 1:n
        ground_flux = (calc_ρ_cp(air_temperature_prev[i]) / resistance_to_ground[i]) * (ground_temperature - air_temperature_prev[i])
        ground_vapor_flux = (ground_vapor_density - vapor_density_prev[i]) / resistance_to_ground[i]
        far_field_accum[i] = (cumulative_sensible_below[i] + ground_flux) * resistance_from_top[i]
        far_field_accum_latent[i] = (cumulative_latent_below[i] + ground_vapor_flux) * resistance_from_top[i]
    end
    return nothing
end

# Near-field small-sample-size correction for :bulk, a function of layer
# count alone; :exact stays eq.-19a/37-faithful and omits it.
@inline _raupach_near_field_scale(::Val{:exact}, n) = 1.0
@inline _raupach_near_field_scale(::Val{:bulk}, n) = 1.0 + 0.894 * exp(-0.01386 * n) + 9.82 * exp(-0.15 * n)

# :exact uses one ρ_cp fixed from ambient data; :bulk recomputes per layer
# from air_temperature_prev[i], matching micropoint's LangrangianOne (see
# _raupach_far_field!'s own comment on why :bulk keeps this).
@inline _raupach_layer_ρ_cp(::Val{:exact}, ρ_cp, _air_temperature_prev, _i) = ρ_cp
@inline _raupach_layer_ρ_cp(::Val{:bulk}, _ρ_cp, air_temperature_prev, i) = calc_ρ_cp(air_temperature_prev[i])

# :exact's backstop is ±max_air_temperature_deviation around T_top/ground_temperature.
# :bulk matches micropoint's tmn/tmx exactly: ±bulk_temperature_margin around
# T_top/ground_temperature, widened (not margined) to cover leaf_temperature's own range.
@inline _raupach_temperature_bounds(::Val{:exact}, T_top, ground_temperature, leaf_temperature, max_air_temperature_deviation, bulk_temperature_margin) =
    (min(T_top, ground_temperature) - max_air_temperature_deviation, max(T_top, ground_temperature) + max_air_temperature_deviation)
function _raupach_temperature_bounds(::Val{:bulk}, T_top, ground_temperature, leaf_temperature, max_air_temperature_deviation, bulk_temperature_margin)
    leaf_lo, leaf_hi = extrema(leaf_temperature)
    bound_lo = min(min(T_top, ground_temperature) - bulk_temperature_margin, leaf_lo)
    bound_hi = max(max(T_top, ground_temperature) + bulk_temperature_margin, leaf_hi)
    return (bound_lo, bound_hi)
end

function canopy_air_profile!(buffers, model::RaupachLTheoryAirProfile, boundary_layer_model;
    canopy_height, displacement_height, friction_velocity, wind_attenuation,
    canopy_top_air_temperature, canopy_top_relative_humidity, ground_temperature, ground_relative_humidity,
    sensible_heat_source, evaporation_mass_flow, obukhov_length, atmospheric_pressure, vapour_pressure_equation=GoffGratch(),
    leaf_temperature=nothing,
)
    (; layer_heights, layer_thickness, vertical_velocity_std, inv_near_field_length, eddy_diffusivity,
       layer_resistance, far_field_accum, far_field_accum_latent, resistance_from_top, resistance_to_ground,
       cumulative_sensible_below, sensible_near_field_weight, cumulative_latent_below, latent_near_field_weight,
       air_temperature, air_temperature_prev, vapor_density, vapor_density_prev, relative_humidity,
       aitken_omega, aitken_omega_latent, aitken_have_prev, aitken_have_prev_latent,
       aitken_residual_prev, aitken_residual_prev_latent,
       ground_heat_conductance, ground_vapor_conductance) = buffers
    n = length(air_temperature)
    length(sensible_heat_source) == n || throw(ArgumentError("sensible_heat_source must have one entry per canopy layer"))
    length(evaporation_mass_flow) == n || throw(ArgumentError("evaporation_mass_flow must have one entry per canopy layer"))

    T_top = canopy_top_air_temperature
    ρv_top = wet_air_properties(T_top, canopy_top_relative_humidity, atmospheric_pressure; vapour_pressure_equation).vapour_density

    # Snapshot before overwriting: ground-flux terms below read last pass's
    # state, not the value being solved this pass (avoids self-reference).
    # Also the relaxation reference for this pass's under-relaxed update.
    air_temperature_prev .= air_temperature
    # Unlike air_temperature/relative_humidity, energy_balance! doesn't seed
    # vapor_density (model-specific state) — seed it from the boundary
    # condition on first use (still exactly zero, its allocate default).
    if all(iszero, vapor_density)
        fill!(vapor_density, ρv_top)
    end
    vapor_density_prev .= vapor_density

    (; ground_velocity_std_factor, canopy_top_velocity_std_factor, min_ground_resistance,
       aitken_omega_min, aitken_omega_max, aitken_weight_bottom, aitken_weight_top, aitken_bottom_emphasis,
       near_field_subdivisions, max_lagrangian_timescale, max_air_temperature_deviation) = model
    a0, a1 = ground_velocity_std_factor, canopy_top_velocity_std_factor  # Raupach's own notation, aliased for the formulas below
    γ = boundary_layer_model.dyer_constant

# Eq. 18 (K_f = σ_w²·T_L) → eddy_diffusivity[i]. T_L is Raupach's Lagrangian integral timescale, computed once above this function 
# (not shown in your snippet) as T_L = min(a2 * canopy_height / friction_velocity, max_lagrangian_timescale), where a2 is calibrated 
# by matching the below-canopy diffusivity continuous with the above-canopy value at canopy top — this is Raupach's own closure for 
# T_L. The paper allows T_L(z') to vary with height inside the integral; the code (and R's LangrangianOne, which does the same: 
# const double TL = windvars.a2 * vegp.hgt / windvars.uf;) holds it constant across the whole canopy — a documented simplification,
# not a discretization artifact. That single T_L then feeds eddy_diffusivity[i] = T_L * σ_w[i]^2 — literally K_f(z_i) from eq. 18,
# with only σ_w varying by layer.
# Eq. 19a's integrand 1/(σ_w²(z')T_L(z')) → layer_resistance[i]. layer_resistance[i] = layer_thickness[i] / eddy_diffusivity[i] is
# Δz_i / K_f(z_i) — a rectangle-rule discretization of dz'/K_f(z') for layer i. This is the "local resistance" of that one layer.
# The bracket [∫₀^z' S(z'')dz'' + F_g] → cumulative_sensible_below[i] + ground_flux. cumulative_sensible_below[i] sums 
# sensible_heat_source[k] for every layer k from i down to the ground (layer n in this file's top-to-bottom indexing) — 
# that's the discretized ∫₀^{z'} S(z'')dz'', cumulative source strength from the ground up to z'. ground_flux is F_g, 
# computed once from the ground-most layer's own resistance/temperature difference, exactly matching the paper's F_g
# being a single fixed input to the integral, not itself a function of z' — which is why it's added identically at every 
# i rather than recomputed.
# The outer integral ∫_z^{z_R}{...}dz' → the prefix-sum recursion. far_field_accum[1] = (...)*layer_resistance[1] starts 
# at layer 1 = canopy top = z_R (the reference height where the boundary condition C_f(z_R) is known), and 
# far_field_accum[i] = far_field_accum[i-1] + (...)*layer_resistance[i] walks the integral downward layer by layer
# — a running Riemann sum of the integrand from z_R down to z. So far_field_accum[i] is C_f(z_i) − C_f(z_R), the 
# left-hand side of eq. 19a, directly.
# C_f(z_R) itself is T_top (expressed as a temperature rather than a concentration — the ρ_cp conversion happens where
# air_temperature[i] = T_top + far_field + ... is assembled, outside this function) — so the caller recovers
#  C_f(z) = C_f(z_R) + (C_f(z) − C_f(z_R)) exactly as eq. 19a's rearranged form.

    # T_L held constant through the canopy (a2·h/u*): a2 matches the
    # below-canopy value σ_w(h)²·T_L to the above-canopy MOST value
    # κu*(h-d)/Φ_h.
    z_eval = max(canopy_height - displacement_height, 1.0e-3u"m")  # numerical floor, not a physical parameter
    # obukhov_length is Inf only on the true neutral fallback now (calc_Φ_h handles both signs).
    Φ_h = isfinite(obukhov_length) ?
        calc_Φ_h(z_eval, γ, obukhov_length, boundary_layer_model.stable_Φ_h_coefficient,
            boundary_layer_model.min_stable_Φ_h, boundary_layer_model.max_stable_Φ_h) : 1.0
    a2 = boundary_layer_model.karman_constant * (1.0 - displacement_height / canopy_height) / (a1^2 * Φ_h)
    T_L = min(a2 * canopy_height / friction_velocity, max_lagrangian_timescale)

    σ_w_mean = (a1 + a0) * 0.5 * friction_velocity
    σ_w_amplitude = (a1 - a0) * 0.5 * friction_velocity
    @inbounds for i in 1:n
        σ_w = σ_w_mean + σ_w_amplitude * cos(π * (1.0 - layer_heights[i] / canopy_height)) # Raupach 1989 eq. 48
        vertical_velocity_std[i] = σ_w
        inv_near_field_length[i] = 1.0 / (σ_w * T_L)
        eddy_diffusivity[i] = T_L * σ_w^2 # Raupach 1989 eq. 11b
        layer_resistance[i] = layer_thickness[i] / eddy_diffusivity[i]
        sensible_near_field_weight[i] = sensible_heat_source[i] / σ_w
        latent_near_field_weight[i] = (evaporation_mass_flow[i] / 1.0u"m^2") / σ_w
    end

    # Shared by ground_heat_conductance/ground_vapor_conductance and _raupach_far_field!'s
    # ground flux -- uses wind_attenuation_profile's shape, not layer_resistance[n]'s σ_w profile.
    eddy_diffusivity_top_wind = boundary_layer_model.karman_constant * friction_velocity * z_eval
    eddy_diffusivity_ground = eddy_diffusivity_top_wind * wind_attenuation[n]
    ground_resistance = max(layer_thickness[n] / eddy_diffusivity_ground, min_ground_resistance)
    # ρ_cp fixed from T_top/ground_temperature, not from the (possibly still
    # drifting) air_temperature_prev iterate: calc_ρ_cp ∝ 1/T, so tying it to
    # each layer's own previous pass makes warming self-amplifying (ΔT ∝
    # flux·T_prev) -- same failure mode as canopy_top_flux_boundary's ρ_cp,
    # here spread across every layer/pass instead of one fixed-point loop.
    ρ_cp = calc_ρ_cp((T_top + ground_temperature) / 2)
    ground_heat_conductance[] = ρ_cp / ground_resistance
    ground_vapor_conductance[] = 1.0 / ground_resistance

    cumulative_sensible_below[n] = sensible_heat_source[n]
    cumulative_latent_below[n] = evaporation_mass_flow[n] / 1.0u"m^2"
    @inbounds for i in (n - 1):-1:1
        cumulative_sensible_below[i] = cumulative_sensible_below[i + 1] + sensible_heat_source[i]
        cumulative_latent_below[i] = cumulative_latent_below[i + 1] + evaporation_mass_flow[i] / 1.0u"m^2"
    end

    ρv_g = wet_air_properties(ground_temperature, ground_relative_humidity,
        atmospheric_pressure; vapour_pressure_equation).vapour_density

    _raupach_far_field!(model.far_field_mode, far_field_accum, far_field_accum_latent, resistance_from_top, resistance_to_ground,
        cumulative_sensible_below, cumulative_latent_below, layer_resistance, layer_thickness,
        ground_temperature, ρv_g, air_temperature_prev, vapor_density_prev, ground_resistance, ρ_cp, n)
    near_field_scale = _raupach_near_field_scale(model.far_field_mode, n)

    # Near-field concentration at the canopy top, subtracted from every
    # layer's own near-field term below. Layers whose height is (numerically)
    # coincident with canopy_height use the subdivided self kernel; every
    # other layer uses a single point-source evaluation. Heat and vapor share
    # one kernel evaluation per layer (same turbulence field).
    near_field_top = zero(eltype(far_field_accum))
    near_field_top_latent = zero(ρv_top)
    @inbounds for i in 1:n
        kernel_weight = _raupach_top_kernel_weight(model.far_field_mode, i, canopy_height, layer_heights, layer_thickness, inv_near_field_length, near_field_subdivisions)
        near_field_top += sensible_near_field_weight[i] * kernel_weight
        near_field_top_latent += latent_near_field_weight[i] * kernel_weight
    end
    near_field_top *= near_field_scale
    near_field_top_latent *= near_field_scale

    @inbounds for i in 1:n
        far_field = far_field_accum[i]
        far_field_latent = far_field_accum_latent[i]

        near_field = zero(near_field_top)
        near_field_latent = zero(ρv_top)
        for j in 1:n
            kernel_weight = _raupach_pair_kernel_weight(model.far_field_mode, i, j, layer_heights, layer_thickness, inv_near_field_length, near_field_subdivisions)
            near_field += sensible_near_field_weight[j] * kernel_weight
            near_field_latent += latent_near_field_weight[j] * kernel_weight
        end
        near_field *= near_field_scale
        near_field_latent *= near_field_scale

        # Heat: T_top is an absolute temperature baseline, not a concentration;
        # ρ_cp converts the source-driven perturbation (a genuine volumetric
        # energy density) into a temperature increment added on top of it.
        ρ_cp_i = _raupach_layer_ρ_cp(model.far_field_mode, ρ_cp, air_temperature_prev, i)
        air_temperature[i] = T_top + (far_field + near_field - near_field_top) / ρ_cp_i  # raw; Aitken-blended below

        # Latent: vapor density is a genuine concentration (kg/m^3), no
        # ρ_cp-style conversion involved -- unaffected by the heat-path fix above.
        vapor_density[i] = ρv_top - near_field_top_latent + far_field_latent + near_field_latent  # raw; Aitken-blended below
    end

    # Aitken needs the whole residual vector to fit omega, so this is a
    # separate pass over the raw solve above, not folded into that loop.
    _aitken_relax!(air_temperature, air_temperature_prev, layer_heights, canopy_height,
        aitken_omega, aitken_have_prev, aitken_residual_prev,
        aitken_omega_min, aitken_omega_max, aitken_weight_bottom, aitken_weight_top, aitken_bottom_emphasis, u"K")
    _aitken_relax!(vapor_density, vapor_density_prev, layer_heights, canopy_height,
        aitken_omega_latent, aitken_have_prev_latent, aitken_residual_prev_latent,
        aitken_omega_min, aitken_omega_max, aitken_weight_bottom, aitken_weight_top, aitken_bottom_emphasis, u"kg/m^3")

    # Absolute backstop anchored to T_top/ground_temperature -- leaf
    # temperature is clamped only relative to this array, so without an
    # absolute bound here the two can climb together, unanchored.
    air_bound_lo, air_bound_hi = _raupach_temperature_bounds(model.far_field_mode, T_top, ground_temperature,
        leaf_temperature, max_air_temperature_deviation, model.bulk_temperature_margin)
    @inbounds for i in 1:n
        air_temperature[i] = clamp(air_temperature[i], air_bound_lo, air_bound_hi)
    end

    all(isfinite, air_temperature) ||
        throw(DomainError(air_temperature, "RaupachLTheoryAirProfile produced a non-finite air temperature"))
    all(isfinite, vapor_density) ||
        throw(DomainError(vapor_density, "RaupachLTheoryAirProfile produced a non-finite vapor density"))

    @inbounds for i in 1:n
        ρv_sat = wet_air_properties(air_temperature[i], 1.0, atmospheric_pressure; vapour_pressure_equation).vapour_density
        relative_humidity[i] = clamp(vapor_density[i] / ρv_sat, 0.0, 1.0)  # same-unit quantities cancel to a bare Float64
    end

    return (; air_temperature, relative_humidity)
end
