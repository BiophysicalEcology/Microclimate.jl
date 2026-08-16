"""
    RaupachLTheoryAirProfile(; ground_velocity_std_factor=0.25, canopy_top_velocity_std_factor=1.25,
                              min_ground_resistance=10.0u"s/m", relaxation=0.7,
                              aitken_omega_min=0.02, aitken_omega_max=0.8,
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
    min_ground_resistance::GR = 10.0u"s/m"
    relaxation::RL = 0.7
    aitken_omega_min::OMIN = 0.02
    aitken_omega_max::OMAX = 0.8
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
    geometry = canopy_layer_heights(heights, canopy_height, n_layers)
    z = geometry.layer_heights
    Δz = geometry.layer_thickness

    ω₀ = clamp(
        Float64(model.relaxation),
        model.aitken_omega_min,
        model.aitken_omega_max,
    )

    σw = zeros(typeof(0.0u"m/s"), n_layers)
    λ⁻¹ = zeros(typeof(1.0 / (1.0u"m/s" * 1.0u"s")), n_layers)
    K = zeros(typeof(0.0u"m^2/s"), n_layers)
    r = zeros(typeof(0.0u"s/m"), n_layers)

    Cff_H = zeros(typeof(0.0u"J/m^3"), n_layers)
    Cff_v = zeros(typeof(0.0u"kg/m^3"), n_layers)

    R_top = zeros(typeof(0.0u"s/m"), n_layers)
    R_ground = zeros(typeof(0.0u"s/m"), n_layers)

    H_below = zeros(typeof(0.0u"W/m^2"), n_layers)
    Hn_weight = zeros(typeof(0.0u"W/m^2" / 1.0u"m/s"), n_layers)

    E_below = zeros(typeof(0.0u"kg/m^2/s"), n_layers)
    En_weight = zeros(typeof(1.0u"kg/m^2/s" / 1.0u"m/s"), n_layers)

    T = zeros(typeof(0.0u"K"), n_layers)
    T_prev = similar(T)

    ρv = zeros(typeof(0.0u"kg/m^3"), n_layers)
    ρv_prev = similar(ρv)

    RH = zeros(n_layers)

    return (;
        layer_heights = z,
        layer_thickness = Δz,
        vertical_velocity_std = σw,
        inv_near_field_length = λ⁻¹,
        eddy_diffusivity = K,
        layer_resistance = r,
        far_field_accum = Cff_H,
        far_field_accum_latent = Cff_v,

        resistance_from_top = R_top,
        resistance_to_ground = R_ground,
        cumulative_sensible_below = H_below,
        sensible_near_field_weight = Hn_weight,
        cumulative_latent_below = E_below,
        latent_near_field_weight = En_weight,

        air_temperature = T,
        air_temperature_prev = T_prev,
        vapor_density = ρv,
        vapor_density_prev = ρv_prev,
        relative_humidity = RH,

        aitken_omega = Ref(ω₀),
        aitken_omega_latent = Ref(ω₀),
        aitken_have_prev = Ref(false),
        aitken_have_prev_latent = Ref(false),
        aitken_residual_prev = zeros(typeof(0.0u"K"), n_layers),
        aitken_residual_prev_latent = zeros(typeof(0.0u"kg/m^3"), n_layers),

        ground_heat_conductance = Ref(0.0u"W/m^2/K"),
        ground_vapor_conductance = Ref(0.0u"m/s"),
    )
end

# Raupach (1989), eq. 34:
#   kₙ(ζ) ≈ c₁ log(1 - exp(-ζ)) + c₂ exp(-ζ),  ζ > 0
# The constants satisfy ∫₀^∞ kₙ(ζ)dζ = 1/2.
const _RAUPACH_KERNEL_C₁ = -1 / sqrt(2π)
const _RAUPACH_KERNEL_C₂ = 0.5 - π^2 / (6sqrt(2π))

@inline function _raupach_kernel(ζ)
    e⁻ζ = exp(-ζ)
    return _RAUPACH_KERNEL_C₁ * log1p(-e⁻ζ) + _RAUPACH_KERNEL_C₂ * e⁻ζ
end

# Odd extension: kₙ(-ζ) = -kₙ(ζ); eq. 34 gives it for ζ>0 only).
@inline _raupach_kernel_signed(ζ) =
    ζ >= zero(ζ) ? _raupach_kernel(ζ) : -_raupach_kernel(-ζ)

# Eq. 37: Cn = ∫ [S(z0)/σw(z0)] · [kn{(z-z0)/λ} + kn{(z+z0)/λ}] dz0 
# Raupach (1989), eq. 37:
#   Cₙ(z) = ∫ S(z₀)/σw(z₀) [kₙ((z-z₀)/λ) + kₙ((z+z₀)/λ)] dz₀
# — direct-path and reflected (ground image-source) terms each get their own 
# kernel evaluation. `source_idx` indexes into
# `layer_heights`/`layer_thickness`/`inv_near_field_length`.
@inline function _raupach_kernel_weight(
    evaluation_height,
    source_index,
    layer_heights,
    inv_near_field_length,
)
    z = evaluation_height
    j = source_index
    z₀ = layer_heights[j]
    λ⁻¹ = inv_near_field_length[j]

    ζd = (z - z₀) * λ⁻¹
    ζr = (z + z₀) * λ⁻¹

    return _raupach_kernel_signed(ζd) + _raupach_kernel_signed(ζr)
end

# Self-term: the kernel is singular at zero distance, so a source layer's
# contribution to a point within its own thickness can't use a single
# point-source evaluation. Subdivides the layer into `subdivisions` sub-points
# and averages — a plain midpoint quadrature of the singularity. Only the
# direct-path term is ever near-singular; the reflected term isn't.
@inline function _raupach_self_kernel_weight(evaluation_height, source_index,
    layer_heights, layer_thickness, inv_near_field_length, subdivisions,)

    z = evaluation_height
    j = source_index
    z₀ = layer_heights[j]
    Δz = layer_thickness[j]
    λ⁻¹ = inv_near_field_length[j]
    N = subdivisions

    total = 0.0

    @inbounds for k in 1:N
        zₖ = z₀ - Δz / 2 + Δz * (k - 0.5) / N
        ζd_signed = (z - zₖ) * λ⁻¹
        s = ζd_signed >= zero(ζd_signed) ? 1.0 : -1.0
        ζd = max(abs(ζd_signed), 1.0e-9)
        ζr = (z + zₖ) * λ⁻¹

        total += s * _raupach_kernel(ζd) + _raupach_kernel_signed(ζr)
    end

    return total / N
end

# :bulk's near-field kernel evaluates kn once (at the direct distance only)
# and multiplies by the raw sum of direct+reflected normalized distances,
# not eq. 37's separately-signed evaluation -- the same shortcut as its
# far-field. It also drops any self term outright rather than resolving it
# by subdivision quadrature.
@inline function _raupach_bulk_kernel_weight(evaluation_height,
    source_index, layer_heights, inv_near_field_length,)
    
    z = evaluation_height
    j = source_index
    z₀ = layer_heights[j]
    λ⁻¹ = inv_near_field_length[j]

    ζd = (z - z₀) * λ⁻¹
    ζr = (z + z₀) * λ⁻¹

    return _raupach_kernel(abs(ζd)) * (ζd + ζr)
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
    x = newv
    x₀ = oldv
    z = layer_heights
    h = canopy_height

    ωref = omega_ref
    have_prev = have_prev_ref
    r_prev = residual_prev

    ωmin = omega_min
    ωmax = omega_max
    w₀ = weight_bottom
    w₁ = weight_top
    β = bottom_emphasis

    n = length(x)
    Δw = w₁ - w₀

    if !have_prev[]
        ω = clamp(ωref[], ωmin, ωmax)

        @inbounds for i in 1:n
            r = x[i] - x₀[i]
            wz = w₀ + Δw * (z[i] / h)^2
            r_prev[i] = r
            x[i] = x₀[i] + ω * wz * r
        end

        ωref[] = ω
        have_prev[] = true
        return nothing
    end

    numerator = 0.0
    denominator = 0.0

    @inbounds for i in 1:n
        r = x[i] - x₀[i]
        Δr = r - r_prev[i]
        ξ = 1.0 - z[i] / h
        w = 1.0 + β * ξ^2

        numerator += w * ustrip(unit, r_prev[i]) * ustrip(unit, Δr)
        denominator += w * ustrip(unit, Δr)^2
    end

    ω = denominator > 0.0 ? -ωref[] * numerator / denominator : ωref[]
    ω = clamp(ω, ωmin, ωmax)

    @inbounds for i in 1:n
        r = x[i] - x₀[i]
        wz = w₀ + Δw * (z[i] / h)^2
        r_prev[i] = r
        x[i] = x₀[i] + ω * wz * r
    end

    ωref[] = ω
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
    Cff_H = far_field_accum
    Cff_v = far_field_accum_latent
    H_below = cumulative_sensible_below
    E_below = cumulative_latent_below
    r = layer_resistance

    Tg = ground_temperature
    ρv_g = ground_vapor_density
    T₀ = air_temperature_prev
    ρv₀ = vapor_density_prev
    Rg = ground_resistance
    ρcp = ρ_cp

    Fg = ρcp * (Tg - T₀[n]) / Rg
    Eg = (ρv_g - ρv₀[n]) / Rg

    Cff_H[1] = (H_below[1] + Fg) * r[1]
    Cff_v[1] = (E_below[1] + Eg) * r[1]

    @inbounds for i in 2:n
        Cff_H[i] = Cff_H[i - 1] + (H_below[i] + Fg) * r[i]
        Cff_v[i] = Cff_v[i - 1] + (E_below[i] + Eg) * r[i]
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

    Cff_H = far_field_accum
    Cff_v = far_field_accum_latent
    Rtop = resistance_from_top
    Rgnd = resistance_to_ground
    H_below = cumulative_sensible_below
    E_below = cumulative_latent_below
    r = layer_resistance

    Tg = ground_temperature
    ρv_g = ground_vapor_density
    T₀ = air_temperature_prev
    ρv₀ = vapor_density_prev
    Rg = ground_resistance

    Rtop[1] = r[1]

    @inbounds for i in 2:n
        Rtop[i] = Rtop[i - 1] + r[i]
    end

    Rgnd[n] = Rg

    @inbounds for i in (n - 1):-1:1
        Rgnd[i] = Rgnd[i + 1] + r[i]
    end

    @inbounds for i in 1:n
        ρcp = calc_ρ_cp(T₀[i])
        Fg = ρcp * (Tg - T₀[i]) / Rgnd[i]
        Eg = (ρv_g - ρv₀[i]) / Rgnd[i]

        Cff_H[i] = (H_below[i] + Fg) * Rtop[i]
        Cff_v[i] = (E_below[i] + Eg) * Rtop[i]
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

    # Mathematical aliases used within the implementation.
    z = layer_heights
    Δz = layer_thickness
    σw_profile = vertical_velocity_std
    λ⁻¹ = inv_near_field_length
    K = eddy_diffusivity
    r = layer_resistance

    Cff_H = far_field_accum
    Cff_v = far_field_accum_latent
    Rtop = resistance_from_top
    Rground = resistance_to_ground

    H_below = cumulative_sensible_below
    Hn_weight = sensible_near_field_weight
    E_below = cumulative_latent_below
    En_weight = latent_near_field_weight

    T = air_temperature
    T₀ = air_temperature_prev
    ρv = vapor_density
    ρv₀ = vapor_density_prev
    RH = relative_humidity

    h = canopy_height
    d = displacement_height
    u★ = friction_velocity
    H = sensible_heat_source
    E = evaporation_mass_flow
    p = atmospheric_pressure

    Ttop = canopy_top_air_temperature
    RHtop = canopy_top_relative_humidity
    Tg = ground_temperature
    RHg = ground_relative_humidity

    n = length(T)

    length(H) == n || throw(
        ArgumentError("sensible_heat_source must have one entry per canopy layer"),
    )

    length(E) == n || throw(
        ArgumentError("evaporation_mass_flow must have one entry per canopy layer"),
    )

    ρv_top = wet_air_properties(
        Ttop,
        RHtop,
        p;
        vapour_pressure_equation,
    ).vapour_density

    # Preserve the previous Picard iterate.
    T₀ .= T

    if all(iszero, ρv)
        fill!(ρv, ρv_top)
    end

    ρv₀ .= ρv

    (;
        ground_velocity_std_factor,
        canopy_top_velocity_std_factor,
        min_ground_resistance,
        aitken_omega_min,
        aitken_omega_max,
        aitken_weight_bottom,
        aitken_weight_top,
        aitken_bottom_emphasis,
        near_field_subdivisions,
        max_lagrangian_timescale,
        max_air_temperature_deviation,
    ) = model

    # Raupach's notation.
    a₀ = ground_velocity_std_factor
    a₁ = canopy_top_velocity_std_factor
    γ = boundary_layer_model.dyer_constant
    κ = boundary_layer_model.karman_constant

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

    # Match K(h) = σw(h)²T_L to the above-canopy MOST diffusivity:
    #
    #   K(h) = κu★(h-d)/Φh
    #   a₂ = κ(1-d/h)/(a₁²Φh)
    #   T_L = min(a₂h/u★, T_L,max)
    z★ = max(h - d, 1.0e-3u"m")

    Φh = isfinite(obukhov_length) ?
        calc_Φ_h(
            z★,
            γ,
            obukhov_length,
            boundary_layer_model.stable_Φ_h_coefficient,
            boundary_layer_model.min_stable_Φ_h,
            boundary_layer_model.max_stable_Φ_h,
        ) : 1.0

    a₂ = κ * (1.0 - d / h) / (a₁^2 * Φh)
    TL = min(a₂ * h / u★, max_lagrangian_timescale)

    σw_mean = (a₁ + a₀) * u★ / 2
    σw_amplitude = (a₁ - a₀) * u★ / 2

    @inbounds for i in 1:n
        σw = σw_mean + σw_amplitude * cos(π * (1.0 - z[i] / h))
        σw_profile[i] = σw
        λ⁻¹[i] = 1 / (σw * TL)
        K[i] = σw^2 * TL
        r[i] = Δz[i] / K[i]
        Hn_weight[i] = H[i] / σw
        En_weight[i] = (E[i] / 1.0u"m^2") / σw
    end

    # Ground exchange uses the wind-attenuation K profile, consistently with
    # the K-theory ground boundary.
    # Shared by ground_heat_conductance/ground_vapor_conductance and _raupach_far_field!'s
    # ground flux -- uses wind_attenuation_profile's shape, not layer_resistance[n]'s σ_w profile.
    Ktop = κ * u★ * z★
    Kg = Ktop * wind_attenuation[n]
    Rg = max(Δz[n] / Kg, min_ground_resistance)

    # ρcp is volumetric heat capacity, with units J m⁻³ K⁻¹.
    # ρ_cp fixed from T_top/ground_temperature, not from the (possibly still
    # drifting) air_temperature_prev iterate: calc_ρ_cp ∝ 1/T, so tying it to
    # each layer's own previous pass makes warming self-amplifying (ΔT ∝
    # flux·T_prev) -- same failure mode as canopy_top_flux_boundary's ρ_cp,
    # here spread across every layer/pass instead of one fixed-point loop.
    ρcp = calc_ρ_cp((Ttop + Tg) / 2)

    # gH = ρcp/Rg [W m⁻² K⁻¹], gv = 1/Rg [m s⁻¹].
    ground_heat_conductance[] = ρcp / Rg
    ground_vapor_conductance[] = 1 / Rg

    H_below[n] = H[n]
    E_below[n] = E[n] / 1.0u"m^2"

    @inbounds for i in (n - 1):-1:1
        H_below[i] = H_below[i + 1] + H[i]
        E_below[i] = E_below[i + 1] + E[i] / 1.0u"m^2"
    end

    ρv_g = wet_air_properties(Tg, RHg, p; vapour_pressure_equation,).vapour_density

    _raupach_far_field!(
        model.far_field_mode,
        Cff_H,
        Cff_v,
        Rtop,
        Rground,
        H_below,
        E_below,
        r,
        Δz,
        Tg,
        ρv_g,
        T₀,
        ρv₀,
        Rg,
        ρcp,
        n,
    )

    μ = _raupach_near_field_scale(model.far_field_mode, n)

    # Near-field concentration at the canopy top, subtracted from every
    # layer's own near-field term below. Layers whose height is (numerically)
    # coincident with canopy_height use the subdivided self kernel; every
    # other layer uses a single point-source evaluation. Heat and vapor share
    # one kernel evaluation per layer (same turbulence field).

    Cn_top_H = zero(eltype(Cff_H))
    Cn_top_v = zero(ρv_top)
    
    @inbounds for j in 1:n
        w = _raupach_top_kernel_weight(
            model.far_field_mode,
            j,
            h,
            z,
            Δz,
            λ⁻¹,
            near_field_subdivisions,
        )

        Cn_top_H += Hn_weight[j] * w
        Cn_top_v += En_weight[j] * w
    end

    Cn_top_H *= μ
    Cn_top_v *= μ

    @inbounds for i in 1:n
        Cn_H = zero(Cn_top_H)
        Cn_v = zero(ρv_top)

        for j in 1:n
            w = _raupach_pair_kernel_weight(
                model.far_field_mode,
                i,
                j,
                z,
                Δz,
                λ⁻¹,
                near_field_subdivisions,
            )

            Cn_H += Hn_weight[j] * w
            Cn_v += En_weight[j] * w
        end

        Cn_H *= μ
        Cn_v *= μ

        # Heat: T₀ is an absolute temperature baseline, not a concentration;
        # ρcp converts the source-driven perturbation (a volumetric
        # energy density) into a temperature increment added on top of it.
        ρcpᵢ = _raupach_layer_ρ_cp(model.far_field_mode, ρcp, T₀, i)

        # ΔT = energy-density perturbation / volumetric heat capacity.
        T[i] = Ttop + (Cff_H[i] + Cn_H - Cn_top_H) / ρcpᵢ

        # Vapor terms already have concentration units, kg m⁻³.
        ρv[i] = ρv_top + Cff_v[i] + Cn_v - Cn_top_v
    end

    # Aitken needs the whole residual vector to fit omega, so this is a
    # separate pass over the raw solve above, not folded into that loop.
    _aitken_relax!(
        T,
        T₀,
        z,
        h,
        aitken_omega,
        aitken_have_prev,
        aitken_residual_prev,
        aitken_omega_min,
        aitken_omega_max,
        aitken_weight_bottom,
        aitken_weight_top,
        aitken_bottom_emphasis,
        u"K",
    )

    _aitken_relax!(
        ρv,
        ρv₀,
        z,
        h,
        aitken_omega_latent,
        aitken_have_prev_latent,
        aitken_residual_prev_latent,
        aitken_omega_min,
        aitken_omega_max,
        aitken_weight_bottom,
        aitken_weight_top,
        aitken_bottom_emphasis,
        u"kg/m^3",
    )

    # Absolute backstop anchored to T_top/ground_temperature -- leaf
    # temperature is clamped only relative to this array, so without an
    # absolute bound here the two can climb together, unanchored.
    Tlo, Thi = _raupach_temperature_bounds(
        model.far_field_mode,
        Ttop,
        Tg,
        leaf_temperature,
        max_air_temperature_deviation,
        model.bulk_temperature_margin,
    )

    @inbounds for i in 1:n
        T[i] = clamp(T[i], Tlo, Thi)
    end

    all(isfinite, T) || throw(
        DomainError(T, "RaupachLTheoryAirProfile produced a non-finite air temperature"),
    )

    all(isfinite, ρv) || throw(
        DomainError(ρv, "RaupachLTheoryAirProfile produced a non-finite vapor density"),
    )

    @inbounds for i in 1:n
        ρv_sat = wet_air_properties(
            T[i],
            1.0,
            p;
            vapour_pressure_equation,
        ).vapour_density

        RH[i] = clamp(ρv[i] / ρv_sat, 0.0, 1.0)
    end

    return (;
        air_temperature = T,
        relative_humidity = RH,
    )

end
