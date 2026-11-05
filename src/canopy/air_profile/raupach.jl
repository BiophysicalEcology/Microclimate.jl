"""
    RaupachLTheoryAirProfile(; ground_velocity_std_factor=0.25, canopy_top_velocity_std_factor=1.25,
                              far_field_mode=Val(:exact), max_lagrangian_timescale=200.0u"s")

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
(`a0·u*`), maximum at canopy top (`a1·u*`). The ground resistance is the ground-most
layer's turbulent resistance plus a quasi-laminar sublayer term (ground-local u★,
Garratt et al. 1988's sublayer_stanton).

The near-field kernel is singular at zero distance, so its self term is the
kernel's exact average over the source layer's own thickness, from its
closed-form integral.

`T_L` (see `canopy_heat_vapor_transport!`) is `∝ 1/friction_velocity`
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
separately-signed kernel evaluation with the self term from the kernel's
closed-form integral. `Val(:bulk)` is an alternate formulation used in 
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
@kwdef struct RaupachLTheoryAirProfile{GVSF,CVSF,FFM,MTL} <: AbstractCanopyAirProfileModel
    ground_velocity_std_factor::GVSF = 0.25
    canopy_top_velocity_std_factor::CVSF = 1.25
    far_field_mode::FFM = Val(:exact)
    max_lagrangian_timescale::MTL = 200.0u"s"
end

function allocate_air_profile(model::RaupachLTheoryAirProfile, canopy_height, displacement_height, plant_area_index, heights, n_layers, boundary_layer_model)
    geometry = canopy_layer_heights(heights, canopy_height, n_layers)
    z = geometry.layer_heights
    Δz = geometry.layer_thickness

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
        canopy_height,
        displacement_height,
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

        ground_heat_conductance = Ref(0.0u"W/m^2/K"),
        ground_vapor_conductance = Ref(0.0u"m/s"),
    )
end

# Raupach (1989), eq. 34, with ζ the source distance over σw·T_L:
#   kₙ(ζ) ≈ c₁ log(1 - exp(-|ζ|)) + c₂ exp(-|ζ|)
# The constants satisfy ∫₀^∞ kₙ(ζ)dζ = 1/2.
const _RAUPACH_KERNEL_C₁ = -1 / sqrt(2π)
const _RAUPACH_KERNEL_C₂ = 0.5 - π^2 / (6sqrt(2π))

function _raupach_kernel(ζ)
    e⁻ζ = exp(-abs(ζ))
    return _RAUPACH_KERNEL_C₁ * log1p(-e⁻ζ) + _RAUPACH_KERNEL_C₂ * e⁻ζ
end

# Raupach (1989), eq. 37:
#   Cₙ(z) = ∫ S(z₀)/σw(z₀) [kₙ((z-z₀)/λ) + kₙ((z+z₀)/λ)] dz₀
# — direct-path and reflected (ground image-source) terms. With
# eᵢ = exp(-zλ⁻¹) and e₀ = exp(-z₀λ⁻¹), the reflected term's exponential is
# eᵢe₀ and the direct term's is eᵢ/e₀ or e₀/eᵢ, so the pair takes one exp
# and one log.
function _raupach_kernel_weight(
    evaluation_height,
    source_index,
    layer_heights,
    inv_near_field_length,
)
    z = evaluation_height
    j = source_index
    z₀ = layer_heights[j]
    λ⁻¹ = inv_near_field_length[j]

    eᵢ = exp(-z * λ⁻¹)
    e₀ = exp(-z₀ * λ⁻¹)
    direct = z > z₀ ? eᵢ / e₀ : e₀ / eᵢ
    reflected = eᵢ * e₀
    return _RAUPACH_KERNEL_C₁ * log((1 - direct) * (1 - reflected)) + _RAUPACH_KERNEL_C₂ * (direct + reflected)
end

# ∫₀^ζ kₙ(t) dt for ζ ≥ 0, from ∫ log(1 - exp(-t)) dt = Li₂(exp(-t)).
function _raupach_kernel_integral(ζ)
    e⁻ζ = exp(-ζ)
    return _RAUPACH_KERNEL_C₁ * (reli2(e⁻ζ) - π^2 / 6) + _RAUPACH_KERNEL_C₂ * (1 - e⁻ζ)
end

# Self-term: the kernel is singular at zero distance, so a source layer's
# contribution to a point within its own thickness is the kernel's average
# over that thickness, from its closed-form integral. The kernel is even, so
# a point inside the layer sums the integrals to either edge, and a point
# outside it takes their difference.
function _raupach_self_kernel_weight(evaluation_height, source_index,
    layer_heights, layer_thickness, inv_near_field_length,)

    z = evaluation_height
    j = source_index
    z₀ = layer_heights[j]
    Δz = layer_thickness[j]
    λ⁻¹ = inv_near_field_length[j]
    bottom = z₀ - Δz / 2
    top = z₀ + Δz / 2

    to_bottom = _raupach_kernel_integral(abs((z - bottom) * λ⁻¹))
    to_top = _raupach_kernel_integral(abs((top - z) * λ⁻¹))
    direct = bottom <= z <= top ? to_bottom + to_top : abs(to_bottom - to_top)
    reflected = _raupach_kernel_integral((z + top) * λ⁻¹) - _raupach_kernel_integral((z + bottom) * λ⁻¹)
    return (direct + reflected) / (Δz * λ⁻¹)
end

# :bulk's near-field kernel evaluates kn once (at the direct distance only)
# and multiplies by the raw sum of direct+reflected normalized distances,
# not eq. 37's separately-signed evaluation -- the same shortcut as its
# far-field. It also drops any self term outright rather than resolving it
# by subdivision quadrature.
function _raupach_bulk_kernel_weight(evaluation_height,
    source_index, layer_heights, inv_near_field_length,)
    
    z = evaluation_height
    j = source_index
    z₀ = layer_heights[j]
    λ⁻¹ = inv_near_field_length[j]

    ζd = (z - z₀) * λ⁻¹
    ζr = (z + z₀) * λ⁻¹

    return _raupach_kernel(ζd) * (ζd + ζr)
end

# Top-boundary near-field weight for source layer i, dispatched on
# far_field_mode. Layer 1 (closest to canopy_height) is always the self term
# -- matching _raupach_pair_kernel_weight's own j==i criterion, not a
# height-coincidence check, since layer_heights[1] need not sit exactly at
# canopy_height and the two criteria disagreeing double-counts/drops the
# singular kernel inconsistently between near_field_top and near_field[1].
function _raupach_top_kernel_weight(::Val{:exact}, i, canopy_height, layer_heights, layer_thickness, inv_near_field_length)
    return i == 1 ?
        _raupach_self_kernel_weight(canopy_height, i, layer_heights, layer_thickness, inv_near_field_length) :
        _raupach_kernel_weight(canopy_height, i, layer_heights, inv_near_field_length)
end
function _raupach_top_kernel_weight(::Val{:bulk}, i, canopy_height, layer_heights, layer_thickness, inv_near_field_length)
    return i == 1 ? 0.0 : _raupach_bulk_kernel_weight(canopy_height, i, layer_heights, inv_near_field_length)
end

# Near-field weight for one (eval layer i, source layer j) pair.
function _raupach_pair_kernel_weight(::Val{:exact}, i, j, layer_heights, layer_thickness, inv_near_field_length)
    return j == i ?
        _raupach_self_kernel_weight(layer_heights[i], j, layer_heights, layer_thickness, inv_near_field_length) :
        _raupach_kernel_weight(layer_heights[i], j, layer_heights, inv_near_field_length)
end
function _raupach_pair_kernel_weight(::Val{:bulk}, i, j, layer_heights, layer_thickness, inv_near_field_length)
    return j == i ? 0.0 : _raupach_bulk_kernel_weight(layer_heights[i], j, layer_heights, inv_near_field_length)
end

# Eq. 19a: Cf(z)-Cf(zR) = ∫[z,zR] F(z')/Kf(z') dz' -- each layer's own local
# flux over its own local resistance, summed; not one flux value times the
# total resistance across the range. Fg is a single ground flux. The bottom layer's
# concentration is its offset (canopy top plus near field) plus Fg across the
# in-canopy resistance R = Σr, so Fg solves through the ground resistance and R in
# series, from the ground to that offset: the same solution as Fg from the bottom
# layer's own concentration, without feeding that back through Fg.
function _raupach_far_field!(::Val{:exact}, far_field_accum, far_field_accum_latent, _resistance_from_top, _resistance_to_ground,
    cumulative_sensible_below, cumulative_latent_below, layer_resistance, _layer_thickness,
    ground_temperature, ground_vapor_density, _air_temperature_prev, _vapor_density_prev,
    ground_resistance_heat, ground_resistance_vapor, ρ_cp, n; bottom_offset_heat, bottom_offset_vapor)

    Cff_H = far_field_accum
    Cff_v = far_field_accum_latent
    H_below = cumulative_sensible_below
    E_below = cumulative_latent_below
    r = layer_resistance

    Tg = ground_temperature
    ρv_g = ground_vapor_density
    Rg_H = ground_resistance_heat
    Rg_v = ground_resistance_vapor
    ρcp = ρ_cp

    R = sum(r)
    offset_H = bottom_offset_heat
    offset_v = bottom_offset_vapor
    for i in 1:n
        offset_H += H_below[i] * r[i]
        offset_v += E_below[i] * r[i]
    end
    Fg = (ρcp * Tg - offset_H) / (Rg_H + R)
    Eg = (ρv_g - offset_v) / (Rg_v + R)

    Cff_H[1] = (H_below[1] + Fg) * r[1]
    Cff_v[1] = (E_below[1] + Eg) * r[1]

    for i in 2:n
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
    ground_temperature, ground_vapor_density, air_temperature_prev, vapor_density_prev,
    _ground_resistance_heat, ground_resistance_vapor, _ρ_cp, n; _...)

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
    # Rgnd is shared heat/vapor here (unlike :exact); vapor resistance for both.
    Rg = ground_resistance_vapor

    Rtop[1] = r[1]

    for i in 2:n
        Rtop[i] = Rtop[i - 1] + r[i]
    end

    Rgnd[n] = Rg

    for i in (n - 1):-1:1
        Rgnd[i] = Rgnd[i + 1] + r[i]
    end

    for i in 1:n
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
_raupach_near_field_scale(::Val{:exact}, n) = 1.0
_raupach_near_field_scale(::Val{:bulk}, n) = 1.0 + 0.894 * exp(-0.01386 * n) + 9.82 * exp(-0.15 * n)

# :exact uses one ρ_cp fixed from ambient data; :bulk recomputes per layer
# from air_temperature_prev[i], matching micropoint's LangrangianOne (see
# _raupach_far_field!'s own comment on why :bulk keeps this).
_raupach_layer_ρ_cp(::Val{:exact}, ρ_cp, _air_temperature_prev, _i) = ρ_cp
_raupach_layer_ρ_cp(::Val{:bulk}, _ρ_cp, air_temperature_prev, i) = calc_ρ_cp(air_temperature_prev[i])

function canopy_heat_vapor_transport!(buffers, model::RaupachLTheoryAirProfile, boundary_layer_model,
    forcing::CanopyAirProfileForcing; site, vapour_pressure_equation=GoffGratch(),
)
    (; friction_velocity, inverse_obukhov_length, wind_speed, canopy_top_air_temperature, canopy_top_vapor_density,
       ground_temperature, ground_relative_humidity, sensible_heat_source, evaporation_mass_flow,
       atmospheric_pressure) = forcing
    ground_roughness_height = site.roughness_height
    (; canopy_height, displacement_height, layer_heights, layer_thickness, vertical_velocity_std, inv_near_field_length, eddy_diffusivity,
       layer_resistance, far_field_accum, far_field_accum_latent, resistance_from_top, resistance_to_ground,
       cumulative_sensible_below, sensible_near_field_weight, cumulative_latent_below, latent_near_field_weight,
       air_temperature, air_temperature_prev, vapor_density, vapor_density_prev,
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

    h = canopy_height
    d = displacement_height
    u★ = friction_velocity
    H = sensible_heat_source
    E = evaporation_mass_flow
    p = atmospheric_pressure

    Ttop = canopy_top_air_temperature
    ρv_top = canopy_top_vapor_density
    Tg = ground_temperature
    RHg = ground_relative_humidity

    n = length(T)

    length(H) == n || throw(
        ArgumentError("sensible_heat_source must have one entry per canopy layer"),
    )

    length(E) == n || throw(
        ArgumentError("evaporation_mass_flow must have one entry per canopy layer"),
    )

    (;
        ground_velocity_std_factor,
        canopy_top_velocity_std_factor,
        max_lagrangian_timescale,
    ) = model

    # Raupach's notation.
    a₀ = ground_velocity_std_factor
    a₁ = canopy_top_velocity_std_factor
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
    Φh = calc_Φ_h(boundary_layer_model.stability_functions,
        stability_parameter(InverseObukhovLength(), z★, inverse_obukhov_length))

    a₂ = κ * (1.0 - d / h) / (a₁^2 * Φh)
    TL = min(a₂ * h / u★, max_lagrangian_timescale)

    σw_mean = (a₁ + a₀) * u★ / 2
    σw_amplitude = (a₁ - a₀) * u★ / 2

    for i in 1:n
        σw = σw_mean + σw_amplitude * cos(π * (1.0 - z[i] / h))
        σw_profile[i] = σw
        λ⁻¹[i] = 1 / (σw * TL)
        K[i] = σw^2 * TL
        r[i] = Δz[i] / K[i]
        Hn_weight[i] = H[i] / σw
        En_weight[i] = (E[i] / 1.0u"m^2") / σw
    end

    # Ground resistance = turbulent (r[n]) + quasi-laminar sublayer, in series. The
    # sublayer uses a ground-local u★ (log law on wind_speed[n] against
    # ground_roughness_height; allocate_canopy checks z[n] is above that roughness).
    u★_ground = κ * wind_speed[n] / log(z[n] / ground_roughness_height)
    kinematic_viscosity = dry_air_properties((Ttop + Tg) / 2, p).kinematic_viscosity
    sublayer_heat_resistance = 1.0 / (u★_ground * sublayer_stanton(ground_roughness_height, u★_ground, kinematic_viscosity))
    Rg_heat = r[n] + sublayer_heat_resistance
    Rg_vapor = r[n] + sublayer_heat_resistance / LEWIS_HEAT_TO_MASS_RATIO

    # ρcp is volumetric heat capacity, with units J m⁻³ K⁻¹.
    # ρ_cp fixed from T_top/ground_temperature, not from the (possibly still
    # drifting) air_temperature_prev iterate: calc_ρ_cp ∝ 1/T, so tying it to
    # each layer's own previous pass makes warming self-amplifying (ΔT ∝
    # flux·T_prev) -- same failure mode as canopy_top_flux_boundary's ρ_cp,
    # here spread across every layer/pass instead of one fixed-point loop.
    ρcp = calc_ρ_cp((Ttop + Tg) / 2)

    ground_heat_conductance[] = ρcp / Rg_heat
    ground_vapor_conductance[] = 1 / Rg_vapor

    H_below[n] = H[n]
    E_below[n] = E[n] / 1.0u"m^2"

    for i in (n - 1):-1:1
        H_below[i] = H_below[i + 1] + H[i]
        E_below[i] = E_below[i + 1] + E[i] / 1.0u"m^2"
    end

    ρv_g = wet_air_properties(Tg, RHg, p; vapour_pressure_equation,).vapour_density

    μ = _raupach_near_field_scale(model.far_field_mode, n)

    # Near-field concentration at the canopy top, subtracted from every
    # layer's own near-field term below. Heat and vapor share one kernel
    # evaluation per layer (same turbulence field).
    Cn_top_H = zero(eltype(Cff_H))
    Cn_top_v = zero(ρv_top)
    for j in 1:n
        w = _raupach_top_kernel_weight(model.far_field_mode, j, h, z, Δz, λ⁻¹)
        Cn_top_H += Hn_weight[j] * w
        Cn_top_v += En_weight[j] * w
    end
    Cn_top_H *= μ
    Cn_top_v *= μ

    bottom = _raupach_near_field(model.far_field_mode, n, z, Δz, λ⁻¹, Hn_weight, En_weight, μ)
    _raupach_far_field!(model.far_field_mode, Cff_H, Cff_v, Rtop, Rground, H_below, E_below, r, Δz,
        Tg, ρv_g, T₀, ρv₀, Rg_heat, Rg_vapor, ρcp, n;
        bottom_offset_heat=Ttop * ρcp + bottom.heat - Cn_top_H, bottom_offset_vapor=ρv_top + bottom.vapor - Cn_top_v)

    for i in 1:n
        Cn = _raupach_near_field(model.far_field_mode, i, z, Δz, λ⁻¹, Hn_weight, En_weight, μ)

        # Heat: T₀ is an absolute temperature baseline, not a concentration;
        # ρcp converts the source-driven perturbation (a volumetric
        # energy density) into a temperature increment added on top of it.
        ρcpᵢ = _raupach_layer_ρ_cp(model.far_field_mode, ρcp, T₀, i)

        # ΔT = energy-density perturbation / volumetric heat capacity.
        T[i] = Ttop + (Cff_H[i] + Cn.heat - Cn_top_H) / ρcpᵢ

        # Vapor terms already have concentration units, kg m⁻³.
        ρv[i] = ρv_top + Cff_v[i] + Cn.vapor - Cn_top_v
    end

    return nothing
end

# Near-field heat and vapor concentration at layer i from every layer's source
function _raupach_near_field(far_field_mode, i, z, Δz, λ⁻¹, Hn_weight, En_weight, μ)
    heat = zero(eltype(Hn_weight))
    vapor = zero(eltype(En_weight))
    for j in eachindex(Hn_weight)
        w = _raupach_pair_kernel_weight(far_field_mode, i, j, z, Δz, λ⁻¹)
        heat += Hn_weight[j] * w
        vapor += En_weight[j] * w
    end
    return (; heat=μ * heat, vapor=μ * vapor)
end
