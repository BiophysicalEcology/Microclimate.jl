"""
    TanhSmoothed(; T_low=-0.45u"°C", T_high=0.4u"°C", window=1.0u"K", smoothing=0.1u"K")

C∞-smooth replacement for [`BonacinaStep`](@ref). The Heaviside indicator of the
step formulation is replaced by a pair of `tanh` ramps centred at the band edges:

    indicator(T) = 0.5 · [tanh((T - T_low)/δ) - tanh((T - T_high)/δ)]
    cp_eff(T)   = cp_base + (L_fusion / window) · indicator(T)

where `δ = smoothing`. `indicator(T)` is ≈ 1 inside the band, ≈ 0 outside, with
smooth ramps of width ~δ at each edge. For `δ ≪ T_high - T_low` the integrated
latent budget matches [`BonacinaStep`](@ref) to within a relative error of
`δ / (T_high - T_low)`.

**No single literature reference.** `tanh`-smoothing of step-function apparent-heat-
capacity is a generic numerical-analysis technique in computational heat transfer;
the double-edge form here was chosen because:

- `cp_eff` is C∞ in `T`, so adaptive RK solvers see no discontinuity in cp or any
  of its derivatives. Step rejection and pinning at the band edges disappear.
- For `δ ≪ window` the plateau interior is preserved, so the integrated latent
  budget equals the step formulation's. Physics calibration is unchanged.
- Each ramp is monotonic and bounded, so no spurious oscillations in cp.

Physically the smooth ramp is also more realistic than a strict step: real snow
melts over a small range because of mixed grain sizes, dissolved impurities
depressing the melting point, and pre-melted quasi-liquid layers on ice grains.
The default `smoothing = 0.1 K` was picked to keep the plateau wide while letting
solvers cross the edges in normal-sized steps.
"""
struct TanhSmoothed{TL,TH,W,D} <: AbstractApparentHeatCapacity
    T_low::TL
    T_high::TH
    window::W
    smoothing::D
end
TanhSmoothed(; T_low=-0.45u"°C", T_high=0.4u"°C", window=1.0u"K", smoothing=0.1u"K") =
    TanhSmoothed(T_low, T_high, window, smoothing)

@inline function apparent_heat_capacity(m::TanhSmoothed, T, cp_base, L_fusion)
    T_c = u"°C"(T)
    # Strip to dimensionless before tanh — Unitful's tanh fallback is fragile for
    # mixed °C/K Quantities.
    δ_K = ustrip(u"K", m.smoothing)
    arg_lo = ustrip(u"K", T_c - m.T_low)  / δ_K
    arg_hi = ustrip(u"K", T_c - m.T_high) / δ_K
    indicator = 0.5 * (tanh(arg_lo) - tanh(arg_hi))
    return cp_base + (L_fusion / m.window) * indicator
end
