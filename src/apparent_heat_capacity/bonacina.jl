"""
    BonacinaStep(; T_low=-0.45u"°C", T_high=0.4u"°C", window=1.0u"K")

Step-function apparent-heat-capacity for ice/water phase change, following the
formulation of Bonacina, Comini, Fasano & Primicerio (1973), *Numerical solution
of phase-change problems*, Int. J. Heat Mass Transfer 16(10), 1825–1832. The
specific heat is augmented by `L_fusion / window` over a fixed band:

    cp_eff(T) = cp_base + L_fusion / window     if T_low < T ≤ T_high
              = cp_base                          otherwise

The defaults (-0.45°C, 0.4°C, 1 K) reproduce NicheMapR Fortran's `SOILPROPS.f`
exactly. The integrated latent budget across the band is
`L_fusion · (T_high - T_low) / window` = `0.85 · L_fusion` with the defaults —
slightly under-counting fusion, but matched to the legacy Fortran for parity.

The discontinuous shape stresses adaptive RK solvers: Tsit5 tends to pin T at the
band edge instead of crossing it, while Adams predictor-correctors (Fortran SFODE)
overshoot uncontrolled. See [`TanhSmoothed`](@ref) for a smoothed variant that any
solver can integrate cleanly.
"""
struct BonacinaStep{TL,TH,W} <: AbstractApparentHeatCapacity
    T_low::TL
    T_high::TH
    window::W
end
BonacinaStep(; T_low=-0.45u"°C", T_high=0.4u"°C", window=1.0u"K") =
    BonacinaStep(T_low, T_high, window)

@inline function apparent_heat_capacity(m::BonacinaStep, T, cp_base, L_fusion)
    T_c = u"°C"(T)
    if T_c > m.T_low && T_c <= m.T_high
        return cp_base + L_fusion / m.window
    end
    return cp_base
end
