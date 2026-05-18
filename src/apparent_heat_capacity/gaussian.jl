"""
    Gaussian(; T_melt=0.0u"°C", width=0.2u"K")

Gaussian-bump apparent heat capacity. The augmentation is a normal distribution
centred at the melting point `T_melt`, with width parameter `σ = width`:

    cp_eff(T) = cp_base + (L_fusion / (σ · √(2π))) · exp(-½ · ((T - T_melt) / σ)²)

The total latent budget `∫(cp_eff - cp_base) dT` integrates to **exactly** `L_fusion`,
since the Gaussian is a normalised PDF. Width `σ` controls the spread of the
latent absorption around `T_melt`.

References: Comini, Del Guidice, Lewis & Zienkiewicz (1974), *Finite element solution
of non-linear heat conduction problems with special reference to phase change*,
Int. J. Numer. Methods Engng 8(3), 613–624 — introduced the use of smooth bumps
for apparent-cp in finite-element phase-change problems. The Gaussian shape in
particular has been used by Civan & Sliepcevich (1984, 1987) and subsequent
computational-heat-transfer literature.

Unlike [`BonacinaStep`](@ref) and [`TanhSmoothed`](@ref), the Gaussian has *infinite
support* — cp is augmented at all temperatures, though by an exponentially small
amount far from T_melt. For typical `width=0.2 K` the augmentation drops below
0.1% of peak at |T - T_melt| > 0.7 K.
"""
struct Gaussian{TM,W} <: AbstractApparentHeatCapacity
    T_melt::TM
    width::W
end
Gaussian(; T_melt=0.0u"°C", width=0.2u"K") = Gaussian(T_melt, width)

@inline function apparent_heat_capacity(m::Gaussian, T, cp_base, L_fusion)
    T_c = u"°C"(T)
    σ_K = ustrip(u"K", m.width)
    x = ustrip(u"K", T_c - m.T_melt) / σ_K
    # peak amplitude = L_fusion / (σ · √(2π)) so the integral over T equals L_fusion
    peak = L_fusion / (m.width * sqrt(2π))
    return cp_base + peak * exp(-0.5 * x * x)
end
