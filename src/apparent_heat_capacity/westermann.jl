"""
    WestermannSigmoid(; T_melt=0.0u"°C", width=0.2u"K")

Phase-fraction-based apparent heat capacity. Model unfrozen-water mass fraction
as a logistic / `tanh` sigmoid centred at `T_melt` with characteristic width `W`:

    f_unfrozen(T) = ½ · (1 + tanh((T - T_melt) / W))

The apparent heat capacity is the derivative-of-fraction times latent heat:

    cp_eff(T) = cp_base + L_fusion · df_unfrozen/dT
              = cp_base + (L_fusion / (2W)) · sech²((T - T_melt) / W)

`sech²` is a bell curve like the Gaussian but with fatter tails. The integral
of `cp_eff - cp_base` over all T equals `L_fusion` exactly (since `f_unfrozen`
goes 0 → 1 monotonically).

This is the formulation used in CryoGrid / Westermann et al. (2011), *Modelling the
impact of wintertime rain events on the thermal regime of permafrost*, The
Cryosphere 5, 945–959; and Westermann et al. (2016), *Simulating the thermal
regime and thaw processes of ice-rich permafrost ground with the land-surface
model CryoGrid 3*, Geosci. Model Dev. 9, 523–546.

Conceptually distinct from [`Gaussian`](@ref) — there the cp shape is *imposed*;
here it *emerges* from a physical model of liquid fraction. For modelling
freeze-thaw in soils with unfrozen water content, this parametrisation is
preferred because `f_unfrozen` itself is a recognised state variable.
"""
struct WestermannSigmoid{TM,W} <: AbstractApparentHeatCapacity
    T_melt::TM
    width::W
end
WestermannSigmoid(; T_melt=0.0u"°C", width=0.2u"K") = WestermannSigmoid(T_melt, width)

@inline function apparent_heat_capacity(m::WestermannSigmoid, T, cp_base, L_fusion)
    T_c = u"°C"(T)
    W_K = ustrip(u"K", m.width)
    x = ustrip(u"K", T_c - m.T_melt) / W_K
    # sech²(x) = 1 - tanh²(x); peak at x=0 of value 1
    sech2 = 1 - tanh(x)^2
    peak = L_fusion / (2 * m.width)
    return cp_base + peak * sech2
end
