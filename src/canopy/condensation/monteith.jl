"""
    MonteithLeafCondensation(; gamma=66.0u"Pa/K")

Leaf dew/frost from Monteith, J. L. (1963). Dew: facts and fallacies. In
The Water Relations of Plants (pp. 37-56). Blackwell. Adapted from the
author's own validated R port (`calc_dew_object.R`), substituting a leaf's
own heat-balance terms for an animal's. A single-term combination equation
-- net radiative loss weighted by the saturation-vapour-pressure-curve
slope, plus an aerodynamic term using the leaf's own heat transfer
coefficient directly (Lewis number ≈ 1, no separate bulk transfer
coefficient). No ground-heat-flux analogue: a thin leaf has negligible heat
storage, matching Garratt & Segal (1988) eq. 6a's own omission of `G` for
canopy dew.

Simpler than [`GarrattSegalCondensation`](@ref): `gamma` is a fixed
psychrometric constant rather than derived each hour, and the slope is
vapour-pressure- rather than specific-humidity-based.
"""
@kwdef struct MonteithLeafCondensation{G} <: AbstractCondensationModel
    gamma::G = 66.0u"Pa/K"
end

function condensation_energy_flux(model::MonteithLeafCondensation;
    leaf_temperature, air_temperature, relative_humidity,
    absorbed_radiation, leaf_emissivity, heat_transfer_coefficient,
    leaf_area, aerodynamic_area,
    vapour_pressure_equation=GoffGratch(), kw...,
)
    T_l = u"K"(leaf_temperature)
    T_1 = u"K"(air_temperature)

    saturated_vapour_pressure_air = vapour_pressure(vapour_pressure_equation, T_1)
    δ = (vapour_pressure(vapour_pressure_equation, T_1 + 1.0u"K") - saturated_vapour_pressure_air) / 1.0u"K"
    γ = model.gamma

    # absorbed_radiation is per unit leaf_area (radiative, self-shaded);
    # the aerodynamic term uses aerodynamic_area (PAI-based, convective) --
    # same split leaf_heat_balance uses. Totalled to Watts then renormalized
    # to ground area, matching sensible_heat_source's own convention.
    emitted_longwave = leaf_emissivity * σ * T_l^4 * leaf_area
    absorbed_total = absorbed_radiation * leaf_area
    R_N = emitted_longwave - absorbed_total  # Monteith convention: +ve = net loss, Watts
    vapour_pressure_deficit = saturated_vapour_pressure_air * (1.0 - relative_humidity)
    aerodynamic_term = heat_transfer_coefficient * aerodynamic_area * vapour_pressure_deficit  # Watts

    net_loss = (δ / (δ + γ)) * R_N - aerodynamic_term / (δ + γ)  # Watts
    return -net_loss / 1.0u"m^2"  # this interface's convention: +ve loss, -ve condensation, per ground area
end
