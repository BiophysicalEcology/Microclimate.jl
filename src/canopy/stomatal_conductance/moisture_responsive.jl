"""
    MoistureResponsiveStomatalConductance(; conductance=LeafEvaporationParameters(...),
                                             stomatal_closure_potential=-1500.0u"J/kg",
                                             stomatal_stability_parameter=10.0)

Stomatal (abaxial/adaxial) conductance closes as leaf water potential drops,
using the same closure curve as `CampbellSoilHydraulics`
(`soil_hydraulics/campbell.jl`, EQ12.28): with
`stomatal_closure_factor = (leaf_water_potential / stomatal_closure_potential)^stomatal_stability_parameter`,
the open fraction is `1 / (1 + stomatal_closure_factor)`. Cuticular
conductance is unaffected (it isn't stomatally controlled). Also shuts fully
at night (`zenith_angle >= 90°`), matching `PrescribedStomatalConductance`.

`leaf_water_potential` is supplied by the caller each call (e.g. from
`CampbellSoilHydraulics`'s own leaf-water-potential solve) — see
[`AbstractStomatalConductanceModel`](@ref).
"""
@kwdef struct MoistureResponsiveStomatalConductance{C,SCP,SSP} <: AbstractStomatalConductanceModel
    conductance::C = LeafEvaporationParameters(;
        abaxial_vapour_conductance=0.3u"mol/m^2/s", adaxial_vapour_conductance=0.0u"mol/m^2/s",
        cuticular_conductance=0.01u"mol/m^2/s",
    )
    stomatal_closure_potential::SCP = -1500.0u"J/kg"
    stomatal_stability_parameter::SSP = 10.0
end

function stomatal_conductance(model::MoistureResponsiveStomatalConductance, zenith_angle, leaf_water_potential)
    (; conductance, stomatal_closure_potential, stomatal_stability_parameter) = model
    if zenith_angle >= 90.0u"°"
        return LeafEvaporationParameters(;
            abaxial_vapour_conductance=0.0u"mol/m^2/s", adaxial_vapour_conductance=0.0u"mol/m^2/s",
            cuticular_conductance=conductance.cuticular_conductance,
        )
    end
    closure_ratio = max(leaf_water_potential / stomatal_closure_potential, 0.0)
    stomatal_closure_factor = closure_ratio^stomatal_stability_parameter
    fraction_open = 1.0 / (1.0 + stomatal_closure_factor)
    return LeafEvaporationParameters(;
        abaxial_vapour_conductance=conductance.abaxial_vapour_conductance * fraction_open,
        adaxial_vapour_conductance=conductance.adaxial_vapour_conductance * fraction_open,
        cuticular_conductance=conductance.cuticular_conductance,
    )
end
