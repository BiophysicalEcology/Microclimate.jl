"""
    MoistureResponsiveStomatalConductance(; conductance=LeafEvaporationParameters(...))

Stomatal (abaxial/adaxial) conductance closes as leaf water potential drops,
using `CampbellSoilHydraulics`'s own `stomatal_closure_potential`/
`stomatal_stability_parameter` (read from the `soil_hydraulic_model` passed
to [`stomatal_conductance`](@ref) — not stored here, so the two models can't
drift out of sync). Cuticular conductance is unaffected. Shuts fully at
night (`zenith_angle >= 90°`).

`leaf_water_potential` is supplied by the caller each call (e.g. from
`CampbellSoilHydraulics`'s own solve). This is the stomatal_model to pair
with `CampbellSoilHydraulics` on a `MultilayerCanopy` if soil-moisture
stress should reach the leaf's transpiration/temperature — pairing with any
other `soil_hydraulic_model` is a `MethodError` by construction.

# References
Campbell, G. S. 1985. “Soil Physics with Basic: Transport Models for Soil-Plant
Systems.” Elsevier.
"""
@kwdef struct MoistureResponsiveStomatalConductance{C} <: AbstractStomatalConductanceModel
    conductance::C = LeafEvaporationParameters(;
        abaxial_vapour_conductance=0.3u"mol/m^2/s", adaxial_vapour_conductance=0.0u"mol/m^2/s",
        cuticular_conductance=0.01u"mol/m^2/s",
    )
end

function stomatal_conductance(model::MoistureResponsiveStomatalConductance, zenith_angle, leaf_water_potential,
    soil_hydraulic_model::CampbellSoilHydraulics,
)
    (; conductance) = model
    (; stomatal_closure_potential, stomatal_stability_parameter) = soil_hydraulic_model
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
