"""
    MoistureResponsiveStomatalConductance(; conductance=LeafEvaporationParameters(...))

Stomata close as leaf water potential drops, using the `stomatal_closure_potential` and
`stomatal_stability_parameter` of the `CampbellSoilHydraulics` passed to
[`stomatal_conductance`](@ref). Cuticular conductance is unaffected. Shut at night
(`zenith_angle >= 90°`). Use with `CampbellSoilHydraulics` so soil-water stress reaches the leaves.

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
