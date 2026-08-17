"""
    PrescribedStomatalConductance(; conductance=LeafEvaporationParameters(...),
                                    closed_conductance=<conductance with zero stomatal terms>)

Stomata open to `conductance` in daylight (`zenith_angle < 90°`) and shut to
`closed_conductance` (cuticular-only, by default) at night. Not coupled to
soil moisture or photosynthesis — see [`AbstractStomatalConductanceModel`](@ref)
for where those variants plug in.

`conductance` defaults are free/tunable, uncited (same values
`MoistureResponsiveStomatalConductance` uses).
"""
@kwdef struct PrescribedStomatalConductance{C,CC} <: AbstractStomatalConductanceModel
    conductance::C = LeafEvaporationParameters(;
        abaxial_vapour_conductance=0.3u"mol/m^2/s", adaxial_vapour_conductance=0.0u"mol/m^2/s",
        cuticular_conductance=0.01u"mol/m^2/s",
    )
    closed_conductance::CC = LeafEvaporationParameters(;
        abaxial_vapour_conductance=0.0u"mol/m^2/s", adaxial_vapour_conductance=0.0u"mol/m^2/s",
        cuticular_conductance=conductance.cuticular_conductance,
    )
end

stomatal_conductance(model::PrescribedStomatalConductance, zenith_angle, leaf_water_potential, soil_hydraulic_model) =
    zenith_angle < 90.0u"°" ? model.conductance : model.closed_conductance
