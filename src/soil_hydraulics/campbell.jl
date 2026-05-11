"""
    CampbellSoilHydraulics(; ...)

Soil hydraulic parameters for Campbell's (1985) soil water balance model.
Holds parameters only — the moisture-solver strategy (`PrescribedSoilMoisture` /
`DynamicSoilMoisture`) and tuning (`moisture_tolerance`, `moisture_max_iterations`,
`moisture_timestep`, `max_surface_pool`) live on `MicroConfig`, since they are
solver/strategy choices, not Campbell-specific.

# References
Campbell, G. S. (1985). Soil Physics with BASIC. Elsevier.
"""
@kwdef struct CampbellSoilHydraulics{AEWP,SHC,CBP,BD,MD,RDen,RRes,SCP,LRes,SSP,RRad} <: AbstractSoilHydraulicsModel
    air_entry_water_potential::AEWP
    saturated_hydraulic_conductivity::SHC
    campbell_b_parameter::CBP
    bulk_density::BD       # per-depth profile (Vector); also consumed by the soil thermal model
    mineral_density::MD    # per-depth profile (Vector); also consumed by the soil thermal model
    root_density::RDen
    root_resistance::RRes
    stomatal_closure_potential::SCP
    leaf_resistance::LRes
    stomatal_stability_parameter::SSP
    root_radius::RRad
end
