"""
    CampbellSoilHydraulics(; ..., mode=PrescribedSoilMoisture())

Soil hydraulic parameters for Campbell's (1985) soil water balance model.
The `mode` field controls how soil moisture is handled during simulation:

- `PrescribedSoilMoisture()` — use prescribed wetness from environment data (default)
- `DynamicSoilMoisture()` — run the full Campbell soil water balance solver

# References
Campbell, G. S. (1985). Soil Physics with BASIC. Elsevier.
"""
@kwdef struct CampbellSoilHydraulics{AEWP,SHC,CBP,BD,MD,RDen,RRes,SCP,LRes,SSP,RRad,ME,MC,MS,MP,Mode<:AbstractSoilMoistureMode} <: AbstractSoilMoistureModel
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
    moist_error::ME
    moist_count::MC
    moist_step::MS
    maxpool::MP
    mode::Mode = PrescribedSoilMoisture()
end
