"""
    MicroParameters(; soil_thermal, soil_hydraulics, snow=NoSnow())

Single home for the simulation's physical-parameter structs. Lives on
`MicroProblem.parameters`. Each field is a swappable formulation:

- `soil_thermal`: e.g. `CampbelldeVriesSoilProperties(...)`
- `soil_hydraulics`: e.g. `CampbellSoilHydraulics(...)` — owns per-depth
  `bulk_density` and `mineral_density` profiles, which the thermal model
  reads via the energy-balance plumbing
- `snow`: `NoSnow()` (default) or `SnowModel(...)`
"""
@kwdef struct MicroParameters{ST,SH,SN}
    soil_thermal::ST
    soil_hydraulics::SH
    snow::SN = NoSnow()
end
