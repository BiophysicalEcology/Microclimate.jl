# TODO are these parameters for a specific named model
@kwdef struct CampbelldeVriesSoilThermal{SF,MC,MD,MHC,BD,SM,RP,RFT} <: AbstractSoilThermalModel
    de_vries_shape_factor::SF
    mineral_conductivity::MC
    mineral_density::MD
    mineral_heat_capacity::MHC
    bulk_density::BD
    saturation_moisture::SM
    recirculation_power::RP
    return_flow_threshold::RFT
end
