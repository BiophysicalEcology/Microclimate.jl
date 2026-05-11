# Bulk density, mineral density, and saturation moisture for the soil are
# properties of the soil profile, not of the thermal formulation — they live on
# the hydraulics model and are passed into `soil_properties` as kwargs.
@kwdef struct CampbelldeVriesSoilThermal{SF,MC,MHC,RP,RFT} <: AbstractSoilThermalModel
    de_vries_shape_factor::SF
    mineral_conductivity::MC
    mineral_heat_capacity::MHC
    recirculation_power::RP
    return_flow_threshold::RFT
end
