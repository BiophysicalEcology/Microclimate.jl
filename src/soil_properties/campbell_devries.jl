# Bulk density, mineral density, mineral conductivity, mineral heat capacity,
# and saturation moisture for the soil are properties of the soil profile, not
# of the thermal formulation — they live on the hydraulics model / soil
# profile and are passed into `soil_properties` as kwargs.
@kwdef struct CampbelldeVriesSoilProperties{SF,RP,RFT} <: AbstractSoilProperties
    de_vries_shape_factor::SF
    recirculation_power::RP
    return_flow_threshold::RFT
end

function example_soil_properties_model(;
    de_vries_shape_factor = 0.1, # de Vries shape factor, 0.33 for organic soils, 0.1 for mineral
    recirculation_power = 4.0, # power for recirculation function
    return_flow_threshold = 0.162, # return-flow cutoff soil moisture, m^3/m^3
)
    CampbelldeVriesSoilProperties(;
        de_vries_shape_factor, recirculation_power, return_flow_threshold,
    )
end

"""
    soil_properties(soil_properties_model; atmospheric_pressure, soil_temperature, soil_moisture, bulk_density, mineral_density, mineral_conductivity, mineral_heat_capacity)

Compute bulk soil properties — thermal conductivity, volumetric heat capacity,
and bulk density — for a given soil layer.

# Arguments

- `soil_properties_model::AbstractSoilProperties`

# Keywords

- `atmospheric_pressure::Quantity`: Atmospheric pressure.
- `soil_temperature::Quantity`: Soil temperature in Kelvin.
- `soil_moisture::Real`: Volumetric soil moisture (m³/m³).
- `bulk_density::Quantity`: Dry soil bulk density (kg/m³). Owned by `SoilProfile`.
- `mineral_density::Quantity`: Soil mineral density (kg/m³). Owned by `SoilProfile`.
- `mineral_conductivity::Quantity`: Soil mineral thermal conductivity (W/m/K). Owned by `SoilProfile`.
- `mineral_heat_capacity::Quantity`: Soil mineral specific heat (J/kg/K). Owned by `SoilProfile`.

# Returns

A named tuple `(bulk_thermal_conductivity, bulk_heat_capacity, bulk_density)`:
- `bulk_thermal_conductivity::Quantity`: Bulk thermal conductivity (W/m/K)
- `bulk_heat_capacity::Quantity`: Bulk mass-specific heat capacity (J/kg/K). Multiply by `bulk_density` for the volumetric form (J/m³/K) used in the surface energy balance.
- `bulk_density::Quantity`: Bulk density (kg/m³)

# Theory

This function calculates soil thermal properties using the Campbell & Norman (1991, 1994)
approach, which accounts for:

1. **Soil composition**: Soil is modeled as a mixture of mineral, water, and air fractions.
   Each fraction contributes to bulk properties depending on its volumetric fraction.

2. **Water content effects**: The function includes a water recirculation factor
   to model how moisture affects heat transfer.

3. **Phase-dependent conductivity**: Effective thermal conductivity of soil is computed
   using a generalization of de Vries' mixing model, considering interactions between
   mineral, liquid, and gas components.

4. **Temperature dependence**: Thermal conductivity of water and air are adjusted with temperature.

5. **Vapor transport**: Small contribution from vapor diffusion is included via a finite-difference
   approximation of saturated vapor pressure.

This method provides an accurate representation of heat transfer in variably wet soils
for land surface or microclimate modeling.

# References

Campbell, G. S., Jungbauer, J. D. Jr., Bidlake, W. R., & Hungerford, R. D. (1994). Predicting
 the effect of temperature on soil thermal conductivity. Soil Science, 158(5), 307–313.

Campbell, G. S., & Norman, J. M. (1998). Environmental Biophysics. Springer.

"""
function soil_properties(soil_properties_model::CampbelldeVriesSoilProperties;
    atmospheric_pressure::Quantity,
    soil_temperature::Quantity,
    soil_moisture::Number,
    bulk_density::Quantity,
    mineral_density::Quantity,
    mineral_conductivity::Quantity,
    mineral_heat_capacity::Quantity,
    vapour_pressure_equation=GoffGratch(),
)
    (; recirculation_power, return_flow_threshold, de_vries_shape_factor) = soil_properties_model

    standard_pressure = Unitful.atm
    shape_factor_c = 1.0 - 2.0 * de_vries_shape_factor

    # generalisation of eq. 8.20, Campbell and Norman (1991)
    weighting_factor(λ_component, λ_fluid) = 2.0 / (3.0 * (1.0 + de_vries_shape_factor * (λ_component / λ_fluid - 1.0))) +
                      1.0 / (3.0 * (1.0 + shape_factor_c * (λ_component / λ_fluid - 1.0)))

    temperature_celsius = ustrip(u"°C", soil_temperature)

    bulk_heat_capacity = bulk_density / mineral_density * mineral_heat_capacity + soil_moisture * water_heat_capacity
    bulk_density_total = soil_moisture * water_density + bulk_density

    # eq. 8 Campbell et al. 1994
    water_thermal_conductivity = (0.554 + 2.24e-3 * temperature_celsius - 9.87e-6 * temperature_celsius^2)u"W/m/K"
    # eq. 9 Campbell et al. 1994
    dry_air_thermal_conductivity = (0.024 + 7.73e-5 * temperature_celsius - 2.6e-8 * temperature_celsius^2)u"W/m/K"

    # p. 309 Campbell et al. 1994
    vapor_diffusivity = water_vapour_diffusivity_stp * (standard_pressure / atmospheric_pressure) * (soil_temperature / 273.15u"K")^1.75
    # p. 309 Campbell et al. 1994
    molar_density = water_vapour_molar_density_stp * (atmospheric_pressure / standard_pressure) * (273.15 / soil_temperature)
    vapor_enthalpy = molar_enthalpy_of_vaporisation(soil_temperature)

    ################################################################
    # This is some of the most expensive code in the package
    # its inlined so most of the work in wet_air_properties is ignored
    vapor_pressure = wet_air_properties(soil_temperature, 0.99, atmospheric_pressure; vapour_pressure_equation).vapour_pressure
    vapor_pressure_minus = wet_air_properties(soil_temperature - 1u"K", 0.99, atmospheric_pressure; vapour_pressure_equation).vapour_pressure
    vapor_pressure_plus = wet_air_properties(soil_temperature + 1u"K", 0.99, atmospheric_pressure; vapour_pressure_equation).vapour_pressure
    ################################################################

    vapor_pressure_gradient = (vapor_pressure_plus - vapor_pressure_minus) / 2.0

    volume_fraction_mineral = bulk_density / mineral_density
    volume_fraction_gas = max(0.0, 1.0 - soil_moisture - volume_fraction_mineral)
    # eq. 3, Campbell et al. 1994
    water_recirculation_factor = 1.0 / (1.0 + (soil_moisture / return_flow_threshold)^(-recirculation_power))

    # eq. 8.18, Campbell and Norman (1991)
    gas_thermal_conductivity = dry_air_thermal_conductivity + vapor_enthalpy * vapor_pressure_gradient * water_recirculation_factor * molar_density * vapor_diffusivity / (atmospheric_pressure - vapor_pressure)
    # eq. 8.19, Campbell and Norman (1991)
    fluid_thermal_conductivity = gas_thermal_conductivity + water_recirculation_factor * (water_thermal_conductivity - gas_thermal_conductivity)

    bulk_thermal_conductivity = (soil_moisture * weighting_factor(water_thermal_conductivity, fluid_thermal_conductivity) * water_thermal_conductivity +
           volume_fraction_mineral * weighting_factor(mineral_conductivity, fluid_thermal_conductivity) * mineral_conductivity +
           volume_fraction_gas * weighting_factor(gas_thermal_conductivity, fluid_thermal_conductivity) * gas_thermal_conductivity) /
          (soil_moisture * weighting_factor(water_thermal_conductivity, fluid_thermal_conductivity) +
           volume_fraction_mineral * weighting_factor(mineral_conductivity, fluid_thermal_conductivity) +
           volume_fraction_gas * weighting_factor(gas_thermal_conductivity, fluid_thermal_conductivity))

    return (; bulk_thermal_conductivity, bulk_heat_capacity, bulk_density=bulk_density_total)
end

# FIXME: this seems like a hack
# Per-layer field selection used by the generic `soil_properties!` loop. Each
# CampbelldeVriesSoilProperties field may be a scalar (shared) or a per-layer
# vector — `maybegetindex` rebuilds a layer-local model for the scalar
# `soil_properties` call.
maybegetindex(obj::CampbelldeVriesSoilProperties, i::Int) =
    CampbelldeVriesSoilProperties(; maybegetindex(ConstructionBase.getproperties(obj), i)...)
