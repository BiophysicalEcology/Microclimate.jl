"""
    soil_properties(soil_thermal; atmospheric_pressure, soil_temperature, soil_moisture)

Compute bulk soil properties — thermal conductivity, volumetric heat capacity,
and bulk density — for a given soil layer.

# Arguments

- `soil_thermal::AbstractSoilThermalModel`

# Keywords

- `atmospheric_pressure::Quantity`: Atmospheric pressure.
- `soil_temperature::Quantity`: Soil temperature in Kelvin.
- `soil_moisture::Real`: Volumetric soil moisture (m³/m³).

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
function soil_properties(soil_thermal::CampbelldeVriesSoilThermal;
    atmospheric_pressure::Quantity,
    soil_temperature::Quantity,
    soil_moisture::Number,
    vapour_pressure_equation=GoffGratch(),
)
    (; bulk_density, mineral_conductivity, mineral_heat_capacity, mineral_density,
       recirculation_power, return_flow_threshold, de_vries_shape_factor) = soil_thermal

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

function allocate_soil_properties(nodes, soil_thermal)
    (; mineral_conductivity, mineral_heat_capacity, mineral_density) = soil_thermal
    num_nodes = length(nodes)

    bulk_thermal_conductivity = fill(mineral_conductivity[1], num_nodes)
    bulk_heat_capacity = fill(mineral_heat_capacity[1], num_nodes)
    bulk_density = fill(mineral_density[1], num_nodes)

    return (; bulk_thermal_conductivity, bulk_heat_capacity, bulk_density)
end

"""
    soil_properties!(buffers, soil_thermal; atmospheric_pressure, soil_temperature, soil_moisture)

Compute soil properties for vectors of soil temperature and moisture using broadcasting.

Returns three arrays: `bulk_thermal_conductivity`, `bulk_heat_capacity`, `bulk_density`.
"""
function soil_properties!(buffers::NamedTuple, soil_thermal;
    atmospheric_pressure::Quantity, soil_temperature::AbstractVector, soil_moisture::AbstractVector,
    vapour_pressure_equation=GoffGratch(),
)
    num_layers = length(soil_temperature)
    @assert length(soil_moisture) == num_layers
    (; bulk_thermal_conductivity, bulk_heat_capacity, bulk_density) = buffers
    soil_props_i(i) = soil_properties(maybegetindex(soil_thermal, i);
        atmospheric_pressure,
        soil_temperature = soil_temperature[i],
        soil_moisture = soil_moisture[i],
        vapour_pressure_equation,
    )

    results = soil_props_i.(1:num_layers)

    bulk_thermal_conductivity .= getindex.(results, 1)
    bulk_heat_capacity .= getindex.(results, 2)
    bulk_density .= uconvert.(unit(bulk_density[1]), getindex.(results, 3))

    return (; bulk_thermal_conductivity, bulk_heat_capacity, bulk_density)
end
