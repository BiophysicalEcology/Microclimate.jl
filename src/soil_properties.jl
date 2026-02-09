"""
    soil_properties(soil_thermal; soil_temperature, soil_moisture, micro_terrain)

Compute bulk soil properties — thermal conductivity (`λ_b`), volumetric heat capacity (`cp_b`), 
and bulk density (`ρ_b`) — for a given soil layer.

# Arguments

- `soil_thermal::AbstractSoilThermalModel`

# Keywords

- `micro_terrain`
- `soil_temperature::Quantity`: Soil temperature in Kelvin.
- `soil_moisture::Real`: Volumetric soil moisture (m³/m³).

# Returns

A named tuple `(λ_b, cp_b, ρ_b)`:
- `λ_b::Quantity`: Bulk thermal conductivity (W/m/K)
- `cp_b::Quantity`: Bulk volumetric heat capacity (J/kg/K)
- `ρ_b::Quantity`: Bulk density (kg/m³)

# Theory

This function calculates soil thermal properties using the Campbell & Norman (1991, 1994) 
approach, which accounts for:

1. **Soil composition**: Soil is modeled as a mixture of mineral, water, and air fractions. 
   Each fraction contributes to bulk properties depending on its volumetric fraction.

2. **Water content effects**: The function includes a water recirculation factor (`f_water`) 
   to model how moisture affects heat transfer.

3. **Phase-dependent conductivity**: Effective thermal conductivity of soil (`λ_b`) is computed 
   using a generalization of de Vries’ mixing model, considering interactions between 
   mineral, liquid, and gas components.

4. **Temperature dependence**: Thermal conductivity of water and air are adjusted with temperature.

5. **Vapor transport**: Small contribution from vapor diffusion is included via a finite-difference 
   approximation of saturated vapor pressure (`∇x`).

This method provides an accurate representation of heat transfer in variably wet soils 
for land surface or microclimate modeling.

# References

Campbell, G. S., Jungbauer, J. D. Jr., Bidlake, W. R., & Hungerford, R. D. (1994). Predicting 
 the effect of temperature on soil thermal conductivity. Soil Science, 158(5), 307–313.

Campbell, G. S., & Norman, J. M. (1998). Environmental Biophysics. Springer.

"""
function soil_properties(soil_thermal::CampbelldeVriesSoilThermal;
    P_atmos::Quantity,
    soil_temperature::Quantity,
    soil_moisture::Number,
)
    # Soil thermal parameters
    ρ_dry = soil_thermal.bulk_density
    λ_mineral = soil_thermal.mineral_conductivity
    cp_mineral = soil_thermal.mineral_heat_capacity
    ρ_mineral = soil_thermal.mineral_density
    recirculation_power = soil_thermal.recirculation_power
    return_flow_threshold = soil_thermal.return_flow_threshold
    shape_factor_a = soil_thermal.deVries_shape_factor

    θ_soil = soil_moisture

    standard_pressure = Unitful.atm
    shape_factor_c = 1.0 - 2.0 * shape_factor_a

    # generalisation of eq. 8.20, Campbell and Norman (1991)
    weighting_factor(λ_component, λ_fluid) = 2.0 / (3.0 * (1.0 + shape_factor_a * (λ_component / λ_fluid - 1.0))) +
                      1.0 / (3.0 * (1.0 + shape_factor_c * (λ_component / λ_fluid - 1.0)))

    temperature_kelvin = soil_temperature
    temperature_celsius = ustrip(u"°C", temperature_kelvin)

    cp_b = ρ_dry / ρ_mineral * cp_mineral + θ_soil * cp_water
    ρ_b = θ_soil * ρ_water + ρ_dry

    # eq. 8 Campbell et al. 1994
    λ_water = (0.554 + 2.24e-3 * temperature_celsius - 9.87e-6 * temperature_celsius^2)u"W/m/K"
    # eq. 9 Campbell et al. 1994
    λ_dry_air = (0.024 + 7.73e-5 * temperature_celsius - 2.6e-8 * temperature_celsius^2)u"W/m/K"

    # p. 309 Campbell et al. 1994
    vapor_diffusivity = D_v0 * (standard_pressure / P_atmos) * (temperature_kelvin / 273.15u"K")^1.75
    # p. 309 Campbell et al. 1994
    molar_density = ρ_hat0 * (P_atmos / standard_pressure) * (273.15 / temperature_kelvin)
    λ_vapor = molar_enthalpy_of_vaporisation(temperature_kelvin)

    ################################################################
    # This is some of the most expensive code in the package
    # its inlined so most of the work in wet_air_properties is ignored
    vapor_pressure = wet_air_properties(temperature_kelvin, 0.99, P_atmos).P_vap
    vapor_pressure_minus = wet_air_properties(temperature_kelvin - 1u"K", 0.99, P_atmos).P_vap
    vapor_pressure_plus = wet_air_properties(temperature_kelvin + 1u"K", 0.99, P_atmos).P_vap
    ################################################################

    vapor_pressure_gradient = (vapor_pressure_plus - vapor_pressure_minus) / 2.0

    volume_fraction_mineral = ρ_dry / ρ_mineral
    volume_fraction_gas = max(0.0, 1.0 - θ_soil - volume_fraction_mineral)
    # eq. 3, Campbell et al. 1994
    water_recirculation_factor = 1.0 / (1.0 + (θ_soil / return_flow_threshold)^(-recirculation_power))

    # eq. 8.18, Campbell and Norman (1991)
    λ_gas = λ_dry_air + λ_vapor * vapor_pressure_gradient * water_recirculation_factor * molar_density * vapor_diffusivity / (P_atmos - vapor_pressure)
    # eq. 8.19, Campbell and Norman (1991)
    λ_fluid = λ_gas + water_recirculation_factor * (λ_water - λ_gas)

    λ_b = (θ_soil * weighting_factor(λ_water, λ_fluid) * λ_water +
           volume_fraction_mineral * weighting_factor(λ_mineral, λ_fluid) * λ_mineral +
           volume_fraction_gas * weighting_factor(λ_gas, λ_fluid) * λ_gas) /
          (θ_soil * weighting_factor(λ_water, λ_fluid) +
           volume_fraction_mineral * weighting_factor(λ_mineral, λ_fluid) +
           volume_fraction_gas * weighting_factor(λ_gas, λ_fluid))

    return (; λ_b, cp_b, ρ_b)
end

function allocate_soil_properties(nodes, soil_thermal)
    (; mineral_conductivity, mineral_heat_capacity, mineral_density) = soil_thermal
    num_nodes = length(nodes)

    λ_b = fill(mineral_conductivity[1], num_nodes)
    cp_b = fill(mineral_heat_capacity[1], num_nodes)
    ρ_b = fill(mineral_density[1], num_nodes)

    return (; λ_b, cp_b, ρ_b)
end

"""
    soil_props_vector!(buffers, soil_thermal; micro_micro_terrain, soil_temperature, soil_moisture)

Compute soil properties for vectors of soil temperature and moisture using broadcasting.

# TODO these should have readable names
Returns three arrays: `λ_b`, `cp_b`, `ρ_b`.
"""
function soil_properties!(buffers::NamedTuple, soil_thermal;
    P_atmos::Quantity, soil_temperature::AbstractVector, soil_moisture::AbstractVector
)
    num_layers = length(soil_temperature)
    @assert length(soil_moisture) == num_layers
    (; λ_b, cp_b, ρ_b) = buffers
    soil_props_i(i) = soil_properties(soil_thermal;
        P_atmos,
        soil_temperature = soil_temperature[i],
        soil_moisture = soil_moisture[i],
    )

    results = soil_props_i.(1:num_layers)

    λ_b .= getindex.(results, 1)
    cp_b .= getindex.(results, 2)
    ρ_b .= getindex.(results, 3)

    return (; λ_b, cp_b, ρ_b)
end
