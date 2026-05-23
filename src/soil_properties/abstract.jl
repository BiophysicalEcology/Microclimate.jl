abstract type AbstractSoilProperties end

"""
    soil_properties(soil_properties_model; atmospheric_pressure, soil_temperature, soil_moisture, bulk_density, mineral_density, vapour_pressure_equation)

Scalar bulk-property computation for one soil layer. Every variant must
accept the listed kwargs and return
`(; bulk_thermal_conductivity, bulk_heat_capacity, bulk_density)` in
W/m/K, J/kg/K, and a density unit respectively. The vectorized
`soil_properties!` loop and the layer-by-layer ODE pipeline rely on that
contract.
"""
function soil_properties end

# Generic vectorized loop over layers — calls the variant's scalar
# `soil_properties` method and the variant's `maybegetindex` accessor.
function soil_properties!(buffers::NamedTuple, soil_properties_model::AbstractSoilProperties;
    atmospheric_pressure::Quantity, soil_temperature::AbstractVector, soil_moisture::AbstractVector,
    bulk_density::AbstractVector,
    mineral_density::AbstractVector,
    vapour_pressure_equation=GoffGratch(),
)
    num_layers = length(soil_temperature)
    @assert length(soil_moisture) == num_layers
    (; bulk_thermal_conductivity, bulk_heat_capacity) = buffers
    out_bulk_density = buffers.bulk_density
    bd_unit = unit(eltype(out_bulk_density))
    @inbounds for i in 1:num_layers
        result = soil_properties(maybegetindex(soil_properties_model, i);
            atmospheric_pressure,
            soil_temperature = soil_temperature[i],
            soil_moisture = soil_moisture[i],
            bulk_density = bulk_density[i],
            mineral_density = mineral_density[i],
            vapour_pressure_equation,
        )
        bulk_thermal_conductivity[i] = result.bulk_thermal_conductivity
        bulk_heat_capacity[i]        = result.bulk_heat_capacity
        out_bulk_density[i]          = uconvert(bd_unit, result.bulk_density)
    end

    return (; bulk_thermal_conductivity, bulk_heat_capacity, bulk_density=out_bulk_density)
end

"""
    allocate_soil_properties(nodes, soil_properties_model, soil_profile)

Allocate the per-layer output buffers for `soil_properties!`. Default
implementation produces the three buffers required by the interface
contract; variants that need additional per-layer storage can override.
"""
function allocate_soil_properties(nodes, ::AbstractSoilProperties, soil_profile)
    num_nodes = length(nodes)
    return (;
        bulk_thermal_conductivity = zeros(typeof(0.0u"W/m/K"), num_nodes),
        bulk_heat_capacity = zeros(typeof(0.0u"J/kg/K"), num_nodes),
        bulk_density = zeros(eltype(soil_profile.mineral_density), num_nodes),
    )
end
