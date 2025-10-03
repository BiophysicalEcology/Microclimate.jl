using Microclimate
using Unitful
using FluidProperties
using Test

soil_mineral_conductivity = 1.25u"W/m/K" # soil minerals thermal conductivity
soil_mineral_density = 2.560u"Mg/m^3" # soil minerals density
soil_mineral_heat_capacity = 870.0u"J/kg/K" # soil minerals specific heat
soil_bulk_density = 1.3u"Mg/m^3" # dry soil bulk density
recirculation_power = 4.0 # power for recirculation function
return_flow_threshold = 0.162 # return-flow cutoff soil moisture, m^3/m^3

elevation = 0u"m"
P_atmos = atmospheric_pressure(elevation)


T_soil = (collect(-10.0:80).+273.15)u"K"
N = length(T_soil)
θ_soil = fill(0.2, N)

soilprops = (;
    ρ_dry=soil_bulk_density,
    λ_mineral=soil_mineral_conductivity,
    cp_mineral=soil_mineral_heat_capacity,
    ρ_mineral=soil_mineral_density,
    q=recirculation_power,
    θ_0=return_flow_threshold,
)

(; ρ_dry, λ_mineral, cp_mineral, ρ_mineral, q, θ_0) = soilprops

λ_b, c_p_b, ρ_b = soil_props(;
    T_soil=T_soil[1], θ_soil=θ_soil[1], soilprops=soilprops, elevation, P_atmos
)

soilprops = (;
    ρ_dry=fill(soil_bulk_density, N),
    λ_mineral=fill(soil_mineral_conductivity, N),
    cp_mineral=fill(soil_mineral_heat_capacity, N),
    ρ_mineral=fill(soil_mineral_density, N),
    q=fill(recirculation_power, N),
    θ_0=fill(return_flow_threshold, N),
)

soil_properties_buffers = allocate_soil_properties(N, soilprops)

λ_b, cp_b, ρ_b = soil_props_vector(soil_properties_buffers; T_soil, θ_soil, soilprops, elevation, P_atmos)