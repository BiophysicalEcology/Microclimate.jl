using Microclimate
using Unitful
using FluidProperties
using Test

soil_mineral_conductivity = 1.25u"W/m/K" # soil minerals thermal conductivity
soil_mineral_density = 2.560u"Mg/m^3" # soil minerals density
soil_mineral_heat_capacity = 870.0u"J/kg/K" # soil minerals specific heat
soil_bulk_density = 1.3u"Mg/m^3" # dry soil bulk density
soil_saturation_moisture = 0.26u"m^3/m^3" # volumetric water content at saturation (0.1 bar matric potential)

elevation = 0u"m"
P_atmos = atmospheric_pressure(elevation)

T_soil = (collect(-10.0:80).+273.15)u"K"
N = length(T_soil)
θ_soil = fill(0.2, N)

soilprops = (;
    ρ_dry=fill(soil_bulk_density, N),
    θ_sat=fill(soil_saturation_moisture, N),
    λ_mineral=fill(soil_mineral_conductivity, N),
    cp_mineral=fill(soil_mineral_heat_capacity, N),
    ρ_mineral=fill(soil_mineral_density, N),
)

# soilprops = (;
#     ρ_dry=soil_bulk_density,
#     θ_sat=soil_saturation_moisture,
#     λ_mineral=soil_mineral_conductivity,
#     cp_mineral=soil_mineral_heat_capacity,
#     ρ_mineral=soil_mineral_density,
# )

(; ρ_dry, θ_sat, λ_mineral, cp_mineral, ρ_mineral) = soilprops

λ_b, c_p_b, ρ_b = soil_props(;
    T_soil=T_soil[1], θ_soil=θ_soil[1], soilprops=soilprops, elevation, P_atmos
)

λ_b, cp_b, ρ_b = soil_props_vector(T_soil, θ_soil, soilprops, elevation, P_atmos)