function allocate_soil_properties(nodes, soilprops)
    (; λ_mineral, cp_mineral, ρ_mineral) = soilprops
    NON = length(nodes)

    λ_b = fill(λ_mineral[1], NON)
    cp_b = fill(cp_mineral[1], NON)
    ρ_b = fill(ρ_mineral[1], NON)

    return (; λ_b, cp_b, ρ_b)
end

"""
    soil_props(; T_soil, θ_soil, soilprops, elevation, P_atmos=atmospheric_pressure(elevation))

Compute bulk soil properties — thermal conductivity (`λ_b`), volumetric heat capacity (`cp_b`), 
and bulk density (`ρ_b`) — for a given soil layer.

# Arguments
- `T_soil::Quantity`: Soil temperature in Kelvin.
- `θ_soil::Real`: Volumetric soil moisture (m³/m³).
- `soilprops::NamedTuple`: Soil physical parameters, with components:
    - `ρ_dry`: Dry soil bulk density (kg/m³)
    - `λ_mineral`: Thermal conductivity of mineral fraction (W/m/K)
    - `cp_mineral`: Heat capacity of mineral fraction (J/kg/K)
    - `ρ_mineral`: Density of mineral fraction (kg/m³)
    - `q`: Empirical parameter for water recirculation (dimensionless), values between
       ~2 and 6, typically lower for coarser particles
    - `θ_0`: Reference water content for recirculation (m³/m³), ~0.05 for coarse sand 
       to 0.25 for heavy clay
    - `g_a`: de Vries shape factor, 0.33 for organic soils, 0.1 for mineral
- `elevation::Quantity`: Elevation above sea level.
- `P_atmos::Quantity` (optional): Atmospheric pressure; defaults to `atmospheric_pressure(elevation)`.

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
function soil_props(; 
    T_soil::T, # vector of soil temperatures, K
    θ_soil::R, # vector of soil moisture, m^3/m^3
    soilprops::S, # soil properties
    elevation::E, # elevation, m
    P_atmos = atmospheric_pressure(elevation),
) where {T<:Quantity, R<:Real, S<:NamedTuple, E<:Quantity}

    # ρ = density, λ = thermal conductivity, cp = specific heat capacity
    (; ρ_dry, λ_mineral, cp_mineral, ρ_mineral, q, θ_0, g_a) = soilprops

    p_a0 = Unitful.atm
    g_c = 1.0 - 2.0 * g_a

    # generalisation of eq. 8.20, Campbell and Norman (1991)
    ϵ(λ_λ, λ_fluid) = 2.0 / (3.0 * (1.0 + g_a * (λ_λ / λ_fluid - 1.0))) +
                      1.0 / (3.0 * (1.0 + g_c * (λ_λ / λ_fluid - 1.0)))

    T_K = T_soil
    T_C = ustrip(u"°C", T_K)

    cp_b = ρ_dry / ρ_mineral * cp_mineral + θ_soil * cp_water
    ρ_b  = θ_soil * ρ_water + ρ_dry

    # eq. 8 Campbell et al. 1994
    λ_water   = (0.554 + 2.24e-3 * T_C - 9.87e-6 * T_C^2)u"W/m/K"
    # eq. 9 Campbell et al. 1994
    λ_dry_air = (0.024 + 7.73e-5 * T_C - 2.6e-8 * T_C^2)u"W/m/K"

    # p. 309 Campbell et al. 1994
    D_v    = D_v0 * (p_a0 / P_atmos) * (T_K / 273.15u"K")^1.75
    # p. 309 Campbell et al. 1994
    ρ_hat  = ρ_hat0 * (P_atmos / p_a0) * (273.15 / T_K)
    λ_vapor = molar_enthalpy_of_vaporisation(T_K)

    ################################################################ 
    # This is some of the most expensive code in the package 
    # its inlined so most of the work in wet_air_properties is ignored 
    e_a = wet_air_properties(T_K; rh=99.0, P_atmos).P_vap 
    e_a1 = wet_air_properties(T_K - 1u"K"; rh=99.0, P_atmos).P_vap 
    e_a2 = wet_air_properties(T_K + 1u"K"; rh=99.0, P_atmos).P_vap 
    ################################################################

    ∇x = (e_a2 - e_a1) / 2.0

    ϕ_mineral = ρ_dry / ρ_mineral
    ϕ_gas     = max(0.0, 1.0 - θ_soil - ϕ_mineral)
    # eq. 3, Campbell et al. 1994
    f_water   = 1.0 / (1.0 + (θ_soil / θ_0)^(-q))
    
    # eq. 8.18, Campell and Norman (1991)
    λ_gas   = λ_dry_air + λ_vapor * ∇x * f_water * ρ_hat * D_v / (P_atmos - e_a)
    # eq. 8.19, Campell and Norman (1991)
    λ_fluid = λ_gas + f_water * (λ_water - λ_gas)

    λ_b = (θ_soil * ϵ(λ_water, λ_fluid) * λ_water +
           ϕ_mineral * ϵ(λ_mineral, λ_fluid) * λ_mineral +
           ϕ_gas * ϵ(λ_gas, λ_fluid) * λ_gas) /
          (θ_soil * ϵ(λ_water, λ_fluid) + 
           ϕ_mineral * ϵ(λ_mineral, λ_fluid) + 
           ϕ_gas * ϵ(λ_gas, λ_fluid))

    return (; λ_b, cp_b, ρ_b)
end


"""
    soil_props_vector(T_soil, θ_soil, soilprops, elevation)

Compute soil properties for vectors of soil temperature and moisture using broadcasting.
Returns three arrays: `λ_b`, `cp_b`, `ρ_b`.
"""
function soil_props_vector(buffers::NamedTuple; T_soil::AbstractVector, θ_soil::AbstractVector, soilprops::NamedTuple, elevation, P_atmos)
    N = length(T_soil)
    @assert length(θ_soil) == N
    (; λ_b, cp_b, ρ_b) = buffers
    soil_props_i(i) = soil_props(
        T_soil = T_soil[i],
        θ_soil = θ_soil[i],
        soilprops = (
            ρ_dry     = soilprops.ρ_dry[i],
            λ_mineral = soilprops.λ_mineral[i],
            cp_mineral = soilprops.cp_mineral[i],
            ρ_mineral = soilprops.ρ_mineral[i],
            q = soilprops.q[i],
            θ_0 = soilprops.θ_0[i],
            g_a = soilprops.g_a[i],
        ),
        elevation = elevation,
        P_atmos = P_atmos
    )

    results = soil_props_i.(1:N)

    λ_b  = getindex.(results, 1)
    cp_b = getindex.(results, 2)
    ρ_b  = getindex.(results, 3)

    return λ_b, cp_b, ρ_b
end

