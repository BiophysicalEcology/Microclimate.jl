# Julia translation of FORTRAN code that calculates variable thermal conductivity and specific heat
# for soil and snow layers, based on Campbell et al. (1994) and Campbell & Norman (1998)

function soil_props(; 
    T_soil::Q, 
    θ_soil::R, 
    soilprops::S, 
    elevation::E, 
    P_atmos = atmospheric_pressure(elevation),
) where {Q<:Quantity, R<:Real, S<:NamedTuple, E<:Quantity}

    (; ρ_dry, λ_mineral, cp_mineral, ρ_mineral) = soilprops

    p_a0 = Unitful.atm
    q = 4.0
    θ_0 = 0.162
    g_a = 0.1
    g_c = 1.0 - 2.0 * g_a

    ϵ(λ_λ, λ_fluid) = 2.0 / (3.0 * (1.0 + g_a * (λ_λ / λ_fluid - 1.0))) + 
                      1.0 / (3.0 * (1.0 + g_c * (λ_λ / λ_fluid - 1.0)))

    T_K = T_soil
    T_C = ustrip(u"°C", T_K)

    cp_b = ρ_dry / ρ_mineral * cp_mineral + θ_soil * cp_water
    ρ_b  = θ_soil * ρ_water + ρ_dry

    λ_water   = (0.554 + 2.24e-3 * T_C - 9.87e-6 * T_C^2)u"W/m/K"
    λ_dry_air = (0.024 + 7.73e-5 * T_C - 2.6e-8 * T_C^2)u"W/m/K"

    D_v    = D_v0 * (p_a0 / P_atmos) * (T_K / 273.15u"K")^1.75
    ρ_hat  = ρ_hat0 * (P_atmos / p_a0) * (273.15 / T_K)
    λ_vapor = molar_enthalpy_of_vaporisation(T_K)

    e_a  = wet_air_properties(T_K; rh=99.0, P_atmos).P_vap
    e_a1 = wet_air_properties(T_K - 1u"K"; rh=99.0, P_atmos).P_vap
    e_a2 = wet_air_properties(T_K + 1u"K"; rh=99.0, P_atmos).P_vap

    ∇x = (e_a2 - e_a1) / 2.0

    ϕ_mineral = ρ_dry / ρ_mineral
    ϕ_gas     = max(0.0, 1.0 - θ_soil - ϕ_mineral)
    f_water   = 1.0 / (1.0 + (θ_soil / θ_0)^(-q))

    λ_gas   = λ_dry_air + λ_vapor * ∇x * f_water * ρ_hat * D_v / (P_atmos - e_a)
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
function soil_props_vector(T_soil::AbstractVector, θ_soil::AbstractVector, soilprops::NamedTuple, elevation, P_atmos)
    N = length(T_soil)
    @assert length(θ_soil) == N

    soil_props_i(i) = soil_props(
        T_soil = T_soil[i],
        θ_soil = θ_soil[i],
        soilprops = (
            ρ_dry     = soilprops.ρ_dry[i],
            λ_mineral = soilprops.λ_mineral[i],
            cp_mineral = soilprops.cp_mineral[i],
            ρ_mineral = soilprops.ρ_mineral[i],
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

