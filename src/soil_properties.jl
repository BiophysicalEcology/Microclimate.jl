# Julia translation of FORTRAN code that calculates variable thermal conductivity and specific heat
# for soil and snow layers, based on Campbell et al. (1994) and Campbell & Norman (1998)

# function allocate_soil_properties(nodes, soilprops)
#     (; ρ_dry, θ_sat, λ_mineral, cp_mineral, ρ_mineral) = soilprops
#     NON = length(nodes)

#     λ_b = fill(λ_mineral[1], NON)
#     cp_b = fill(cp_mineral[1], NON)
#     ρ_b = fill(ρ_mineral[1], NON)

#     return (; λ_b, cp_b, ρ_b)
# end


# soil_properties(; T_soil, θ_soil, nodes::AbstractVector{<:Real}, soilprops, kw...) = 
#     soil_properties!(allocate_soil_properties(nodes, soilprops); T_soil, θ_soil, nodes, soilprops, kw...)
# function soil_properties!(buffers::NamedTuple;
#     T_soil::AbstractVector{<:Quantity},
#     θ_soil::AbstractVector{<:Real},
#     nodes::AbstractVector{<:Real},
#     soilprops::NamedTuple,
#     elevation::Quantity,
#     P_atmos = atmospheric_pressure(elevation),
#     runmoist::Bool,
#     runsnow::Bool,
# )
#     numtyps = findfirst(==(0.0), nodes) - 1
#     NON = length(nodes)
#     (; λ_b, cp_b, ρ_b) = buffers
#     (; ρ_dry, θ_sat, λ_mineral, cp_mineral, ρ_mineral) = soilprops

#     p_a0 = Unitful.atm
#     q = 4.0 # make a parameter with default, or q_0 * (T_soil / 303) ^ 2, q_0 = ~2 to 6 (power for recirculation function)
#     θ_0 = 0.162 # m3/m3 return-flow cutoff water content (~0.05 for coarse sand to 0.25 for heavy clay), p. 121 Campbell & Norman 1991 TODO make a parameter with default
#     g_a = 0.1 # make a parameter (~0.07 to 1.1), de Vries shape factor, 0.33 for organic soils, 0.1 for mineral
#     g_c = 1.0 - 2.0 * g_a


#     ϵ(λ_λ, λ_fluid) = 2.0 / (3.0 * (1.0 + g_a * (λ_λ / λ_fluid - 1.0))) + 1.0 / (3.0 * (1.0 + g_c * (λ_λ / λ_fluid - 1.0)))

#     ii, ij = runsnow ? (9, 8) : (1, 0)
#     for i in ii:NON
#         j = searchsortedlast(nodes, i - ij)
#         j = min(j, numtyps)

#         T_K = T_soil[i]
#         T_C = ustrip(u"°C", T_K)

#         m = runsnow ? θ_soil[i - 8] : θ_soil[i]
#         θj = θ_soil[j]

#         if runmoist
#             cp_b[i] = ρ_dry[j] / ρ_mineral[j] * cp_mineral[j] + m * cp_water
#             ρ_b[i] = m * ρ_water + ρ_dry[j]
#         else
#             cp_b[i] = ρ_dry[j] / ρ_mineral[j] * cp_mineral[j] + θ_sat[j] * θj * cp_water
#             ρ_b[i] = θ_sat[j] * θj * ρ_water + ρ_dry[j]
#         end

#         λ_water = (0.554 + 2.24e-3 * T_C - 9.87e-6 * T_C^2)u"W/m/K" # eq. 8 Campbell et al. 1994
#         λ_dry_air = (0.024 + 7.73e-5 * T_C - 2.6e-8 * T_C^2)u"W/m/K" # eq. 9 Campbell et al. 1994

#         D_v = D_v0 * (p_a0 / P_atmos) * (T_K / 273.15u"K")^1.75 # p. 309 Campbell et al. 1994
#         ρ_hat = ρ_hat0 * (P_atmos / p_a0) * (273.15 / T_K) # p. 309 Campbell et al. 1994
#         λ_vapor = molar_enthalpy_of_vaporisation(T_K)

#         ################################################################
#         # This is some of the most expensive code in the package
#         # its inlined so most of the work in wet_air_properties is ignored
#         e_a = wet_air_properties(T_K; rh=99.0, P_atmos).P_vap
#         e_a1 = wet_air_properties(T_K - 1u"K"; rh=99.0, P_atmos).P_vap
#         e_a2 = wet_air_properties(T_K + 1u"K"; rh=99.0, P_atmos).P_vap
#         ################################################################

#         ∇x = (e_a2 - e_a1) / 2.0

#         ϕ_mineral = ρ_dry[j] / ρ_mineral[j]
#         θ = runmoist ? m : θj * θ_sat[j]
#         ϕ_gas = max(0.0, 1.0 - θ - ϕ_mineral)
#         f_water = 1.0 / (1.0 + (θ / θ_0)^(-q)) # eq. 3, Campbell et al. 1994

#         λ_gas = λ_dry_air + λ_vapor * ∇x * f_water * ρ_hat * D_v / (P_atmos - e_a)
#         λ_fluid = λ_gas + f_water * (λ_water - λ_gas) 

#         λ_b[i] = (θ * ϵ(λ_water, λ_fluid) * λ_water + ϕ_mineral * ϵ(λ_mineral[j], λ_fluid) * λ_mineral[j] + ϕ_gas * ϵ(λ_gas, λ_fluid) * λ_gas) /
#                  (θ * ϵ(λ_water, λ_fluid) + ϕ_mineral * ϵ(λ_mineral[j], λ_fluid) + ϕ_gas * ϵ(λ_gas, λ_fluid))
#     end

#     #HTOFN = 333500.0  # J/kg
#     # if runsnow
#     #     snowdens = 0.0
#     #     if densfun[1] > 0
#     #         if densfun[3] > 0
#     #             snowdens = (densfun[1] - densfun[2]) * (1.0 - exp(-densfun[3] * cursnow - densfun[4] * snowage)) + densfun[2]
#     #         else
#     #             snowdens = min(0.9167, densfun[1] * snowage + densfun[2])
#     #         end
#     #     end

#     #     if cursnow >= minsnow
#     #         e_a = wet_air_properties(T_L; rh=100, P_atmos = P_atmos).r_w
#     #         cpsnow = (2100.0 * snowdens + (1.005 + 1.82 * (RW / (1.0 + RW))) * 1000.0 * (1.0 - snowdens))
#     #         snowcond2 = (0.00395 + 0.00084 * snowdens * 1000.0 - 0.0000017756 * (snowdens * 1000.0)^2 +
#     #                      0.00000000380635 * (snowdens * 1000.0)^3) / 418.6 * 60.0

#     #         for i in 1:8
#     #             λ_b[i] = snowcond > 0.0 ? snowcond : snowcond2
#     #             cp_b[i] = (T_C > -0.45 && T_C <= 0.4) ? (cpsnow + HTOFN) / 4186.0 : cpsnow / 4186.0
#     #             ρ_b[i] = (T_C > 0.0 || (snode[i] < 1e-8 && snode[min(8, i + 1)] > 0.0)) ? snowdens : 0.0
#     #         end
#     #     else
#     #         for i in 1:8
#     #             λ_b[i] = λ_b[9]
#     #             cp_b[i] = cp_b[9]
#     #             ρ_b[i] = 0.0
#     #         end
#     #     end
#     # end

#     return (; λ_b, cp_b, ρ_b)
# end

function soil_props(; 
    T_soil::Q, 
    θ_soil::R, 
    soilprops::S, 
    elevation::E, 
    P_atmos = atmospheric_pressure(elevation),
) where {Q<:Quantity, R<:Real, S<:NamedTuple, E<:Quantity}

    (; ρ_dry, θ_sat, λ_mineral, cp_mineral, ρ_mineral) = soilprops

    p_a0 = 101325u"Pa"
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
            θ_sat     = soilprops.θ_sat[i],
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

