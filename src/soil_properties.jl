# Julia translation of FORTRAN code that calculates variable thermal conductivity and specific heat
# for soil and snow layers, based on Campbell et al. (1994) and Campbell & Norman (1998)

function soil_properties(
    T_soil::AbstractVector{<:Quantity},
    θ_soil::AbstractVector{<:Real},
    nodes::AbstractVector{<:Real},
    soilprops::NamedTuple,
    elevation::Quantity,
    runmoist::Bool,
    runsnow::Bool,
)
    NON = length(nodes)
    numtyps = findfirst(nodes .== 0.0) - 1
    (; ρ_dry, θ_sat, λ_m, cp_m, ρ_m) = soilprops
    λ_b = fill(λ_m[1], NON)
    cp_b = fill(cp_m[1], NON)
    ρ_b = fill(ρ_m[1], NON)

    cp_water = 4184.0u"J/kg/K"
    ρ_water = 1000.0u"kg/m^3"
    ρ_hat0 = 44.65u"mol/m^3"
    D_v0 = 2.12e-5u"m^2/s"
    q = 4.0
    θ_0 = 0.162
    g_a = 0.1
    g_c = 1.0 - 2.0 * g_a
    p_a0 = 101325.0u"Pa"
    p_a = atmospheric_pressure(elevation)

    ϵ(λ_λ, λ_f) = 2.0 / (3.0 * (1.0 + g_a * (λ_λ / λ_f - 1.0))) + 1.0 / (3.0 * (1.0 + g_c * (λ_λ / λ_f - 1.0)))

    ii, ij = runsnow ? (9, 8) : (1, 0)
    for i in ii:NON
        j = searchsortedlast(nodes, i - ij)
        j = min(j, numtyps)

        T_K = T_soil[i]
        T_C = ustrip(u"°C", T_K)

        m = runsnow ? θ_soil[i - 8] : θ_soil[i]
        θj = θ_soil[j]

        if runmoist
            cp_b[i] = ρ_dry[j] / ρ_m[j] * cp_m[j] + m * cp_water
            ρ_b[i] = m * ρ_water + ρ_dry[j]
        else
            cp_b[i] = ρ_dry[j] / ρ_m[j] * cp_m[j] + θ_sat[j] * θj * cp_water
            ρ_b[i] = θ_sat[j] * θj * ρ_water + ρ_dry[j]
        end

        λ_a = (0.024 + 7.73e-5 * T_C - 2.6e-8 * T_C^2)u"W/m/K"
        λ_w = (0.554 + 2.24e-3 * T_C - 9.87e-6 * T_C^2)u"W/m/K"

        D_v = D_v0 * (p_a0 / p_a) * (T_K / 273.15u"K")^1.75
        ρ_hat = ρ_hat0 * (p_a / p_a0) * (273.15 / T_K)
        λ_vap = (45144.0 - 48.0 * T_C)u"J/mol"

        e_a = wet_air_properties(T_K; rh=99.0, P_atmos=p_a).P_vap
        e_a1 = wet_air_properties(T_K - 1u"K"; rh = 99.0, P_atmos = p_a).P_vap
        e_a2 = wet_air_properties(T_K + 1u"K"; rh = 99.0, P_atmos = p_a).P_vap
        ∇x = (e_a2 - e_a1) / 2.0

        ϕ_m = ρ_dry[j] / ρ_m[j]
        θ = runmoist ? m : θj * θ_sat[j]
        ϕ_g = max(0.0, 1.0 - θ - ϕ_m)
        f_w = 1.0 / (1.0 + (θ / θ_0)^(-q))

        λ_g = λ_a + λ_vap * ∇x * f_w * ρ_hat * D_v / (p_a - e_a)
        λ_f = λ_g + f_w * (λ_w - λ_g)

        λ_b[i] = (θ * ϵ(λ_w, λ_f) * λ_w + ϕ_m * ϵ(λ_m[j], λ_f) * λ_m[j] + ϕ_g * ϵ(λ_g, λ_f) * λ_g) /
                 (θ * ϵ(λ_w, λ_f) + ϕ_m * ϵ(λ_m[j], λ_f) + ϕ_g * ϵ(λ_g, λ_f))
    end

    #HTOFN = 333500.0  # J/kg
    # if runsnow
    #     snowdens = 0.0
    #     if densfun[1] > 0
    #         if densfun[3] > 0
    #             snowdens = (densfun[1] - densfun[2]) * (1.0 - exp(-densfun[3] * cursnow - densfun[4] * snowage)) + densfun[2]
    #         else
    #             snowdens = min(0.9167, densfun[1] * snowage + densfun[2])
    #         end
    #     end

    #     if cursnow >= minsnow
    #         e_a = wet_air_properties(T_L; rh=100, P_atmos = p_a).r_w
    #         cpsnow = (2100.0 * snowdens + (1.005 + 1.82 * (RW / (1.0 + RW))) * 1000.0 * (1.0 - snowdens))
    #         snowcond2 = (0.00395 + 0.00084 * snowdens * 1000.0 - 0.0000017756 * (snowdens * 1000.0)^2 +
    #                      0.00000000380635 * (snowdens * 1000.0)^3) / 418.6 * 60.0

    #         for i in 1:8
    #             λ_b[i] = snowcond > 0.0 ? snowcond : snowcond2
    #             cp_b[i] = (T_C > -0.45 && T_C <= 0.4) ? (cpsnow + HTOFN) / 4186.0 : cpsnow / 4186.0
    #             ρ_b[i] = (T_C > 0.0 || (snode[i] < 1e-8 && snode[min(8, i + 1)] > 0.0)) ? snowdens : 0.0
    #         end
    #     else
    #         for i in 1:8
    #             λ_b[i] = λ_b[9]
    #             cp_b[i] = cp_b[9]
    #             ρ_b[i] = 0.0
    #         end
    #     end
    # end

    return λ_b, cp_b, ρ_b
end
