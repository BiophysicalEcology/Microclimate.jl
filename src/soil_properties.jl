# Julia translation of FORTRAN code that calculates variable thermal conductivity and specific heat
# for soil and snow layers, based on Campbell et al. (1994) and Campbell & Norman (1998)

function soil_properties(T_soil, θ_soil, nodes, soilprops, elev, runmoist, runsnow)
    NON = length(nodes)
    #runmoist = false # to do
    #runsnow = false # to do
    numtyps = findfirst(nodes .== 0.0)-1
    θ_sat = soilprops[:, 2] # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    λ_m = soilprops[:, 3] # mineral thermal conductivity, W/m/K
    cp_m = soilprops[:, 4] # mineral specific heat capacity, J/kg/K
    ρ_m = u"kg/m^3".(soilprops[:, 5]) # mineral density, kg/m^3
    ρ_dry = u"kg/m^3".(soilprops[:, 1]) # dry density, kg/m^3
    λ_b = fill(λ_m[1], NON) # bulk thermal conductivity, W/m/K
    cp_b = fill(cp_m[1], NON) # bulk specific heat capacity, J/kg/K
    ρ_b = fill(ρ_m[1], NON) # bulk density, kg/m^3

    ii, ij = runsnow ? (9, 8) : (1, 0)
    j = 1
    for i in ii:NON # ii ensures that it starts at soil depth if snow on top 
        if i >= nodes[j] + ij
            j = min(j + 1, numtyps)
        end
        # don't make it volumetric, but rather mass-specific, so don't multiply by kg/m3 converters (constants from Campbell and Norman 1998, Table 8.2)
        if runmoist
            m = runsnow ? θ_soil[i-8] : θ_soil[i]
            cp_b[i] = ρ_dry[j] / ρ_m[j] * cp_m[j] + m * 4184.0u"J/kg/K"
        else
            cp_b[i] = ρ_dry[j] / ρ_m[j] * cp_m[j] + θ_sat[j] * θ_soil[j] * 4184.0u"J/kg/K"
        end
        # constants from Campbell and Norman 1998, Table 8.2
        if runmoist
            m = runsnow ? θ_soil[i-8] : θ_soil[i]
            ρ_b[i] = m * 1000.0u"kg/m^3" + ρ_dry[j]
        else
            ρ_b[i] = θ_sat[j] * θ_soil[j] * 1000.0u"kg/m^3" + ρ_dry[j]
        end

        p_a0 = 101325.0u"Pa" # standard sea level air pressure, Pa
        p_a = get_pressure(elev)

        T_K = T_soil[i]
        T_C = u"°C"(T_K)
        λ_a = (0.024 + 7.73e-5 * ustrip(T_C) - 2.6e-8 * ustrip(T_C)^2)u"W/m/K" # thermal conductivity of dry air (equation 9 in Campbell et al. 1994)
        λ_w = (0.554 + 2.24e-3 * ustrip(T_C) - 9.87e-6 * ustrip(T_C)^2)u"W/m/K" # thermal conductivity of water (equation 8 in Campbell et al. 1994)

        hr = 1.0 # relative humidity
        D_v0 = 2.12e-5u"m^2/s" # vapour diffusivity in air (m2/s), standard value at 0 deg C and sea level pressure (Campbell et al. 1994)
        ρ_hat0 = 44.65u"mol/m^3" # molar density of air (mol/m3), standard value at 0 deg C and sea level pressure (Campbell et al. 1994)

        D_v = D_v0 * (p_a0 / p_a) * (T_K / 273.15u"K")^1.75 # temperature/pressure-corrected vapour diffusivity in air (m2/s) (p. 309 in Campbell et al. 1994)
        ρ_hat = ρ_hat0 * (p_a / p_a0) * (273.15 / T_K) # temperature/pressure-corrected molar density of air (mol/m3) (p. 309 in Campbell et al. 1994)
        λ_vap = (45144.0 - 48.0 * ustrip(T_C))u"J/mol" # J/mol latent heat of vaporization (Cambell et al. 1994, p. 309)

        rh = 99.0
        e_a = wet_air(T_K; rh=rh, P_atmos = p_a).P_vap
        e_a1 = wet_air(T_K-1u"K"; rh=rh, P_atmos = p_a).P_vap
        e_a2 = wet_air(T_K+1u"K"; rh=rh, P_atmos = p_a).P_vap

        ∇x = (e_a2 - e_a1) / 2.0 # slope of the vapour pressure function centred at focal temperature

        # these could vary with soil texture but the relationship isn't strong
        # power for liquid recirculation, mean in Table 2 of Campell et al. 1994, excluding peat moss value
        # q_0=4
        # q=q_0*(T_L/303.0u"K")^2
        q = 4.0 # using a typical value for 'q' of 4 - program becomes unstable if this is temperature dependent
        θ_0 = 0.162 # mean in Table 2 of Campell et al. 1994, excluding peat moss value

        ϕ_m = ρ_dry[j] / ρ_m[j] # volume fraction of minerals
        θ = runmoist == 1 ? (runsnow == 1 ? θ_soil[i-8] : θ_soil[i]) : θ_soil[j] * θ_sat[j] # volume fraction of water

        ϕ_g = max(0.0, 1.0 - θ - ϕ_m) # volume fraction of gas

        f_w = 1.0 / (1.0 + (θ / θ_0)^(-q)) # eq 8.17 Campbell and Norman 1988, using temperature-specific q

        λ_g = λ_a + λ_vap * ∇x * hr * f_w * ρ_hat * D_v / (p_a - e_a) # eq 8.18 Campbell and Norman 1988
        λ_f = λ_g + f_w * (λ_w - λ_g) # eq 8.19 Campbell and Norman 1988

        g_a = 0.1 # 0.1 for mineral soils, 0.33 for organic, p 125, Campbell and Norman 1988
        g_c = 1.0 - 2.0 * g_a # p 125, Campbell and Norman

        # equation 8.20 in Campbell and Norman 1988
        ϵ_g = 2.0 / (3.0 * (1.0 + g_a * (λ_g / λ_f - 1.0))) + 1.0 / (3.0 * (1.0 + g_c * (λ_g / λ_f - 1.0))) 
        ϵ_w = 2.0 / (3.0 * (1.0 + g_a * (λ_w / λ_f - 1.0))) + 1.0 / (3.0 * (1.0 + g_c * (λ_w / λ_f - 1.0)))
        ϵ_m = 2.0 / (3.0 * (1.0 + g_a * (λ_m[j] / λ_f - 1.0))) + 1.0 / (3.0 * (1.0 + g_c * (λ_m[j] / λ_f - 1.0)))
        
        # equation 8.13 in Campbell and Norman 1988
        λ_b[i] = (θ * ϵ_w * λ_w + ϕ_m * ϵ_m * λ_m[j] + ϕ_g * ϵ_g * λ_g) /
                       (θ * ϵ_w + ϕ_m * ϵ_m + ϕ_g * ϵ_g)

        #convert to cal/min/g
        #λ_b[i] = λ_b[i] / 418.6 * 60.0
        #cp_b[i] /= 4186.0
        #ρ_b[i] /= 1000.0
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
    #         e_a = wet_air(T_L; rh=100, P_atmos = p_a).r_w
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