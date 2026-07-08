"""
    MatricFluxPotentialAlgorithm()

Campbell (1985) Program 8.2: solves the same infiltration physics as
[`MatricPotentialAlgorithm`](@ref) but using matric flux potential
Φ = ∫K dψ (`kg/(m*s)`, i.e. `Pa*s`) as the dependent variable, per Bittelli,
Campbell & Tomei, *Soil Physics with Python* (Ch. 8.4-8.5):

    Φ_e = Ks·ψ_e/(1-n),  n = 2+3/b       (8.46) — Φ at saturation
    Φ = Φ_e·(θ/θ_s)^(b+3)                (8.50)
    C = dθ/dΦ = θ/((b+3)·Φ)              (8.51)

Converges in fewer iterations and is more stable at high water content than
`MatricPotentialAlgorithm`, at the cost of being harder to resolve for dry
soil (θ(Φ) is far more convex than θ(ψ)).
"""
struct MatricFluxPotentialAlgorithm <: AbstractInfiltrationAlgorithm end

function initial_state_from_water_content!(::MatricFluxPotentialAlgorithm, buffers, campbell_b, saturated_conductivity, num_layers)
    (; matric_flux_potential, water_content_new, air_entry_flux_potential, saturation_water_content,
       campbell_flux_exponent, hydraulic_conductivity) = buffers
    for i in 2:num_layers
        matric_flux_potential[i] = air_entry_flux_potential[i] * (water_content_new[i] / saturation_water_content[i])^(campbell_b[i] + 3.0)
        hydraulic_conductivity[i] = saturated_conductivity[i] * (matric_flux_potential[i] / air_entry_flux_potential[i])^campbell_flux_exponent[i]
    end
    matric_flux_potential[num_layers+1] = air_entry_flux_potential[num_layers]
    return nothing
end

function algorithm_hydraulic_conductivity!(::MatricFluxPotentialAlgorithm, buffers, i, saturated_conductivity)
    (; hydraulic_conductivity, matric_flux_potential, air_entry_flux_potential, campbell_flux_exponent) = buffers
    hydraulic_conductivity[i] = saturated_conductivity[i] * (matric_flux_potential[i] / air_entry_flux_potential[i])^campbell_flux_exponent[i]
    return nothing
end

function state_to_water_potential!(::MatricFluxPotentialAlgorithm, buffers, campbell_b, num_layers)
    (; water_potential, water_content_new, air_entry_potential, saturation_water_content) = buffers
    for i in 2:num_layers
        water_potential[i] = air_entry_potential[i] * (saturation_water_content[i] / water_content_new[i])^campbell_b[i]
    end
    water_potential[num_layers+1] = air_entry_potential[num_layers]
    return nothing
end

# vapor_flux_derivative isn't added into flux_diagonal (unlike MatricPotentialAlgorithm's
# diagonal) since it's d(vapor flux)/dψ — converting to d/dΦ needs a dψ/dΦ=1/K factor with no
# reference to check it against. flux_diagonal is just the Newton Jacobian, not the residual,
# so mass_balance_residual (which does include vapor_flux) still converges correctly without it.
function assemble_tridiagonal_system!(::MatricFluxPotentialAlgorithm, buffers, i, dt, campbell_b)
    (; matric_flux_potential, hydraulic_conductivity, water_content_new, water_content, depth,
       layer_water_mass, campbell_flux_exponent,
       vapor_flux, root_water_uptake,
       flux_sub_diagonal, flux_diagonal, flux_super_diagonal, mass_balance_residual) = buffers
    b3 = campbell_flux_exponent[i]
    hydraulic_capacitance = layer_water_mass[i] * water_content_new[i] / ((campbell_b[i] + 3.0) * matric_flux_potential[i] * dt)
    flux_sub_diagonal[i] = -1.0 / (depth[i] - depth[i-1]) - Unitful.gn * b3 * hydraulic_conductivity[i-1] / matric_flux_potential[i-1]
    flux_super_diagonal[i] = -1.0 / (depth[i+1] - depth[i])
    flux_diagonal[i] = 1.0 / (depth[i] - depth[i-1]) + 1.0 / (depth[i+1] - depth[i]) + hydraulic_capacitance + Unitful.gn * b3 * hydraulic_conductivity[i] / matric_flux_potential[i]
    mass_balance_residual[i] = (matric_flux_potential[i] - matric_flux_potential[i-1]) / (depth[i] - depth[i-1]) - (matric_flux_potential[i+1] - matric_flux_potential[i]) / (depth[i+1] - depth[i]) + layer_water_mass[i] * (water_content_new[i] - water_content[i]) / dt - Unitful.gn * (hydraulic_conductivity[i-1] - hydraulic_conductivity[i]) + vapor_flux[i-1] - vapor_flux[i] + root_water_uptake[i]
    return nothing
end

function apply_thomas_elimination!(::MatricFluxPotentialAlgorithm, buffers, num_layers)
    (; flux_sub_diagonal, flux_diagonal, flux_super_diagonal, normalized_super_diagonal, mass_balance_residual, flux_normalized_residual) = buffers
    thomas_forward_elimination!(flux_sub_diagonal, flux_diagonal, flux_super_diagonal, normalized_super_diagonal, mass_balance_residual, flux_normalized_residual, num_layers)
    return nothing
end

function apply_lower_boundary_clamp!(::MatricFluxPotentialAlgorithm, buffers, num_layers)
    (; flux_potential_change, mass_balance_residual, flux_diagonal, matric_flux_potential, air_entry_flux_potential) = buffers
    flux_potential_change[num_layers] = mass_balance_residual[num_layers] / flux_diagonal[num_layers]
    matric_flux_potential[num_layers] -= flux_potential_change[num_layers]
    clamp_flux_potential!(matric_flux_potential, air_entry_flux_potential, num_layers)
    return nothing
end

function apply_interior_clamp!(::MatricFluxPotentialAlgorithm, buffers, i)
    (; flux_potential_change, flux_normalized_residual, normalized_super_diagonal, matric_flux_potential, air_entry_flux_potential) = buffers
    flux_potential_change[i] = flux_normalized_residual[i] - normalized_super_diagonal[i] * flux_potential_change[i+1]
    matric_flux_potential[i] -= flux_potential_change[i]
    clamp_flux_potential!(matric_flux_potential, air_entry_flux_potential, i)
    return nothing
end

# Φ can't exceed its saturated value Φ_e; the lower floor keeps θ(Φ)'s fractional power finite.
function clamp_flux_potential!(matric_flux_potential, air_entry_flux_potential, i)
    if matric_flux_potential[i] < zero(matric_flux_potential[i])
        matric_flux_potential[i] = 1e-6 * oneunit(matric_flux_potential[i])
    elseif matric_flux_potential[i] > air_entry_flux_potential[i]
        matric_flux_potential[i] = air_entry_flux_potential[i]
    end
    return nothing
end

function recover_water_content!(::MatricFluxPotentialAlgorithm, buffers, campbell_b, num_layers)
    (; water_content_new, matric_flux_potential, air_entry_flux_potential, saturation_water_content) = buffers
    for i in 2:num_layers
        water_content_new[i] = max(saturation_water_content[i] * (matric_flux_potential[i] / air_entry_flux_potential[i])^(1.0 / (campbell_b[i] + 3.0)), 1e-7)
    end
    return nothing
end
