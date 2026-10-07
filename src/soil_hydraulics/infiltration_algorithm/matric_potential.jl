"""
    MatricPotentialAlgorithm()

Campbell (1985) Program 8.1/11.1: solves the infiltration core using matric
potential ψ (J/kg) as the dependent variable. Well-behaved for dry soil; needs
an averaging clamp near saturation to stay stable at high water content (see
[`MatricFluxPotentialAlgorithm`](@ref) for the alternative that avoids this).
"""
struct MatricPotentialAlgorithm <: AbstractInfiltrationAlgorithm end

function initial_state_from_water_content!(::MatricPotentialAlgorithm, buffers, campbell_b, saturated_conductivity, num_layers)
    (; water_potential, water_content_new, air_entry_potential, saturation_water_content,
       campbell_exponent, hydraulic_conductivity) = buffers
    for i in 2:num_layers
        water_potential[i] = air_entry_potential[i] * (saturation_water_content[i] / water_content_new[i])^campbell_b[i] # matric water potential, EQ5.9
        hydraulic_conductivity[i] = saturated_conductivity[i] * (air_entry_potential[i] / water_potential[i])^campbell_exponent[i] # hydraulic conductivity, EQ6.14
    end
    water_potential[num_layers+1] = air_entry_potential[num_layers] # saturated lower boundary
    return nothing
end

function algorithm_hydraulic_conductivity!(::MatricPotentialAlgorithm, buffers, i, saturated_conductivity)
    (; hydraulic_conductivity, water_potential, air_entry_potential, campbell_exponent) = buffers
    hydraulic_conductivity[i] = saturated_conductivity[i] * (air_entry_potential[i] / water_potential[i])^campbell_exponent[i]
    return nothing
end

state_to_water_potential!(::MatricPotentialAlgorithm, buffers, campbell_b, num_layers) = nothing

function assemble_tridiagonal_system!(::MatricPotentialAlgorithm, buffers, i, dt, campbell_b)
    (; water_potential, hydraulic_conductivity, water_content_new, water_content, depth,
       layer_water_mass, campbell_exponent, campbell_exponent_complement,
       vapor_flux, vapor_flux_derivative, rainfall_flux, root_water_uptake,
       sub_diagonal, diagonal, super_diagonal, mass_balance_residual) = buffers
    hydraulic_capacitance = -1.0 * layer_water_mass[i] * water_content_new[i] / (campbell_b[i] * water_potential[i] * dt)
    sub_diagonal[i] = -1.0 * hydraulic_conductivity[i-1] / (depth[i] - depth[i-1]) + Unitful.gn * campbell_exponent[i] * hydraulic_conductivity[i-1] / water_potential[i-1]
    super_diagonal[i] = -1.0 * hydraulic_conductivity[i+1] / (depth[i+1] - depth[i])
    diagonal[i] = hydraulic_conductivity[i] / (depth[i] - depth[i-1]) + hydraulic_conductivity[i] / (depth[i+1] - depth[i]) + hydraulic_capacitance - Unitful.gn * campbell_exponent[i] * hydraulic_conductivity[i] / water_potential[i] + vapor_flux_derivative[i-1] + vapor_flux_derivative[i]
    mass_balance_residual[i] = ((water_potential[i] * hydraulic_conductivity[i] - water_potential[i-1] * hydraulic_conductivity[i-1]) / (depth[i] - depth[i-1]) - (water_potential[i+1] * hydraulic_conductivity[i+1] - water_potential[i] * hydraulic_conductivity[i]) / (depth[i+1] - depth[i])) / campbell_exponent_complement[i] + layer_water_mass[i] * (water_content_new[i] - water_content[i]) / dt - Unitful.gn * (hydraulic_conductivity[i-1] - hydraulic_conductivity[i]) + vapor_flux[i-1] - vapor_flux[i] + rainfall_flux[i-1] - rainfall_flux[i] + root_water_uptake[i]
    return nothing
end

function apply_thomas_elimination!(::MatricPotentialAlgorithm, buffers, num_layers)
    (; sub_diagonal, diagonal, super_diagonal, normalized_super_diagonal, mass_balance_residual, normalized_residual) = buffers
    thomas_forward_elimination!(sub_diagonal, diagonal, super_diagonal, normalized_super_diagonal, mass_balance_residual, normalized_residual, num_layers)
    return nothing
end

function apply_lower_boundary_clamp!(::MatricPotentialAlgorithm, buffers, num_layers)
    (; potential_change, mass_balance_residual, diagonal, water_potential, air_entry_potential) = buffers
    potential_change[num_layers] = mass_balance_residual[num_layers] / diagonal[num_layers]
    water_potential[num_layers] -= potential_change[num_layers]
    water_potential[num_layers] = min(water_potential[num_layers], air_entry_potential[num_layers])
    return nothing
end

function apply_interior_clamp!(::MatricPotentialAlgorithm, buffers, i)
    (; potential_change, normalized_residual, normalized_super_diagonal, water_potential, air_entry_potential) = buffers
    potential_change[i] = normalized_residual[i] - normalized_super_diagonal[i] * potential_change[i+1] # change in matric potential in an iteration step, J/kg
    water_potential[i] -= potential_change[i]
    if water_potential[i] > air_entry_potential[i]
        water_potential[i] = (water_potential[i] + potential_change[i] + air_entry_potential[i]) / 2.0
    end
    return nothing
end

function recover_water_content!(::MatricPotentialAlgorithm, buffers, campbell_b, num_layers)
    (; water_content_new, air_entry_potential, water_potential, campbell_b_inverse, saturation_water_content) = buffers
    for i in 2:num_layers
        water_content_new[i] = max(saturation_water_content[i] * (air_entry_potential[i] / water_potential[i])^campbell_b_inverse[i], 1e-7)
        water_potential[i] = air_entry_potential[i] * (saturation_water_content[i] / water_content_new[i])^campbell_b[i]
    end
    return nothing
end
