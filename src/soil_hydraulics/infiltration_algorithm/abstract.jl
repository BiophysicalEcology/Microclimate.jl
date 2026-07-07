abstract type AbstractInfiltrationAlgorithm end

"""
    initial_state_from_water_content!(algorithm, buffers, campbell_b, saturated_conductivity, num_layers)

Fill the algorithm's state variable and `hydraulic_conductivity` from the
initial `water_content_new`, for layers `2:num_layers` plus the saturated
lower boundary.
"""
function initial_state_from_water_content! end

"""
    algorithm_hydraulic_conductivity!(algorithm, buffers, i, saturated_conductivity)

Set `hydraulic_conductivity[i]` from the algorithm's current state variable.
"""
function algorithm_hydraulic_conductivity! end

"""
    state_to_water_potential!(algorithm, buffers, campbell_b, num_layers)

Bridge back to real matric potential (`water_potential`, `J/kg`), which every
downstream consumer (vapor flux, soil humidity, root/plant water uptake,
diagnostics) reads regardless of algorithm. No-op for
[`MatricPotentialAlgorithm`](@ref).
"""
function state_to_water_potential! end

"""
    assemble_tridiagonal_system!(algorithm, buffers, i, dt, campbell_b)

Fill the tridiagonal entries and `mass_balance_residual[i]` for layer `i`,
including the shared vapor-flux and root-water-uptake source terms.
"""
function assemble_tridiagonal_system! end

"""
    apply_thomas_elimination!(algorithm, buffers, num_layers)

Run the Thomas forward-elimination sweep over whichever tridiagonal arrays
the algorithm assembled.
"""
function apply_thomas_elimination! end

"""
    apply_lower_boundary_clamp!(algorithm, buffers, num_layers)

Back-substitute the state-variable change at the deepest interior layer and
apply the algorithm's stability clamp.
"""
function apply_lower_boundary_clamp! end

"""
    apply_interior_clamp!(algorithm, buffers, i)

Back-substitute the state-variable change at layer `i` and apply the
algorithm's stability clamp.
"""
function apply_interior_clamp! end

"""
    recover_water_content!(algorithm, buffers, campbell_b, num_layers)

Recompute `water_content_new` from the just-updated state variable.
"""
function recover_water_content! end

# Generic Thomas-algorithm forward elimination, shared by every AbstractInfiltrationAlgorithm
# via apply_thomas_elimination! — doesn't depend on what the arrays physically represent.
function thomas_forward_elimination!(sub_diagonal, diagonal, super_diagonal,
        normalized_super_diagonal, mass_balance_residual, normalized_residual, num_layers)
    for i in 2:num_layers-1
        normalized_super_diagonal[i] = super_diagonal[i] / diagonal[i]
        if abs(normalized_super_diagonal[i]) < 1e-8
            normalized_super_diagonal[i] = zero(normalized_super_diagonal[i])
        end
        normalized_residual[i] = mass_balance_residual[i] / diagonal[i]
        diagonal[i+1] -= sub_diagonal[i+1] * normalized_super_diagonal[i]
        mass_balance_residual[i+1] -= sub_diagonal[i+1] * normalized_residual[i]
    end
    return nothing
end
