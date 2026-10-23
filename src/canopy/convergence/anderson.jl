"""
    AndersonCanopyConvergence(; tolerance=0.1u"K", history=5, damping=1.0, maxiters=100,
        fallback=NonlinearCanopyConvergence(), fallback_threshold=10.0)

Jacobian-free solve of [`evaluate_canopy!`](@ref) by Anderson-accelerated fixed-point passes:
air and vapor from their transport, each leaf from a Newton step on its own balance.
Stops when no temperature changes by more than `tolerance`; else keeps the steadiest pass,
or solves the hour with `fallback` if that pass is more than `fallback_threshold`
tolerances from settled.

# References
- Walker, H. F., & Ni, P. (2011). Anderson acceleration for fixed-point iterations.
  *SIAM Journal on Numerical Analysis*, 49(4), 1715-1735.
"""
@kwdef struct AndersonCanopyConvergence{T,D,F,FT} <: AbstractCanopyConvergenceModel
    tolerance::T = 0.1u"K"
    history::Int = 5
    damping::D = 1.0
    maxiters::Int = 100
    fallback::F = NonlinearCanopyConvergence()
    fallback_threshold::FT = 10.0
end

function allocate_canopy_convergence(convergence_model::AndersonCanopyConvergence, buffers, model, boundary_layer_model;
    site, vapour_pressure_equation, soil_hydraulic_model,
)
    n_state = canopy_state_length(buffers)
    n_layers = length(buffers.leaf.leaf_temperature)
    m = convergence_model.history
    # Row weights: K, g/m³, and 1/m × 100 m
    scale = vcat(ones(3n_layers), 100.0)
    return (; forcing=Ref(initial_residual_forcing(model, soil_hydraulic_model)),
        state=zeros(n_state), residual=zeros(n_state), mapped=zeros(n_state), step=zeros(n_state),
        best=zeros(n_state), previous_mapped=zeros(n_state), previous_step=zeros(n_state),
        mapped_differences=zeros(n_state, m), step_differences=zeros(n_state, m),
        normal_matrix=zeros(m, m), normal_rhs=zeros(m), scale,
        fallback=allocate_canopy_convergence(convergence_model.fallback, buffers, model, boundary_layer_model;
            site, vapour_pressure_equation, soil_hydraulic_model))
end

# One fixed-point pass: `mapped` from `state`
function canopy_fixed_point!(mapped, residual, state, buffers, model, boundary_layer_model, forcing;
    site, vapour_pressure_equation,
)
    evaluate_canopy!(residual, state, buffers, model, boundary_layer_model, forcing; site, vapour_pressure_equation)
    (; leaf_temperature, net_balance) = buffers.leaf
    n_layers = length(leaf_temperature)
    δ = 0.1u"K"
    for i in 1:n_layers
        slope = 1.0
        if !iszero(1.0 - buffers.longwave.layer_transmission[i])
            (; balance, leaf_area) = layer_leaf_balance(buffers, model, forcing, i, leaf_temperature[i] + δ;
                vapour_pressure_equation)
            # Net balance falls as the leaf warms
            slope = min(ustrip(u"W/m^2/K", (balance.net - net_balance[i]) / (δ * leaf_area)), -1.0)
            slope = -slope
        end
        mapped[i] = state[i] + residual[i] / slope
    end
    for j in (n_layers + 1):length(state)
        mapped[j] = state[j] - residual[j]
    end
    return mapped
end

function converge_canopy!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs,
    initial_wind, cm::AndersonCanopyConvergence; site, vapour_pressure_equation, soil_hydraulic_model,
)
    (; forcing, state, residual, mapped, step, best, previous_mapped, previous_step, mapped_differences,
       step_differences, normal_matrix, normal_rhs, scale) = buffers.convergence
    n_layers = length(buffers.leaf.leaf_temperature)
    canopy_hour_setup!(state, forcing, buffers, model, inputs, initial_wind;
        site, vapour_pressure_equation, soil_hydraulic_model)

    tolerance = ustrip(u"K", cm.tolerance)
    β = cm.damping
    m = cm.history
    lo, hi = ustrip.(u"K", CANOPY_TEMPERATURE_BOUNDS)
    stored = 0
    best_change = Inf
    iterations = 0
    for k in 1:cm.maxiters
        iterations = k
        canopy_fixed_point!(mapped, residual, state, buffers, model, boundary_layer_model, forcing[];
            site, vapour_pressure_equation)
        change = 0.0
        for j in eachindex(step)
            step[j] = mapped[j] - state[j]
            change = max(change, abs(step[j]) * scale[j])
        end
        isfinite(change) || break
        if change < best_change
            best_change = change
            copyto!(best, mapped)
        end
        change <= tolerance && break

        if k > 1
            column = stored < m ? stored + 1 : mod1(k - 1, m)
            for j in eachindex(step)
                step_differences[j, column] = (step[j] - previous_step[j]) * scale[j]
                mapped_differences[j, column] = mapped[j] - previous_mapped[j]
            end
            stored = min(stored + 1, m)
        end
        copyto!(previous_step, step)
        copyto!(previous_mapped, mapped)

        # Least squares for the history weights, by normal equations
        γ = view(normal_rhs, 1:stored)
        if stored > 0
            A = view(normal_matrix, 1:stored, 1:stored)
            for a in 1:stored
                γ[a] = 0.0
                for j in eachindex(step)
                    γ[a] += step_differences[j, a] * step[j] * scale[j]
                end
                for b in 1:stored
                    A[a, b] = 0.0
                    for j in eachindex(step)
                        A[a, b] += step_differences[j, a] * step_differences[j, b]
                    end
                end
                A[a, a] *= 1.0 + 1.0e-10
            end
            factorization = cholesky!(Symmetric(A); check=false)
            if issuccess(factorization)
                ldiv!(factorization, γ)
            else
                stored = 0
            end
        end

        # x ← x + βf − Σ (ΔG - (1-β)ΔF)γ
        for j in eachindex(state)
            next = state[j] + β * step[j]
            for a in 1:stored
                next -= (mapped_differences[j, a] - (1.0 - β) * step_differences[j, a] / scale[j]) * γ[a]
            end
            state[j] = next
        end
        for i in 1:2n_layers
            state[i] = clamp(state[i], lo, hi)
        end
        for i in (2n_layers + 1):3n_layers
            state[i] = max(state[i], 0.0)
        end
    end
    if best_change > cm.fallback_threshold * tolerance
        fallback_buffers = merge(buffers, (; convergence=buffers.convergence.fallback))
        result = converge_canopy!(fallback_buffers, model, boundary_layer_model, inputs, initial_wind, cm.fallback;
            site, vapour_pressure_equation, soil_hydraulic_model)
        return Base.setindex(result, iterations + result[2], 2)
    end
    copyto!(state, best)
    return canopy_solution!(residual, state, buffers, model, boundary_layer_model, forcing[], inputs, iterations;
        site, vapour_pressure_equation, soil_hydraulic_model)
end
