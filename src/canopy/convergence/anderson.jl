"""
    AndersonCanopyConvergence(; leaf_balance_tolerance=2e-5u"W/m^2", air_temperature_tolerance=1e-6u"K",
        vapor_density_tolerance=1e-6u"g/m^3", inverse_obukhov_length_tolerance=1e-8u"m^-1",
        history=5, damping=1.0, maxiters=100, backend=AutoEnzyme(; mode=Enzyme.Forward))

Jacobian-free solve of [`evaluate_canopy!`](@ref) by Anderson-accelerated fixed-point passes:
air and vapor from their transport, each leaf from a Newton step on its own balance,
with the balance's derivative from `backend` through DifferentiationInterface.

Converges when every residual row is within the tolerance for its kind, as for
[`NonlinearCanopyConvergence`](@ref), and errors if it doesn't within `maxiters` passes.

# References
- Walker, H. F., & Ni, P. (2011). Anderson acceleration for fixed-point iterations.
  *SIAM Journal on Numerical Analysis*, 49(4), 1715-1735.
"""
@kwdef struct AndersonCanopyConvergence{LB,AT,VD,IL,D,B} <: AbstractCanopyConvergenceModel
    leaf_balance_tolerance::LB = 2e-5u"W/m^2"
    air_temperature_tolerance::AT = 1e-6u"K"
    vapor_density_tolerance::VD = 1e-6u"g/m^3"
    inverse_obukhov_length_tolerance::IL = 1e-8u"m^-1"
    history::Int = 5
    damping::D = 1.0
    maxiters::Int = 100
    backend::B = AutoEnzyme(; mode=Enzyme.Forward)
end

function allocate_canopy_convergence(convergence_model::AndersonCanopyConvergence, buffers, model, boundary_layer_model;
    site, vapour_pressure_equation, soil_hydraulic_model,
)
    n_state = canopy_state_length(buffers)
    n_layers = length(buffers.leaf.leaf_temperature)
    m = convergence_model.history
    rows = canopy_state_rows(n_layers)
    tolerances = canopy_tolerances(convergence_model)
    # Each pass's step over its tolerance; a leaf step is a temperature
    scale = map(1:n_state) do row
        row in rows.leaf_temperature ? 1 / tolerances.air_temperature : 1 / canopy_row_tolerance(tolerances, buffers, row)
    end
    forcing = Ref(initial_residual_forcing(model, soil_hydraulic_model))
    (; backend) = convergence_model
    leaf_preparation = prepare_derivative(layer_leaf_balance_density, backend, 288.0,
        Constant(buffers), Constant(model), Constant(forcing[]), Constant(1))
    return (; forcing, tolerances, leaf_buffers=buffers, leaf_preparation, backend,
        state=zeros(n_state), residual=zeros(n_state), mapped=zeros(n_state), step=zeros(n_state),
        previous_mapped=zeros(n_state), previous_step=zeros(n_state),
        mapped_differences=zeros(n_state, m), step_differences=zeros(n_state, m),
        normal_matrix=zeros(m, m), normal_rhs=zeros(m), scale)
end

# Layer `layer`'s leaf energy balance per unit leaf area (W/m²) at `leaf_temperature` (K)
function layer_leaf_balance_density(leaf_temperature, buffers, model, forcing, layer)
    (; balance, leaf_area) = layer_leaf_balance(buffers, model, forcing, layer, leaf_temperature * u"K")
    return ustrip(u"W/m^2", balance.net / leaf_area)
end

# One fixed-point pass: `mapped` from `state`
function canopy_fixed_point!(mapped, residual, state, buffers, model, boundary_layer_model, forcing;
    site, vapour_pressure_equation,
)
    (; leaf_buffers, leaf_preparation, backend) = buffers.convergence
    evaluate_canopy!(residual, state, buffers, model, boundary_layer_model, forcing; site, vapour_pressure_equation)
    rows = canopy_state_rows(length(buffers.leaf.leaf_temperature))
    for (i, j) in enumerate(rows.leaf_temperature)
        if iszero(layer_leaf_areas(buffers, i).exchange_fraction)
            # The residual is T_leaf - T_air: the leaf takes the air temperature
            mapped[j] = state[j] - residual[j]
            continue
        end
        slope = derivative(layer_leaf_balance_density, leaf_preparation, backend, state[j],
            Constant(leaf_buffers), Constant(model), Constant(forcing), Constant(i))
        mapped[j] = state[j] - residual[j] / slope
    end
    # Air, vapor and stability rows already read state minus transport
    for j in first(rows.air_temperature):rows.inverse_obukhov_length
        mapped[j] = state[j] - residual[j]
    end
    return mapped
end

function converge_canopy!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs,
    initial_wind, cm::AndersonCanopyConvergence; site, vapour_pressure_equation, soil_hydraulic_model,
)
    (; forcing, tolerances, state, residual, mapped, step, previous_mapped, previous_step, mapped_differences,
       step_differences, normal_matrix, normal_rhs, scale) = buffers.convergence
    rows = canopy_state_rows(length(buffers.leaf.leaf_temperature))
    canopy_hour_setup!(state, forcing, buffers, model, inputs, initial_wind;
        site, vapour_pressure_equation, soil_hydraulic_model)

    β = cm.damping
    m = cm.history
    lo, hi = ustrip.(u"K", CANOPY_TEMPERATURE_BOUNDS)
    stored = 0
    oldest = 1
    for k in 1:cm.maxiters
        canopy_fixed_point!(mapped, residual, state, buffers, model, boundary_layer_model, forcing[];
            site, vapour_pressure_equation)
        all(isfinite, mapped) || canopy_convergence_error(cm, SciMLBase.ReturnCode.Unstable, inputs)
        converged = all(eachindex(residual)) do row
            abs(residual[row]) <= canopy_row_tolerance(tolerances, buffers, row)
        end
        converged && return canopy_solution!(residual, state, buffers, model, boundary_layer_model, forcing[],
            inputs, k; site, vapour_pressure_equation, soil_hydraulic_model)
        for j in eachindex(step)
            step[j] = mapped[j] - state[j]
        end

        if k > 1
            # Fill free columns, then overwrite the oldest
            if stored < m
                stored += 1
                column = stored
            else
                column = oldest
                oldest = mod1(oldest + 1, m)
            end
            for j in eachindex(step)
                step_differences[j, column] = (step[j] - previous_step[j]) * scale[j]
                mapped_differences[j, column] = mapped[j] - previous_mapped[j]
            end
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
                oldest = 1
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
        for i in first(rows.leaf_temperature):last(rows.air_temperature)
            state[i] = clamp(state[i], lo, hi)
        end
        for i in rows.vapor_density
            state[i] = max(state[i], 0.0)
        end
    end
    canopy_convergence_error(cm, SciMLBase.ReturnCode.MaxIters, inputs)
end
