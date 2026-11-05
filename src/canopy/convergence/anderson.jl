"""
    AndersonCanopyConvergence(; leaf_balance_tolerance=2e-5u"W/m^2", air_temperature_tolerance=1e-6u"K",
        vapor_density_tolerance=1e-6u"g/m^3", inverse_obukhov_length_tolerance=1e-8u"m^-1",
        history=5, damping=1.0, maxiters=100, on_failure=ThrowError(), backend=AutoEnzyme(; mode=Enzyme.Forward))

Solves the canopy for one hour by improving a guess until nothing changes.

Each pass updates the air and vapor from how they move through the canopy, and each leaf
from a Newton step on its energy balance. Anderson acceleration then uses the last `history`
passes to take a larger step towards the answer. `damping` is the share of each new pass
that is used, and `backend` differentiates the leaf energy balance.

The hour is solved when every value is within its tolerance. A layer without leaves uses
`air_temperature_tolerance` for its leaf.

A pass can overshoot to temperatures outside [`CANOPY_TEMPERATURE_BOUNDS`](@ref). The solver
then moves halfway back to its last good guess and forgets its history. `on_failure` decides
what happens if the hour starts out of bounds or hasn't settled after `maxiters` passes.

# References
- Walker, H. F., & Ni, P. (2011). Anderson acceleration for fixed-point iterations.
  *SIAM Journal on Numerical Analysis*, 49(4), 1715-1735.
"""
@kwdef struct AndersonCanopyConvergence{LB,AT,VD,IL,D,F<:AbstractCanopyFailure,B} <: AbstractCanopyConvergenceModel
    leaf_balance_tolerance::LB = 2e-5u"W/m^2"
    air_temperature_tolerance::AT = 1e-6u"K"
    vapor_density_tolerance::VD = 1e-6u"g/m^3"
    inverse_obukhov_length_tolerance::IL = 1e-8u"m^-1"
    history::Int = 5
    damping::D = 1.0
    maxiters::Int = 100
    on_failure::F = ThrowError()
    backend::B = AutoEnzyme(; mode=Enzyme.Forward)
end

function allocate_canopy_convergence(convergence_model::AndersonCanopyConvergence, buffers, model, boundary_layer_model;
    site, vapour_pressure_equation, soil_hydraulic_model,
)
    state_length = canopy_state_length(buffers)
    layer_count = length(buffers.leaf.leaf_temperature)
    (; history, backend) = convergence_model
    rows = canopy_state_rows(layer_count)
    tolerances = canopy_tolerances(convergence_model)
    # Puts each row's step in units of its tolerance. A leaf's step is a temperature.
    row_scale = map(1:state_length) do row
        row in rows.leaf_temperature ? 1 / tolerances.air_temperature : 1 / canopy_row_tolerance(tolerances, buffers, row)
    end
    forcing = Ref(initial_residual_forcing(model, soil_hydraulic_model))
    leaf_preparation = prepare_derivative(layer_leaf_balance_density, backend, 288.0,
        Constant(buffers), Constant(model), Constant(forcing[]), Constant(1))
    return (; forcing, tolerances, leaf_buffers=buffers, leaf_preparation, backend, row_scale,
        state=zeros(state_length), best_state=zeros(state_length), last_good_state=zeros(state_length),
        residual=zeros(state_length), next_state=zeros(state_length), step=zeros(state_length),
        previous_next_state=zeros(state_length), previous_step=zeros(state_length),
        next_state_changes=zeros(state_length, history), step_changes=zeros(state_length, history),
        history_matrix=zeros(history, history), history_weights=zeros(history))
end

# One layer's leaf energy balance per unit leaf area (W/m²) at `leaf_temperature` (K)
function layer_leaf_balance_density(leaf_temperature, buffers, model, forcing, layer)
    (; balance, leaf_area) = layer_leaf_balance(buffers, model, forcing, layer, leaf_temperature * u"K")
    return ustrip(u"W/m^2", balance.net / leaf_area)
end

# One pass: writes the state this pass would move to into `next_state`
function canopy_fixed_point!(next_state, residual, state, buffers, model, boundary_layer_model, forcing;
    site, vapour_pressure_equation,
)
    (; leaf_buffers, leaf_preparation, backend) = buffers.convergence
    evaluate_canopy!(residual, state, buffers, model, boundary_layer_model, forcing; site, vapour_pressure_equation)
    rows = canopy_state_rows(length(buffers.leaf.leaf_temperature))
    for (layer, row) in enumerate(rows.leaf_temperature)
        if iszero(layer_leaf_areas(buffers, layer).exchange_fraction)
            # A layer without leaves takes the air temperature
            next_state[row] = state[row] - residual[row]
            continue
        end
        slope = derivative(layer_leaf_balance_density, leaf_preparation, backend, state[row],
            Constant(leaf_buffers), Constant(model), Constant(forcing), Constant(layer))
        next_state[row] = state[row] - residual[row] / slope
    end
    # The residuals of air, vapor and stability are already the difference to move by
    for row in first(rows.air_temperature):rows.inverse_obukhov_length
        next_state[row] = state[row] - residual[row]
    end
    return next_state
end

# Finds the weight of each stored pass that best cancels the current step, by least squares.
# Returns false if the stored passes are too alike to separate.
function solve_history_weights!(history_weights, history_matrix, step_changes, step, row_scale, stored_count)
    weights = view(history_weights, 1:stored_count)
    matrix = view(history_matrix, 1:stored_count, 1:stored_count)
    for first_column in 1:stored_count
        weights[first_column] = 0.0
        for row in eachindex(step)
            weights[first_column] += step_changes[row, first_column] * step[row] * row_scale[row]
        end
        for second_column in 1:stored_count
            matrix[first_column, second_column] = 0.0
            for row in eachindex(step)
                matrix[first_column, second_column] += step_changes[row, first_column] * step_changes[row, second_column]
            end
        end
        matrix[first_column, first_column] *= 1.0 + 1.0e-10
    end
    factorization = cholesky!(Symmetric(matrix); check=false)
    issuccess(factorization) || return false
    ldiv!(factorization, weights)
    return true
end

function converge_canopy!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs,
    initial_wind, convergence_model::AndersonCanopyConvergence; site, vapour_pressure_equation, soil_hydraulic_model,
)
    (; forcing, tolerances, state, best_state, last_good_state, residual, next_state, step, previous_next_state,
       previous_step, next_state_changes, step_changes, history_matrix, history_weights, row_scale) = buffers.convergence
    (; damping, history, maxiters, on_failure) = convergence_model
    keywords = (; site, vapour_pressure_equation, soil_hydraulic_model)
    rows = canopy_state_rows(length(buffers.leaf.leaf_temperature))
    canopy_hour_setup!(state, forcing, buffers, model, boundary_layer_model, inputs, initial_wind; keywords...)

    lowest_temperature, highest_temperature = ustrip.(u"K", CANOPY_TEMPERATURE_BOUNDS)
    stored_count = 0
    oldest_column = 1
    smallest_residual = Inf
    copyto!(best_state, state)
    for pass in 1:maxiters
        canopy_fixed_point!(next_state, residual, state, buffers, model, boundary_layer_model, forcing[];
            site, vapour_pressure_equation)
        if !all(isfinite, next_state)
            pass == 1 && return canopy_failure!(on_failure, convergence_model, SciMLBase.ReturnCode.Unstable, best_state,
                buffers, model, boundary_layer_model, forcing[], inputs; keywords...)
            # Overshot: go halfway back to the last good state, without the history that led here
            for row in eachindex(state)
                state[row] = (state[row] + last_good_state[row]) / 2
            end
            stored_count = 0
            oldest_column = 1
            continue
        end
        copyto!(last_good_state, state)
        scaled_residual = canopy_scaled_residual(residual, tolerances, buffers)
        scaled_residual <= 1 && return canopy_solution!(residual, state, buffers, model, boundary_layer_model, forcing[],
            inputs, pass; keywords...)
        if scaled_residual < smallest_residual
            smallest_residual = scaled_residual
            copyto!(best_state, state)
        end
        for row in eachindex(step)
            step[row] = next_state[row] - state[row]
        end

        if pass > 1
            # Fill the free columns, then overwrite the oldest
            if stored_count < history
                stored_count += 1
                column = stored_count
            else
                column = oldest_column
                oldest_column = mod1(oldest_column + 1, history)
            end
            for row in eachindex(step)
                step_changes[row, column] = (step[row] - previous_step[row]) * row_scale[row]
                next_state_changes[row, column] = next_state[row] - previous_next_state[row]
            end
        end
        copyto!(previous_step, step)
        copyto!(previous_next_state, next_state)

        if stored_count > 0 && !solve_history_weights!(history_weights, history_matrix, step_changes, step, row_scale, stored_count)
            stored_count = 0
            oldest_column = 1
        end

        # Take the damped step, then subtract the weighted stored passes
        for row in eachindex(state)
            updated = state[row] + damping * step[row]
            for column in 1:stored_count
                updated -= (next_state_changes[row, column] -
                    (1.0 - damping) * step_changes[row, column] / row_scale[row]) * history_weights[column]
            end
            state[row] = updated
        end
        for row in first(rows.leaf_temperature):last(rows.air_temperature)
            state[row] = clamp(state[row], lowest_temperature, highest_temperature)
        end
        for row in rows.vapor_density
            state[row] = max(state[row], 0.0)
        end
    end
    return canopy_failure!(on_failure, convergence_model, SciMLBase.ReturnCode.MaxIters, best_state,
        buffers, model, boundary_layer_model, forcing[], inputs; keywords...)
end
