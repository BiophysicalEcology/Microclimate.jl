"""
    AndersonCanopyConvergence(; leaf_balance_tolerance=1.5u"W/m^2", air_temperature_tolerance=0.1u"K",
        vapor_density_tolerance=0.07u"g/m^3", ground_temperature_tolerance=0.1u"K", inverse_obukhov_length_tolerance=1e-2u"m^-1",
        history=5, damping=1.0, maxiters=100, on_failure=ThrowError(), backend=AutoEnzyme(; mode=Enzyme.Forward))

Solves the canopy for one hour by improving a guess until nothing changes.

Each pass updates the air and vapor from how they move through the canopy, each leaf from a
Newton step on its energy balance, and the ground from what the soil reaches by the end of the
hour. Anderson acceleration then uses the last `history`
passes to take a larger step towards the answer. `damping` is the share of each new pass
that is used, and `backend` differentiates the leaf energy balance.

The hour is solved when every value is within its tolerance. A layer without leaves uses
`air_temperature_tolerance` for its leaf.

A pass can overshoot past [`CANOPY_TEMPERATURE_BOUNDS`](@ref). The
solver then moves halfway back to its last good guess and forgets its history. `on_failure`
decides what happens if the hour starts out of bounds or hasn't settled after `maxiters` passes.

# References
- Walker, H. F., & Ni, P. (2011). Anderson acceleration for fixed-point iterations.
  *SIAM Journal on Numerical Analysis*, 49(4), 1715-1735.
"""
@kwdef struct AndersonCanopyConvergence{LB,AT,VD,GT,IL,D,F<:AbstractCanopyFailure,B} <: AbstractCanopyConvergenceModel
    leaf_balance_tolerance::LB = 1.5u"W/m^2"
    air_temperature_tolerance::AT = 0.1u"K"
    vapor_density_tolerance::VD = 0.07u"g/m^3"
    ground_temperature_tolerance::GT = 0.1u"K"
    inverse_obukhov_length_tolerance::IL = 1e-2u"m^-1"
    history::Int = 5
    damping::D = 1.0
    maxiters::Int = 100
    on_failure::F = ThrowError()
    backend::B = AutoEnzyme(; mode=Enzyme.Forward)
end

function allocate_canopy_convergence(convergence_model::AndersonCanopyConvergence, buffers, model, boundary_layer_model;
    site, vapour_pressure_equation, soil_hydraulic_model,
)
    layer_count = length(buffers.leaf.leaf_temperature)
    state_length = canopy_state_length(layer_count)
    (; history, backend, air_temperature_tolerance, vapor_density_tolerance, ground_temperature_tolerance,
       inverse_obukhov_length_tolerance) = convergence_model
    history >= 1 || throw(ArgumentError("history must be at least 1, got $history"))
    tolerances = canopy_tolerances(convergence_model)
    # Puts each step in units of its tolerance. A leaf's step is a temperature.
    step_scale = 1 ./ canopy_state_filled(layer_count, air_temperature_tolerance, air_temperature_tolerance,
        vapor_density_tolerance, ground_temperature_tolerance, inverse_obukhov_length_tolerance)
    lowest_temperature, highest_temperature = CANOPY_TEMPERATURE_BOUNDS
    lowest_state = canopy_state_filled(layer_count, lowest_temperature, lowest_temperature, 0.0u"g/m^3", lowest_temperature, -Inf * u"m^-1")
    highest_state = canopy_state_filled(layer_count, highest_temperature, highest_temperature, Inf * u"g/m^3", highest_temperature, Inf * u"m^-1")
    forcing = Ref(initial_residual_forcing(model, soil_hydraulic_model))
    # Prepared with a placeholder temperature and the first layer, to fix the argument types
    leaf_preparation = prepare_derivative(layer_leaf_balance_density, backend, 288.0,
        Constant(buffers), Constant(model), Constant(forcing[]), Constant(1))
    return (; forcing, tolerances, leaf_buffers=buffers, leaf_preparation, backend, step_scale, lowest_state, highest_state,
        residual=allocate_canopy_residual(layer_count), leaf_temperature_error=zeros(typeof(0.0u"K"), layer_count),
        state=zeros(state_length), best_state=zeros(state_length), last_good_state=zeros(state_length),
        next_state=zeros(state_length), step=zeros(state_length),
        previous_next_state=zeros(state_length), previous_step=zeros(state_length),
        next_state_changes=zeros(state_length, history), step_changes=zeros(state_length, history),
        history_matrix=zeros(history, history), history_weights=zeros(history))
end

# One layer's leaf energy balance per unit leaf area (W/m²) at `leaf_temperature` (K)
function layer_leaf_balance_density(leaf_temperature, buffers, model, forcing, layer)
    (; balance, leaf_area) = layer_leaf_balance(buffers, model, forcing, layer, leaf_temperature * CANOPY_STATE_UNITS.temperature)
    return in_state_units(balance.net_heat_flow / leaf_area)
end

# One pass: writes the state this pass would move to into `next_state`, and the move into `step`
function canopy_fixed_point!(next_state, step, state, buffers, model, boundary_layer_model, forcing;
    site, vapour_pressure_equation, ground_temperature_after_hour,
)
    (; residual, tolerances, leaf_temperature_error, leaf_buffers, leaf_preparation, backend) = buffers.convergence
    (; ground_temperature, inverse_obukhov_length) = unpack_canopy_state!(buffers, state)
    residual.ground_temperature[] = zero(ground_temperature)
    result = evaluate_canopy!(residual, inverse_obukhov_length, ground_temperature, buffers, model, boundary_layer_model, forcing;
        site, vapour_pressure_equation)
    # The soil can't take a canopy that overshot
    if isfinite(canopy_scaled_residual(residual, tolerances))
        residual.ground_temperature[] = ground_temperature -
            ground_temperature_after_hour(ground_temperature, result.longwave_result.ground_absorbed_longwave)
    end
    for layer in eachindex(leaf_temperature_error)
        if iszero(layer_leaf_areas(buffers, layer).exchange_fraction)
            # A layer without leaves takes the air temperature
            leaf_temperature_error[layer] = residual.leafless_temperature_difference[layer]
            continue
        end
        slope = derivative(layer_leaf_balance_density, leaf_preparation, backend, in_state_units(buffers.leaf.leaf_temperature[layer]),
            Constant(leaf_buffers), Constant(model), Constant(forcing), Constant(layer))
        leaf_temperature_error[layer] = in_state_units(residual.leaf_balance[layer]) / slope * CANOPY_STATE_UNITS.temperature
    end
    # The residuals of air, vapor, ground and stability are already the distance to move back
    pack_canopy_state!(step, leaf_temperature_error, residual.air_temperature, residual.vapor_density,
        residual.ground_temperature[], residual.inverse_obukhov_length[])
    next_state .= state .- step
    step .= .-step
    return next_state
end

# Finds, by least squares, the weight of each stored pass that best cancels the current step.
# Solves in-place by Cholesky, overwriting `matrix` and `weights`; returns false if the
# stored passes are too alike to separate (the system is not positive definite).
function solve_history_weights!(history_weights, history_matrix, step_changes, step, step_scale, stored_count)
    weights = view(history_weights, 1:stored_count)
    matrix = view(history_matrix, 1:stored_count, 1:stored_count)
    for first_column in 1:stored_count
        weights[first_column] = 0.0
        for index in eachindex(step)
            weights[first_column] += step_changes[index, first_column] * step[index] * step_scale[index]
        end
        for second_column in 1:stored_count
            matrix[first_column, second_column] = 0.0
            for index in eachindex(step)
                matrix[first_column, second_column] += step_changes[index, first_column] * step_changes[index, second_column]
            end
        end
        # A tiny increase on the diagonal keeps nearly identical passes solvable
        matrix[first_column, first_column] *= 1.0 + 1.0e-10
    end
    factorization = cholesky!(Symmetric(matrix); check=false)
    issuccess(factorization) || return false
    ldiv!(factorization, weights)
    return true
end

function converge_canopy!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs,
    initial_wind, convergence_model::AndersonCanopyConvergence; site, vapour_pressure_equation, soil_hydraulic_model,
    ground_temperature_after_hour::G,
) where G
    (; forcing, tolerances, state, best_state, last_good_state, residual, next_state, step, previous_next_state,
       previous_step, next_state_changes, step_changes, history_matrix, history_weights, step_scale,
       lowest_state, highest_state) = buffers.convergence
    (; damping, history, maxiters, on_failure) = convergence_model
    keywords = (; site, vapour_pressure_equation, soil_hydraulic_model)
    canopy_hour_setup!(state, forcing, buffers, model, boundary_layer_model, inputs, initial_wind; keywords...)

    stored_count = 0
    oldest_column = 1
    smallest_residual = Inf
    copyto!(best_state, state)
    for pass in 1:maxiters
        canopy_fixed_point!(next_state, step, state, buffers, model, boundary_layer_model, forcing[];
            site, vapour_pressure_equation, ground_temperature_after_hour)
        if !all(isfinite, next_state)
            pass == 1 && return canopy_failure!(on_failure, convergence_model, SciMLBase.ReturnCode.Unstable, best_state, pass,
                buffers, model, boundary_layer_model, forcing[], inputs; keywords...)
            # Overshot: go halfway back to the last good state, without the history that led here
            state .= (state .+ last_good_state) ./ 2
            stored_count = 0
            oldest_column = 1
            continue
        end
        copyto!(last_good_state, state)
        scaled_residual = canopy_scaled_residual(residual, tolerances)
        scaled_residual <= 1 && return canopy_solution!(residual, state, buffers, model, boundary_layer_model, forcing[],
            inputs, pass; keywords...)
        if scaled_residual < smallest_residual
            smallest_residual = scaled_residual
            copyto!(best_state, state)
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
            for index in eachindex(step)
                step_changes[index, column] = (step[index] - previous_step[index]) * step_scale[index]
                next_state_changes[index, column] = next_state[index] - previous_next_state[index]
            end
        end
        copyto!(previous_step, step)
        copyto!(previous_next_state, next_state)

        if stored_count > 0 && !solve_history_weights!(history_weights, history_matrix, step_changes, step, step_scale, stored_count)
            stored_count = 0
            oldest_column = 1
        end

        # Take the damped step, then subtract the weighted stored passes
        for index in eachindex(state)
            updated = state[index] + damping * step[index]
            for column in 1:stored_count
                updated -= (next_state_changes[index, column] -
                    (1.0 - damping) * step_changes[index, column] / step_scale[index]) * history_weights[column]
            end
            state[index] = clamp(updated, lowest_state[index], highest_state[index])
        end
    end
    return canopy_failure!(on_failure, convergence_model, SciMLBase.ReturnCode.MaxIters, best_state, maxiters,
        buffers, model, boundary_layer_model, forcing[], inputs; keywords...)
end
