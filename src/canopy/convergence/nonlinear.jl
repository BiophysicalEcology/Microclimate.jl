"""
    NonlinearCanopyConvergence(; algorithm=NewtonRaphson(), fallback_algorithm=NewtonRaphson(; linesearch=BackTracking()),
        abstol=1e-8, maxiters=50, backend=AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward)))

Solves the hourly canopy residual ([`evaluate_canopy!`](@ref)) for leaf and
air temperature, vapor density and inverse Obukhov length together, with a
NonlinearSolve `algorithm` and the residual's Jacobian from `backend`
through DifferentiationInterface. The solver cache is built once by
`allocate_canopy` and re-initialized each hour.

`algorithm` is any in-place NonlinearSolve algorithm, e.g. `NewtonRaphson()`
or `Broyden(; init_jacobian=Val(:true_jacobian))` to compute the Jacobian
once per hour. When it fails, `fallback_algorithm` is tried from the previous
hour's solution, then from the reference conditions; its line search shortens
steps that leave [`CANOPY_TEMPERATURE_BOUNDS`](@ref).
"""
@kwdef struct NonlinearCanopyConvergence{A,F,T,B} <: AbstractCanopyConvergenceModel
    algorithm::A = NewtonRaphson()
    fallback_algorithm::F = NewtonRaphson(; linesearch=BackTracking())
    abstol::T = 1e-8
    maxiters::Int = 50
    # TODO drop set_runtime_activity with the reference-height pin in atmospheric_surface_profile! (#137)
    backend::B = AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward))
end

function allocate_canopy_convergence(convergence_model::NonlinearCanopyConvergence, buffers, model, boundary_layer_model;
    site, vapour_pressure_equation, soil_hydraulic_model,
)
    n_state = canopy_state_length(buffers)
    n_layers = length(buffers.leaf.leaf_temperature)
    forcing = Ref(initial_residual_forcing(model, soil_hydraulic_model))
    state = vcat(fill(288.0, 2n_layers), fill(5.0, n_layers), 0.0)
    residual = zeros(n_state)
    (; backend) = convergence_model
    cache = Cache(buffers)
    constants = Constant((model, boundary_layer_model, forcing, site, vapour_pressure_equation))
    preparation = prepare_jacobian(canopy_residual!, residual, backend, state, cache, constants)
    parameters = Ref((; buffers, cache, constants, preparation, backend, residual))
    problem = NonlinearProblem(
        NonlinearFunction(canopy_nonlinear_residual!; jac=canopy_nonlinear_jacobian!, jac_prototype=zeros(n_state, n_state)),
        state, parameters)
    (; abstol, maxiters) = convergence_model
    solver = SciMLBase.init(problem, convergence_model.algorithm; abstol, maxiters)
    fallback_solver = SciMLBase.init(problem, convergence_model.fallback_algorithm; abstol, maxiters)
    return (; solver, fallback_solver, parameters, forcing, state, warm_start=similar(state), residual, has_jacobian=Ref(false))
end

canopy_nonlinear_residual!(residual, state, parameters) =
    canopy_residual!(residual, state, parameters[].buffers, parameters[].constants.data)

function canopy_nonlinear_jacobian!(jacobian, state, parameters)
    (; buffers, cache, constants, preparation, backend, residual) = parameters[]
    jacobian!(canopy_residual!, residual, jacobian, preparation, backend, state, cache, constants)
    return nothing
end

function converge_canopy!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs,
    initial_wind, cm::NonlinearCanopyConvergence; site, vapour_pressure_equation, soil_hydraulic_model,
)
    (; atmospheric_pressure, reference_temperature, reference_humidity) = inputs.environment_instant
    (; solver, fallback_solver, parameters, forcing, state, warm_start, residual, has_jacobian) = buffers.convergence
    n_layers = length(buffers.leaf.leaf_temperature)

    canopy_hour_setup!(state, forcing, buffers, model, inputs, initial_wind;
        site, vapour_pressure_equation, soil_hydraulic_model)

    copyto!(warm_start, state)
    SciMLBase.reinit!(solver, state; p=parameters)
    # First step with the previous hour's Jacobian; restart with a fresh one if it fails
    if has_jacobian[]
        SciMLBase.step!(solver; recompute_jacobian=false)
        solver.retcode == SciMLBase.ReturnCode.Default || SciMLBase.successful_retcode(solver.retcode) ||
            SciMLBase.reinit!(solver, state; p=parameters)
    end
    SciMLBase.solve!(solver)
    has_jacobian[] = true
    iterations = solver.stats.nsteps
    if SciMLBase.successful_retcode(solver.retcode)
        copyto!(state, solver.u)
    else
        # Fallback from the warm start, then from the reference conditions
        SciMLBase.reinit!(fallback_solver, warm_start; p=parameters)
        SciMLBase.solve!(fallback_solver)
        iterations += fallback_solver.stats.nsteps
        if !SciMLBase.successful_retcode(fallback_solver.retcode)
            reference_vapor_density = ustrip(u"g/m^3", wet_air_properties(reference_temperature, reference_humidity,
                atmospheric_pressure; vapour_pressure_equation).vapour_density)
            for i in 1:n_layers
                state[i] = ustrip(u"K", reference_temperature)
                state[n_layers + i] = ustrip(u"K", reference_temperature)
                state[2n_layers + i] = reference_vapor_density
            end
            state[3n_layers + 1] = 0.0
            SciMLBase.reinit!(fallback_solver, state; p=parameters)
            SciMLBase.solve!(fallback_solver)
            iterations += fallback_solver.stats.nsteps
            SciMLBase.successful_retcode(fallback_solver.retcode) ||
                @warn "NonlinearCanopyConvergence did not converge" fallback_solver.retcode maxlog=5
        end
        copyto!(state, fallback_solver.u)
    end

    return canopy_solution!(residual, state, buffers, model, boundary_layer_model, forcing[], inputs, iterations;
        site, vapour_pressure_equation, soil_hydraulic_model)
end
