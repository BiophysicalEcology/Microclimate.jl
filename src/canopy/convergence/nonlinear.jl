"""
    NonlinearCanopyConvergence(; algorithm=NewtonRaphson(), fallback_algorithm=NewtonRaphson(; linesearch=BackTracking()),
        leaf_balance_tolerance=2e-5u"W/m^2", air_temperature_tolerance=1e-6u"K",
        vapor_density_tolerance=1e-6u"g/m^3", inverse_obukhov_length_tolerance=1e-8u"m^-1",
        maxiters=50, backend=AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward)))

Solves the hourly canopy residual ([`evaluate_canopy!`](@ref)) for leaf and
air temperature, vapor density and inverse Obukhov length together, with a
NonlinearSolve `algorithm` and the residual's Jacobian from `backend`
through DifferentiationInterface. The solver cache is built once by
`allocate_canopy` and re-initialized each hour.

The solve converges when every residual row is within the tolerance for its kind:
leaf energy balance, air temperature, vapor density and inverse Obukhov length.
The leaf row of a layer without leaf area is a temperature, so it takes
`air_temperature_tolerance`.

`algorithm` is any in-place NonlinearSolve algorithm, e.g. `NewtonRaphson()`
or `Broyden(; init_jacobian=Val(:true_jacobian))` to compute the Jacobian
once per hour. When it fails, `fallback_algorithm` is tried from the previous
hour's solution, then from the reference conditions; its line search shortens
steps that leave [`CANOPY_TEMPERATURE_BOUNDS`](@ref). If that fails too, it errors.
"""
@kwdef struct NonlinearCanopyConvergence{A,F,LB,AT,VD,IL,B} <: AbstractCanopyConvergenceModel
    algorithm::A = NewtonRaphson()
    fallback_algorithm::F = NewtonRaphson(; linesearch=BackTracking())
    leaf_balance_tolerance::LB = 2e-5u"W/m^2"
    air_temperature_tolerance::AT = 1e-6u"K"
    vapor_density_tolerance::VD = 1e-6u"g/m^3"
    inverse_obukhov_length_tolerance::IL = 1e-8u"m^-1"
    maxiters::Int = 50
    # TODO drop set_runtime_activity with the reference-height pin in atmospheric_surface_profile! (#137)
    backend::B = AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward))
end

function allocate_canopy_convergence(convergence_model::NonlinearCanopyConvergence, buffers, model, boundary_layer_model;
    site, vapour_pressure_equation, soil_hydraulic_model,
)
    n_state = canopy_state_length(buffers)
    n_layers = length(buffers.leaf.leaf_temperature)
    forcing = Ref(CanopyResidualForcing(;
        atmospheric_pressure=101325.0u"Pa", reference_temperature=288.0u"K", reference_wind_speed=1.0u"m/s",
        reference_humidity=0.5, incoming_longwave=300.0u"W/m^2", ground_temperature=288.0u"K",
        ground_emissivity=0.95, ground_relative_humidity=0.5, leaf_water_potential=0.0u"J/kg",
        dry_conductance=stomatal_conductance(model.stomatal_model, 0.0u"°", 0.0u"J/kg", soil_hydraulic_model),
        canopy_top_air_temperature=288.0u"K", canopy_top_relative_humidity=0.5))
    state = vcat(fill(288.0, 2n_layers), fill(5.0, n_layers), 0.0)
    residual = zeros(n_state)
    (; backend) = convergence_model
    cache = Cache(buffers)
    constants = Constant((model, boundary_layer_model, forcing, site, vapour_pressure_equation))
    preparation = prepare_jacobian(canopy_residual!, residual, backend, state, cache, constants)
    # In the state's units (see canopy_state_length)
    tolerances = (; leaf_balance=ustrip(u"W/m^2", convergence_model.leaf_balance_tolerance),
        air_temperature=ustrip(u"K", convergence_model.air_temperature_tolerance),
        vapor_density=ustrip(u"g/m^3", convergence_model.vapor_density_tolerance),
        inverse_obukhov_length=ustrip(u"m^-1", convergence_model.inverse_obukhov_length_tolerance))
    parameters = Ref((; buffers, cache, constants, preparation, backend, residual, tolerances))
    problem = NonlinearProblem(
        NonlinearFunction(canopy_nonlinear_residual!; jac=canopy_nonlinear_jacobian!, jac_prototype=zeros(n_state, n_state)),
        state, parameters)
    # The solver sees each row over its tolerance, so 1 is every row within tolerance
    (; maxiters) = convergence_model
    solver = SciMLBase.init(problem, convergence_model.algorithm; abstol=1.0, maxiters)
    fallback_solver = SciMLBase.init(problem, convergence_model.fallback_algorithm; abstol=1.0, maxiters)
    return (; solver, fallback_solver, parameters, forcing, state, warm_start=similar(state), residual, has_jacobian=Ref(false))
end

# The tolerance of residual row `row`, in the state's units
function canopy_row_tolerance(tolerances, buffers, row)
    n_layers = length(buffers.leaf.leaf_temperature)
    if row <= n_layers
        # A layer without leaf area gives T_leaf - T_air
        iszero(1.0 - buffers.longwave.layer_transmission[row]) && return tolerances.air_temperature
        return tolerances.leaf_balance
    elseif row <= 2n_layers
        return tolerances.air_temperature
    elseif row <= 3n_layers
        return tolerances.vapor_density
    else
        return tolerances.inverse_obukhov_length
    end
end

function canopy_nonlinear_residual!(residual, state, parameters)
    (; buffers, constants, tolerances) = parameters[]
    canopy_residual!(residual, state, buffers, constants.data)
    for row in eachindex(residual)
        residual[row] /= canopy_row_tolerance(tolerances, buffers, row)
    end
    return nothing
end

function canopy_nonlinear_jacobian!(jacobian, state, parameters)
    (; buffers, cache, constants, preparation, backend, residual, tolerances) = parameters[]
    jacobian!(canopy_residual!, residual, jacobian, preparation, backend, state, cache, constants)
    for row in axes(jacobian, 1)
        tolerance = canopy_row_tolerance(tolerances, buffers, row)
        for column in axes(jacobian, 2)
            jacobian[row, column] /= tolerance
        end
    end
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
    iterations = 0
    # First step with the previous hour's Jacobian; restart with a fresh one if it fails
    if has_jacobian[]
        SciMLBase.step!(solver; recompute_jacobian=false)
        if !(solver.retcode == SciMLBase.ReturnCode.Default || SciMLBase.successful_retcode(solver.retcode))
            iterations += solver.stats.nsteps
            SciMLBase.reinit!(solver, state; p=parameters)
        end
    end
    SciMLBase.solve!(solver)
    has_jacobian[] = SciMLBase.successful_retcode(solver.retcode)
    iterations += solver.stats.nsteps
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
        end
        SciMLBase.successful_retcode(fallback_solver.retcode) ||
            error("NonlinearCanopyConvergence did not converge from the previous hour or the reference conditions: ",
                fallback_solver.retcode, " (reference temperature ", reference_temperature, ", wind ", reference_wind_speed,
                ", ground temperature ", ground_temperature, ", zenith angle ", zenith_angle, ")")
        copyto!(state, fallback_solver.u)
    end

    return canopy_solution!(residual, state, buffers, model, boundary_layer_model, forcing[], inputs, iterations;
        site, vapour_pressure_equation, soil_hydraulic_model)
end
