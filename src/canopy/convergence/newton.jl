"""
    NewtonCanopyConvergence(; algorithm=NewtonRaphson(), abstol=1e-8, maxiters=50,
        backend=AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward)))

Solves the hourly canopy residual ([`evaluate_canopy!`](@ref)) for leaf and
air temperature, vapor density and inverse Obukhov length together, with a
NonlinearSolve `algorithm` and the residual's Jacobian from `backend`
through DifferentiationInterface. The solver cache is built once by
`allocate_canopy` and re-initialised each hour.
"""
@kwdef struct NewtonCanopyConvergence{A,T,B} <: AbstractCanopyConvergenceModel
    algorithm::A = NewtonRaphson()
    abstol::T = 1e-8
    maxiters::Int = 50
    # TODO drop set_runtime_activity with the reference-height pin in atmospheric_surface_profile! (#137)
    backend::B = AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward))
end

function allocate_canopy_convergence(convergence_model::NewtonCanopyConvergence, buffers, model, boundary_layer_model;
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
    parameters = Ref((; buffers, cache, constants, preparation, backend, residual))
    problem = NonlinearProblem(
        NonlinearFunction(canopy_nonlinear_residual!; jac=canopy_nonlinear_jacobian!, jac_prototype=zeros(n_state, n_state)),
        state, parameters)
    solver = SciMLBase.init(problem, convergence_model.algorithm; abstol=convergence_model.abstol, maxiters=convergence_model.maxiters)
    return (; solver, parameters, forcing, state, residual)
end

canopy_nonlinear_residual!(residual, state, parameters) =
    canopy_residual!(residual, state, parameters[].buffers, parameters[].constants.data)

function canopy_nonlinear_jacobian!(jacobian, state, parameters)
    (; buffers, cache, constants, preparation, backend, residual) = parameters[]
    jacobian!(canopy_residual!, residual, jacobian, preparation, backend, state, cache, constants)
    return nothing
end

function converge_canopy!(buffers, model::MultilayerCanopy, boundary_layer_model, inputs::CanopyEnergyBalanceInputs,
    initial_wind, cm::NewtonCanopyConvergence; site, vapour_pressure_equation, soil_hydraulic_model,
)
    (; environment_instant, zenith_angle, ground_temperature, ground_emissivity, ground_relative_humidity,
       leaf_water_potential) = inputs
    (; atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity) = environment_instant
    (; solver, parameters, forcing, state, residual) = buffers.convergence
    (; leaf_temperature, sensible_heat_source, potential_evaporation_mass_flow, leaf_body) = buffers.leaf
    air = buffers.air_profile
    (; layer_transmission) = buffers.longwave
    (; layer_plant_area_index_boundaries) = buffers.shortwave
    n_layers = length(leaf_temperature)

    forcing[] = CanopyResidualForcing(; atmospheric_pressure, reference_temperature, reference_wind_speed,
        reference_humidity, ground_temperature, ground_emissivity, ground_relative_humidity, leaf_water_potential,
        incoming_longwave=canopy_incoming_longwave(model; site, environment_instant, vapour_pressure_equation),
        dry_conductance=stomatal_conductance(model.stomatal_model, zenith_angle, leaf_water_potential, soil_hydraulic_model),
        canopy_top_air_temperature=initial_wind.canopy_top_air_temperature,
        canopy_top_relative_humidity=initial_wind.canopy_top_relative_humidity)

    # Warm start from last hour's solution; seed vapor density on the first call.
    if all(iszero, air.vapor_density)
        fill!(air.vapor_density, wet_air_properties(initial_wind.canopy_top_air_temperature,
            initial_wind.canopy_top_relative_humidity, atmospheric_pressure; vapour_pressure_equation).vapour_density)
    end
    for i in 1:n_layers
        state[i] = ustrip(u"K", leaf_temperature[i])
        state[n_layers + i] = ustrip(u"K", air.air_temperature[i])
        state[2n_layers + i] = ustrip(u"g/m^3", air.vapor_density[i])
    end
    state[3n_layers + 1] = ustrip(u"m^-1", initial_wind.inverse_obukhov_length)

    SciMLBase.reinit!(solver, state; p=parameters)
    SciMLBase.solve!(solver)
    SciMLBase.successful_retcode(solver.retcode) ||
        @warn "NewtonCanopyConvergence did not converge" solver.retcode maxlog=5

    # Leave buffers at the solution: the residual's own intermediates, with
    # the air profile at the solved state rather than its transport.
    result = evaluate_canopy!(residual, solver.u, buffers, model, boundary_layer_model, forcing[];
        site, vapour_pressure_equation)
    unstressed_dry_conductance = stomatal_conductance(model.stomatal_model, zenith_angle, 0.0u"J/kg", soil_hydraulic_model)
    for i in 1:n_layers
        air.air_temperature[i] = air.air_temperature_prev[i]
        air.vapor_density[i] = air.vapor_density_prev[i]
        saturation = wet_air_properties(air.air_temperature[i], 1.0, atmospheric_pressure; vapour_pressure_equation).vapour_density
        air.relative_humidity[i] = clamp(air.vapor_density[i] / saturation, 0.0, 1.0)

        # Unstressed (water_potential=0) reference transpiration, for a
        # soil-hydraulics demand term — not the leaf's own energy balance.
        exchange_fraction = 1.0 - layer_transmission[i]
        if iszero(exchange_fraction)
            potential_evaporation_mass_flow[i] = zero(potential_evaporation_mass_flow[i])
            continue
        end
        layer_pai = layer_plant_area_index_boundaries[i + 1] - layer_plant_area_index_boundaries[i]
        aerodynamic_area = 2.0 * layer_pai * 1.0u"m^2"
        conv = leaf_convection(model.leaf_convection_model, leaf_body, aerodynamic_area, air.air_temperature[i],
            leaf_temperature[i], buffers.wind.wind_speed[i], atmospheric_pressure)
        atmos = AtmosphericConditions(air.relative_humidity[i], buffers.wind.wind_speed[i], atmospheric_pressure)
        potential_evaporation_mass_flow[i] = HeatExchange.evaporation(unstressed_dry_conductance,
            conv.mass_transfer_coefficient, atmos, aerodynamic_area, leaf_temperature[i], air.air_temperature[i];
            water_potential=0.0u"J/kg").transpiration_mass_flow
    end

    return result.longwave_result, solver.stats.nsteps, result.canopy_top_air_temperature,
        result.canopy_top_relative_humidity, result.friction_velocity, result.inverse_obukhov_length
end
