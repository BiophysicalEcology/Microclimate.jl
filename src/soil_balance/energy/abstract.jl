"""
    SoilHeatTransportModel

Supertype for soil column energy balance models. [`SoilHeatTransport1D`](@ref) is the
only implementation.
"""
abstract type SoilHeatTransportModel end

"""
    soil_energy_balance(model, temperature_state, p, t)

ODE right hand side for the soil temperature column: returns the per-node `dT/dt`
SVector. Variants supply the energy budget formulation (conduction +
surface boundary condition: net radiation + convection + evaporation + phase corrections).
"""
function soil_energy_balance end

"""
    allocate_soil_energy_balance(model, num_nodes)

Allocate the buffers for the energy ODE.
"""
function allocate_soil_energy_balance end

"""
    AbstractGroundTemperatureModel

How the soil gives the ground temperature at the end of an hour for a surface boundary.
Concrete models implement [`ground_temperature`](@ref):
[`IntegratedGroundTemperature`](@ref), [`PredictedGroundTemperature`](@ref) and [`HourStartGroundTemperature`](@ref).
"""
abstract type AbstractGroundTemperatureModel end

"""
    HourStartGroundTemperature()

Keeps the ground at its temperature at the start of the hour. The canopy is solved at that
temperature, and the soil is integrated once, with the canopy's solved boundary.
"""
struct HourStartGroundTemperature <: AbstractGroundTemperatureModel end

"""
    IntegratedGroundTemperature()

Integrates the whole soil column over the hour.
"""
struct IntegratedGroundTemperature <: AbstractGroundTemperatureModel end

"""
    PredictedGroundTemperature(; nodes=3, solver=NewtonRaphson(), tolerance=1e-3u"K")

One backward Euler step over the hour for the surface node and the `nodes - 1` nodes below it.
The deeper nodes keep their temperatures from the start of the hour. `solver`, a NonlinearSolve
algorithm, solves the step to within `tolerance`.
"""
struct PredictedGroundTemperature{N,S,T} <: AbstractGroundTemperatureModel
    solver::S
    tolerance::T
end
function PredictedGroundTemperature(; nodes::Int=3, solver=NewtonRaphson(), tolerance=1e-3u"K")
    nodes >= 1 || throw(ArgumentError("PredictedGroundTemperature needs at least 1 node, got $nodes"))
    return PredictedGroundTemperature{nodes,typeof(solver),typeof(tolerance)}(solver, tolerance)
end

"""
    ground_temperature(model::AbstractGroundTemperatureModel, solver, integrator, hour_state, hour_index; surface_node)

The temperature of node `surface_node` at the end of hour `hour_index`, starting from
`hour_state`. The surface boundary is the one in `integrator`'s inputs. `solver` is
`model`'s own solver, from `allocate_ground_temperature`.
"""
function ground_temperature end

"""
    GroundAirExchange(; air_temperature, air_relative_humidity, wind_speed, height,
        heat_conductance, vapor_conductance)

The air at `height` just above the ground, and how the ground exchanges heat and water
vapor with it. A canopy gives this to the soil, and the soil uses these conductances in place of its own surface exchange.

- Sensible heat flux up from the ground = `heat_conductance * (ground temperature - air_temperature)`.
- Water vapor flux up from the ground = `vapor_conductance * (ground vapor density - air vapor density)`.
"""
@kwdef struct GroundAirExchange{AT,RH,WS,H,HC,VC}
    air_temperature::AT
    air_relative_humidity::RH
    wind_speed::WS
    height::H
    heat_conductance::HC
    vapor_conductance::VC
end

"""
    GroundSurface(; temperature, emissivity, reflectance, wetness)

The ground surface a canopy sees: its temperature, longwave emissivity, shortwave
reflectance and wetness. The reflectance is one number, or a [`GroundShortwaveReflectance`](@ref)
for the unscattered beam and diffuse light. The `wetness` is the fraction of the surface that
evaporates, as the soil's own evaporation takes it: water vapor flux up from the ground =
`wetness * vapor_conductance * (saturated vapor density at the ground - air vapor density)`.
"""
@kwdef struct GroundSurface{T,E,R,W}
    temperature::T
    emissivity::E
    reflectance::R
    wetness::W
end
