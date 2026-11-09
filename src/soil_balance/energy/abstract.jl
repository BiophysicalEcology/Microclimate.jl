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

Allocate per-layer scratch buffers (`layer_depths`, `heat_capacity`,
`thermal_conductance`, warm-started Obukhov length) for the energy ODE.
"""
function allocate_soil_energy_balance end

"""
    GroundAirExchange(; air_temperature, air_relative_humidity, wind_speed, height,
        heat_conductance, vapor_conductance)

The air at `height` just above the ground, and how the ground exchanges heat and water
vapor with it. A canopy gives this to the soil. Use the conductances as they are.

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
