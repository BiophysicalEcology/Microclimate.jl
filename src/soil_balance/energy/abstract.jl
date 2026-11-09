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
[`SoilIntegration`](@ref) and [`SurfaceNodePrediction`](@ref).
"""
abstract type AbstractGroundTemperatureModel end

"""
    SoilIntegration()

Integrates the whole soil column over the hour.
"""
struct SoilIntegration <: AbstractGroundTemperatureModel end

"""
    SurfaceNodePrediction()

One backward Euler step over the hour for the surface node alone. The nodes below it keep
their temperatures from the start of the hour.
"""
struct SurfaceNodePrediction <: AbstractGroundTemperatureModel end

"""
    ground_temperature(model::AbstractGroundTemperatureModel, integrator, hour_state, hour_index; surface_node)

The temperature of node `surface_node` at the end of hour `hour_index`, starting from
`hour_state`. The surface boundary is the one in `integrator`'s inputs.
"""
function ground_temperature end

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

"""
    GroundSurface(; temperature, emissivity, reflectance, relative_humidity)

The ground surface a canopy sees: its temperature, longwave emissivity, shortwave
reflectance and the relative humidity of the air at the surface.
"""
@kwdef struct GroundSurface{T,E,R,RH}
    temperature::T
    emissivity::E
    reflectance::R
    relative_humidity::RH
end
