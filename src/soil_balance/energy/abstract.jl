abstract type AbstractSoilEnergyScheme end

"""
    soil_energy_balance(scheme, temperature_state, p, t)

ODE rhs for the soil temperature column: returns the per-node `dT/dt`
SVector. Variants supply the energy budget formulation (conduction +
surface BC: net radiation + convection + evaporation + phase corrections).
"""
function soil_energy_balance end

"""
    allocate_soil_energy_balance(scheme, num_nodes)

Allocate per-layer scratch buffers (`layer_depths`, `heat_capacity`,
`thermal_conductance`, warm-started Obukhov length) for the energy ODE.
"""
function allocate_soil_energy_balance end
