"""
    AbstractLeafTemperatureSolver

How [`leaf_temperature`](@ref) resolves leaf temperature from the leaf
energy balance (see `leaf_heat_balance`).

Concrete variants: [`LinearizedLeafTemperature`](@ref) (default; cheap,
single linearised step, no bracket) and [`RootFindLeafTemperature`](@ref)
(exact root-find each call, more expensive).

Mirrors HeatExchange.jl's own `SurfaceSolveStrategy` (on an unmerged
branch there) — kept local until that branch merges to main.
"""
abstract type AbstractLeafTemperatureSolver end

"""
    leaf_temperature(solver, absorbed_radiation, air_temperature, relative_humidity,
                      wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance,
                      leaf_water_potential, body; leaf_area, leaf_temperature_guess)

Leaf temperature from the energy balance in `leaf_heat_balance`, solved by
`solver::AbstractLeafTemperatureSolver`.
"""
function leaf_temperature end
