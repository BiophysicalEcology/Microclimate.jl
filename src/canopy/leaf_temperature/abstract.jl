"""
    AbstractLeafTemperatureSolver

How [`leaf_temperature`](@ref) resolves leaf temperature from the leaf
energy balance (absorbed radiation vs. emitted longwave, convection, and
transpiration — see `leaf_heat_balance`). A canopy model holds one of these
in its `leaf_temperature_solver` field.

Concrete variants: [`LinearizedLeafTemperature`](@ref) (default; cheap,
single linearised step, no bracket) and [`RootFindLeafTemperature`](@ref)
(exact root-find each call, more expensive).

Mirrors HeatExchange.jl's own `SurfaceSolveStrategy`/`LinearizedSurface`/
`RootFindSurface` (defined on an unmerged branch there) — kept local for now
to avoid depending on unreleased API; TODO: switch to HeatExchange.jl's
types once that branch merges to main.
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
