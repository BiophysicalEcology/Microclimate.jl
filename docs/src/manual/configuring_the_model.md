# Configuring the model

Every physical process in Microclimate.jl — soil thermal properties, soil hydraulics,
radiation, snow, canopy, evaporation, condensation, convergence — is a separate,
swappable model type. A simulation is built by choosing one implementation for each
process and composing them into a [`MicroModel`](@ref).

## Composability and memory management, a cooking metaphor

Microclimate calculations are computationally intensive and Microclimate.jl has been 
designed to be highly performant. Part of the performance comes from memory managemnt,
reducing the need for repeated memory allocations. [`init`](@ref) builds a [`MicroCache`](@ref) once: 
it allocates every buffer the solve loop touches — soil ODE workspace, snow scratch, 
canopy buffers, the solar radiation output — sized for the chosen model. 

A useful metaphor is setting up a kitchen for cooking cakes. The mixing bowls, oven, pans, to 
hold the required ingredients are all laid out once via the cache. [`solve!`](@ref) then 
bakes a cake of the type you requested with that tailored kitchen. [`reinit!`](@ref) 
swaps in new [`MicroInputs`](@ref) — a different site, a different year of weather — 
and bakes another cake without resetting up the kitchen:

```julia
cache = init(problem)
out1 = solve!(cache)          # first cake

reinit!(cache, other_inputs)
out2 = solve!(cache)          # same kitchen, different cake
```

Changing the *model* — say, swapping `snow_model = NoSnow()` for `SnowModel(...)` (frosting!), 
or `soil_hydraulic_model` for a different infiltration algorithm — changes what buffers
are needed, so that goes through a fresh `MicroProblem`/`init`, not `reinit!`.

## Why compositional

`MicroModel` is deliberately compositional: each physical process is a type
implementing a small interface
(an `Abstract*Model` supertype and a handful of methods), and swapping one out is a
single keyword argument. This is meant to let the package grow — a new snow-density
formula, a new infiltration algorithm, a new leaf-temperature solver is a new type
implementing the existing interface, not a rewrite. It also makes model *choices*
cheap to interrogate: because Julia's multiple dispatch selects behaviour from the
argument types rather than from branches inside one function, comparing "what changes
if I use `MatricFluxPotentialAlgorithm` instead of `MatricPotentialAlgorithm`" or
"what does `KTheoryAirProfile` predict versus `RaupachLTheoryAirProfile`" is a
one-line change to a `MicroModel`/`MicroConfig`, not a fork of the simulation.
[Soil hydraulics algorithms](../tutorials/soil_hydraulics_algorithms.md) demonstrates
exactly this kind of comparison.

## `MicroModel`

[`MicroModel`](@ref) is the constant-across-runs description of the simulation: the
site, and one model per part of the system, each holding its own sub-models. Not all
slots have alternatives at present but the idea is to allow a community of developers to
build new functionality over time where needed.

| Slot | Holds | Default |
|:---|:---|:---|
| `hours` | hours of the day solved | 0–23 |
| `site` | [`Site`](@ref): location, terrain, albedo, roughness | — (required) |
| `soil` | [`SoilModel`](@ref): `depths`, `profile`, `properties_model`, `hydraulic_model`, `energy_model`, `moisture_strategy` | — (required) |
| `surface` | [`SurfaceModel`](@ref): `evaporation_model`, `condensation_model`, `max_surface_pool` | `SurfaceModel()` |
| `atmosphere` | [`AtmosphereModel`](@ref): `heights`, `radiation`, `boundary_layer_model` | `AtmosphereModel()` |
| `snow` | [`NoSnow`](@ref) or [`SnowModel`](@ref) | `NoSnow()` |
| `canopy` | [`NoCanopy`](@ref) or [`MultilayerCanopy`](@ref) | `NoCanopy()` |
| `vapour_pressure_equation` | `GoffGratch`, `Teten` or `Huang` (from FluidProperties.jl) | `GoffGratch()` |
| `config` | [`MicroConfig`](@ref): iteration and data delivery | `MicroConfig()` |

The soil's `properties_model` ([`CampbelldeVriesSoilProperties`](@ref)) and
`hydraulic_model` ([`CampbellSoilHydraulics`](@ref)) are required; its `energy_model`
defaults to [`SoilHeatTransport1D`](@ref) and its `moisture_strategy` to
[`PrescribedSoilMoisture`](@ref) (from `environment_daily`), with
[`DynamicSoilMoisture`](@ref) solving it from rainfall. The surface's
`condensation_model` defaults to [`GarrattSegalCondensation`](@ref), with
[`BulkTransferCondensation`](@ref) and [`NoCondensation`](@ref) as alternatives.
`radiation`'s own three slots (`solar_radiation_model`, `longwave_model`,
`shortwave_model`) are documented in [Radiation](radiation.md), and the canopy's in
[Canopy](canopy.md).

## `MicroConfig`

[`MicroConfig`](@ref) is the "how to iterate, how data is delivered" side of the
model, separate from the physical-process choices above:

- `convergence` — a reusable [`FixedIterationConvergence`](@ref)/
  [`IterationToleranceConvergence`](@ref) strategy. How much it actually controls
  depends on `MicroProblem.time_mode`: for [`NonConsecutiveDayMode`](@ref) (the
  default) it's ignored — that mode's own `iterations_per_day` field sets the
  per-day iteration count instead; for [`ConsecutiveDayMode`](@ref) it only applies
  on day 1, and only if `spinup_first_day=true`.
- `rainfall_schedule` — [`DailyRainfall`](@ref) or [`HourlyRainfall`](@ref)
- `canopy_soil_convergence` — how the canopy and soil-heat ODE are coupled each hour:
  [`IterationToleranceConvergence`](@ref) (the default) iterates them to agreement;
  [`SurfaceNodePrediction`](@ref) integrates the soil once per hour
  (see [Canopy](canopy.md#Solving-each-hour))
- `canopy_soil_relaxation` — under-relaxation of the iterated coupling

## How the pieces fit together each hour

The physical processes above are not all solved as one monolithic coupled system but 
with different aspects converging in different ways via different algorithsm:

- **Soil temperature** is an ODE system ([`SoilHeatTransport1D`](@ref), handed
  to a SciML integrator) — continuous within the hour, not a single explicit step.
- **Canopy**, when present, is jointly converged against that same ODE each hour.
  With [`IterationToleranceConvergence`](@ref) the canopy-soil iteration repeatedly
  (a) solves the canopy's shortwave, longwave, wind and leaf-temperature state (with
  the canopy's `convergence_model`) for the current guess of ground temperature,
  (b) feeds the resulting ground-surface boundary condition into the soil ODE and
  advances it through the hour, (c) compares the new ground temperature to the
  previous guess, and repeats until `canopy_soil_convergence` is satisfied (or its
  iteration limit is hit). [`SurfaceNodePrediction`](@ref) instead predicts the
  surface temperature between canopy solves and integrates the ODE once. Either way,
  canopy and soil *temperature* reach a mutually consistent state each hour.
- **Soil moisture**, ground dew/frost, and snow growth/melt are then computed once
  per hour, downstream of that hour's already-converged temperature — not solved
  simultaneously with it. Soil moisture in particular assumes the current hour's
  temperature as given, rather than being part of the same implicit system; the
  *next* hour's temperature ODE then reads back the moisture-dependent soil thermal
  properties ([Soil thermal properties](soil_thermal_properties.md)), so temperature
  and moisture are coupled at the hourly cadence, just not simultaneously within it.

## `MicroInputs` and `MicroProblem`

[`MicroInputs`](@ref) is the per-run data: [`Site`](@ref), [`SoilProfile`](@ref),
environment forcings and initial conditions. [`MicroProblem`](@ref) pairs a
`MicroModel` with a `MicroInputs`, the days to simulate, and a
[`AbstractTimeMode`](@ref) ([`NonConsecutiveDayMode`](@ref) for independent
representative days, [`ConsecutiveDayMode`](@ref) for a continuous run where state
carries from day to day).

## References

See the [References](references.md) page for the full bibliography.
