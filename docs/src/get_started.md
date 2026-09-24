# Get started

Microclimate.jl computes hourly above- and below-ground microclimates — air
temperature, wind speed, humidity, radiation, soil temperature, soil moisture, snow, dew, frost and
canopy state — from site, weather and soil inputs. Inputs and outputs are
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantities. It is part of the [BiophysicalEcology org](https://github.com/BiophysicalEcology), an ecosystem of packages for computing environment-organism interactions.

```julia
using Pkg
Pkg.add(url = "https://github.com/BiophysicalEcology/Microclimate.jl")
```

## The building blocks

A Microclimate.jl run needs four things: a [`MicroModel`](@ref) describing which physical process
models to use, a [`Site`](@ref), a [`SoilProfile`](@ref), and a set of
environment forcings. Microclimate.jl ships `example_*` constructors for all of these, built
around a default site (Madison, Wisconsin, USA):

```@example get_started
using Microclimate, Unitful

site = example_site()
soil_profile = example_soil_profile()
environment_minmax = example_monthly_weather()
environment_daily = example_daily_environment()
environment_hourly = example_hourly_environment()
nothing # hide
```

The model is presently restricted to solve on an hourly time step, so every forcing is resolved 
to this cadence, regardless of its input resolution: `environment_minmax`'s daily min/max values are
expanded into an hourly diel cycle (see [Diel curves and time handling](manual/diel_time_handling.md)),
`environment_daily`'s one-value-per-day forcings (shade, rainfall, leaf area index, …)
are held constant across that day's 24 hours, and `environment_hourly` supplies
already-hourly values directly.

## The model

A [`MicroModel`](@ref) collects the physical process models — including soil thermal
properties, soil hydraulics, radiation, snow, canopy algorithms — and the hours, depths and heights
across which to compute. By default it uses
[`CampbelldeVriesSoilProperties`](@ref) for soil thermal properties,
[`CampbellSoilHydraulics`](@ref) for water infiltration, [`NoSnow`](@ref) and
[`NoCanopy`](@ref):

```@example get_started
model = MicroModel(;
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model = example_soil_hydraulic_model(),
)
nothing # hide
```

See [Configuring the model](manual/configuring_the_model.md) for what each of
`MicroModel`'s slots does and how to swap in a different model.

## Assembling and solving a problem

[`MicroInputs`](@ref) pairs the site, soil profile and environment with initial
conditions; [`MicroProblem`](@ref) pairs a model with a set of inputs, the days to
simulate and how each day is connected (middle day of each month solved non-consecutively 
by default):

```@example get_started
depths = model.depths
inputs = MicroInputs(;
    site, soil_profile, environment_minmax, environment_daily, environment_hourly,
    initial_soil_temperature = nothing,  # each representative day resets to its own mean air temperature
    initial_soil_moisture = fill(0.105, length(depths)),
)
problem = MicroProblem(model, inputs)
out = solve(problem)
```

`out` is a `MicroResult` — matrices and vectors of hourly, per-depth state.
The simplest way to get all of the above at once is [`example_microclimate_problem`](@ref):

```@example get_started
out = solve(example_microclimate_problem())
nothing # hide
```

## Reading the output

Soil temperature at the surface and at the deepest simulated depth, for the first
representative day:

```@example get_started
first_day = @view out.soil_temperature[1:24, :]
(surface = first_day[:, 1], deep = first_day[:, end])
```

Air temperature and wind speed at each of `model.heights`:

```@example get_started
out.profile.air_temperature[1:3, :]
```

## Re-solving with different inputs

`init`/`solve!` split memory allocation from solving, so the same [`MicroCache`](@ref) can be
re-solved for different input data without reallocating:

```@example get_started
cache = init(problem)
out1 = solve!(cache)
inputs = MicroInputs(;
    site, soil_profile, environment_minmax, environment_daily, environment_hourly,
    initial_soil_temperature = nothing,  # must match the type used to build `problem`/`cache`
    initial_soil_moisture = fill(0.2, length(depths)),
)
reinit!(cache, inputs)
out2 = solve!(cache)
nothing # hide
```

`reinit!` reuses the cache's existing memory, so the new `MicroInputs` must have the same
field types as the one `cache` was built from — switching `initial_soil_temperature`
between `nothing` and an explicit vector between calls isn't allowed; `initial_soil_moisture`
can vary freely here since its type didn't change.

## Next steps

- [Introduction](manual/introduction.md) — what a microclimate is, and how this
  package relates to NicheMapR and other microclimate/land-surface models.
- [Configuring the model](manual/configuring_the_model.md) — the compositional
  design behind `MicroModel`.
- [Monthly workflow](tutorials/monthly_workflow.md) — driving a full year of monthly
  representative days from climate-normal data.
- [Daily/hourly workflow](tutorials/daily_hourly_workflow.md) — driving the model
  from hourly weather-station data.
