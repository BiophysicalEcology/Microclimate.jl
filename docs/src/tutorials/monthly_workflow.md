# Monthly workflow

One way to run the model is with long-term monthly climate normals:
[`MonthlyMinMaxEnvironment`](@ref) supplies one representative day per month (day of
year 15, 46, 74, …, `DEFAULT_DAYS`), each expanded to an hourly cycle from
that month's minimum/maximum values (see [Diel curves and time handling](../manual/diel_time_handling.md)).
This is driven here by [`example_monthly_weather`](@ref)'s Madison, Wisconsin climate
normals (from the CRU CL v2.0 dataset).

## Building the run

```@example monthly
using Microclimate, Unitful, CairoMakie

site = example_site()
soil_profile = example_soil_profile()
model = MicroModel(;
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model = example_soil_hydraulic_model(),
)
depths = model.depths

inputs = MicroInputs(;
    site, soil_profile,
    environment_minmax = example_monthly_weather(),
    environment_daily = example_daily_environment(),
    environment_hourly = example_hourly_environment(),
    initial_soil_temperature = fill(u"K"(7.74u"°C"), length(depths)),
    initial_soil_moisture = fill(0.42 * 0.25, length(depths)),
)
problem = MicroProblem(model, inputs)
out = solve(problem)
size(out.soil_temperature)
```

Each of the 12 representative days is iterated 3 times by default
([`NonConsecutiveDayMode`](@ref)'s own `iterations_per_day`), starting fresh each
pass from a uniform initial profile (the day's mean air temperature), to reach a
steady-periodic solution before the output is kept — see
[Configuring the model](../manual/configuring_the_model.md) and
[Time modes](../manual/diel_time_handling.md#Time-modes).

## Soil temperature by depth, across the year

```@example monthly
fig = Figure(size = (700, 400))
ax = Axis(fig[1, 1]; xlabel = "Hour of representative day (by month)", ylabel = "Soil temperature (°C)")
depth_labels = ["$(round(Int, ustrip(u"cm", d))) cm" for d in depths]
for (col, label) in zip((1, 4, 7, 10), depth_labels[[1, 4, 7, 10]])
    lines!(ax, ustrip.(u"°C", out.soil_temperature[:, col]); label)
end
axislegend(ax; position = :rt)
fig
```

The surface node tracks the diurnal cycle closely; deeper nodes are progressively
damped and lagged, as the [heat equation](../manual/soil_thermal_properties.md)
predicts.

## Substrate type changes the amplitude, not just the level

Swapping [`SoilProfile`](@ref) for a rock rather than a soil substrate lets more heat
penetrate deeper, reducing the amplitude near the surface while increasing it at
depth. The main lever for this is porosity, i.e. `bulk_density` relative to
`mineral_density`: soil's default (`bulk_density = 1.3 Mg/m^3`,
`mineral_density = 2.56 Mg/m^3`) is about half pore space, while solid rock has
almost none, so `bulk_density` needs to sit close to `mineral_density`. Rock's own
mineral conductivity and heat capacity are also higher than typical soil minerals,
and with almost no pore space, there's almost no room for soil moisture either:

```@example monthly
rock_profile = example_soil_profile(;
    bulk_density = 2.55u"Mg/m^3", mineral_density = 2.6u"Mg/m^3",
    mineral_conductivity = 2.5u"W/m/K", mineral_heat_capacity = 750.0u"J/kg/K",
)
rock_inputs = MicroInputs(;
    site, soil_profile = rock_profile,
    environment_minmax = example_monthly_weather(),
    environment_daily = example_daily_environment(),
    environment_hourly = example_hourly_environment(),
    initial_soil_temperature = fill(u"K"(7.74u"°C"), length(depths)),
    initial_soil_moisture = fill(0.01, length(depths)),
)
rock_out = solve(MicroProblem(model, rock_inputs))

fig2 = Figure(size = (700, 400))
ax2 = Axis(fig2[1, 1]; xlabel = "Hour (July, day 196)", ylabel = "Surface temperature (°C)")
july = (7 - 1) * 24 + 1:7 * 24
lines!(ax2, ustrip.(u"°C", out.soil_temperature[july, 1]); label = "soil")
lines!(ax2, ustrip.(u"°C", rock_out.soil_temperature[july, 1]); label = "rock")
axislegend(ax2)
fig2
```

## Next steps

- [Daily/hourly workflow](daily_hourly_workflow.md) — driving the model from real
  hourly weather-station data instead of monthly normals.
- [Configuring the model](../manual/configuring_the_model.md) — every other slot on
  `MicroModel` this tutorial left at its default.
