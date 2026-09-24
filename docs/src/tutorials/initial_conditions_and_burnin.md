# Initial conditions and spin-up

`initial_soil_temperature`/`initial_soil_moisture` on [`MicroInputs`](@ref) matter
differently depending on [time mode](../manual/diel_time_handling.md#Time-modes), and
getting this wrong silently produces results that look plausible but aren't what was
intended — this tutorial demonstrates both failure modes.

## Monthly mode resets every representative day

Under [`NonConsecutiveDayMode`](@ref) (the default — one independent representative
day per month), **every** day resets both soil temperature and soil moisture back to
their starting values before its spin-up iterations begin: temperature to
`initial_soil_temperature` if supplied, or otherwise to that day's own mean reference
air temperature; moisture to `initial_soil_moisture`, unconditionally. There is no
month-to-month memory — July's initial moisture is exactly what you passed in, not
whatever June's simulation produced.

```@example burnin
using Microclimate, Unitful, CairoMakie

n = length(Microclimate.DEFAULT_DEPTHS)
config = MicroConfig(; soil_moisture_strategy = DynamicSoilMoisture())

wet_start = solve(example_microclimate_problem(; config, initial_soil_moisture = fill(0.35, n)))
dry_start = solve(example_microclimate_problem(; config, initial_soil_moisture = fill(0.05, n)))

# December (the 12th representative day) is identical either way -- the very
# different starting moisture never had a chance to persist between months.
last_month = (12 - 1) * 24 + 1:12 * 24
maximum(abs.(wet_start.soil_moisture[last_month, :] .- dry_start.soil_moisture[last_month, :]))
```

That last value is (numerically) zero: January's moisture never survives to February's
spin-up. Monthly mode with `DynamicSoilMoisture` can't show a realistic wet-up/dry-down
trajectory — soil moisture's memory is weeks to months long, far longer than a
representative day's handful of spin-up iterations.

## Burn-in iterations within a single representative day

What `NonConsecutiveDayMode`'s `iterations_per_day` *does* control is how well each
day's own steady-periodic cycle is resolved. The *deepest* node isn't a useful probe
for this — it's a Dirichlet boundary condition, pinned to `deep_soil_temperature`
every hour regardless of iteration count — but the node just above it has real
thermal inertia and needs more than one pass to stop drifting from its (arbitrary)
reset value toward the day's real diurnal equilibrium:

```@example burnin
function deep_temperature_trace(iterations_per_day)
    problem = example_microclimate_problem()
    time_mode = NonConsecutiveDayMode(; iterations_per_day)
    out = solve(MicroProblem(problem.model, problem.inputs; problem.days, time_mode))
    out.soil_temperature[1:24, end - 1]  # second-deepest node (150 cm), first representative day
end

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "Hour", ylabel = "150 cm node temperature (°C)")
for iterations_per_day in (1, 3, 10, 30)
    lines!(ax, ustrip.(u"°C", deep_temperature_trace(iterations_per_day)); label = "$iterations_per_day iterations")
end
axislegend(ax)
fig
```

At 150 cm, 1 and 3 iterations are indistinguishable — the signal hasn't arrived yet. By
30, the per-iteration change is still *growing*, not shrinking: nowhere near converged.
`iterations_per_day = 3` (the default) isn't enough this deep, whatever the starting guess.

## The starting guess matters as much as the iteration count

Leaving `initial_soil_temperature` unset (`nothing`) resets each day to *its own* mean
air temperature instead of one fixed value for all twelve — much closer for shallow/mid
depths, where a single fixed value can be many degrees off for any given month:

```@example burnin
function temperature_trace(col, initial_soil_temperature)
    problem = example_microclimate_problem(; initial_soil_temperature)
    out = solve(MicroProblem(problem.model, problem.inputs; problem.days))  # default iterations_per_day = 3
    out.soil_temperature[1:24, col]
end

fig3 = Figure(size = (700, 500))
ax3a = Axis(fig3[1, 1]; ylabel = "Surface temperature (°C)", title = "0 cm")
ax3b = Axis(fig3[2, 1]; xlabel = "Hour", ylabel = "150 cm temperature (°C)", title = "150 cm")
for (init, label) in ((fill(u"K"(7.74u"°C"), n), "uniform annual mean"), (nothing, "day's own mean"))
    lines!(ax3a, ustrip.(u"°C", temperature_trace(1, init)); label)
    lines!(ax3b, ustrip.(u"°C", temperature_trace(n - 1, init)); label)
end
axislegend(ax3a)
fig3
```

`nothing` converges faster at the surface; the annual mean wins at 150 cm, since deep
soil tracks the annual cycle rather than any one day's mean. No uniform value suits
every depth — that needs a profile blending day's-mean-at-surface to annual-mean-at-depth,
which the package doesn't build automatically.

## Running daily instead, to let soil moisture actually evolve

[`MonthlyMinMaxEnvironment`](@ref) is built for independent representative
days — its per-day value series doesn't line up with a consecutive calendar run.
[`DailyMinMaxEnvironment`](@ref) is the consecutive-day analogue: one min/max pair
per actual calendar day, and it's automatically treated as consecutive (so day `i+1`'s
diel curve blends toward day `i+2`'s values instead of wrapping). Combined with
[`ConsecutiveDayMode`](@ref) — which does not reset moisture or temperature between
days — a run over many consecutive days lets `DynamicSoilMoisture` show a real
trajectory, including, with `spinup_first_day=true`, a proper spin-up on day 1 before
the consecutive run begins:

```@example burnin
days = collect(1:60)
site = example_site()
warming = range(0.0, 15.0; length = length(days))  # a slow synthetic spring warm-up
environment_minmax = DailyMinMaxEnvironment(; forcings = minmax_forcings(;
    reference_temperature_min = (268.15 .+ warming)u"K",
    reference_temperature_max = (278.15 .+ warming)u"K",
    reference_wind_speed_min = fill(0.5, length(days))u"m/s",
    reference_wind_speed_max = fill(3.0, length(days))u"m/s",
    reference_humidity_min = fill(0.4, length(days)),
    reference_humidity_max = fill(0.9, length(days)),
    cloud_cover_min = fill(0.2, length(days)),
    cloud_cover_max = fill(0.5, length(days)),
))

consecutive_inputs = MicroInputs(;
    site, soil_profile = example_soil_profile(),
    environment_minmax,
    environment_daily = example_daily_environment(days; rainfall = fill(0.0u"kg/m^2", length(days))),
    environment_hourly = example_hourly_environment(days),
    initial_soil_temperature = fill(u"K"(-2.0u"°C"), n),
    initial_soil_moisture = fill(0.35, n),
)
model = MicroModel(;
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model = example_soil_hydraulic_model(),
    config,
)
daily_out = solve(MicroProblem(model, consecutive_inputs; days,
    time_mode = ConsecutiveDayMode(; spinup_first_day = true)))

fig2 = Figure()
ax2 = Axis(fig2[1, 1]; xlabel = "Day", ylabel = "Surface soil moisture (m³/m³)")
lines!(ax2, daily_out.soil_moisture[1:24:end, 1])
fig2
```

Unlike the monthly comparison above, this trajectory depends on its starting value and
the rainfall/evapotranspiration history, not just the current day's forcing.

## Next steps

- [Configuring the model](../manual/configuring_the_model.md) for `MicroConfig`'s
  `convergence` field and why it doesn't control `NonConsecutiveDayMode`'s iteration
  count.
- [Daily/hourly workflow](daily_hourly_workflow.md) for driving a consecutive run from
  real weather data instead of interpolated monthly normals.
