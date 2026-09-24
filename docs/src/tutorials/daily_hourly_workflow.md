# Daily/hourly workflow

[`HourlyTimeseries`](@ref) drives the model directly from measured hourly weather —
air temperature, humidity, wind speed, cloud cover, solar radiation, rainfall —
instead of interpolating from monthly or daily min/max values. Kearney and Maino
(2018)'s continental-scale tests use the same approach at daily resolution with
gridded weather data (see [Soil moisture](../manual/soil_moisture.md#Background)).

Passing real hourly data only takes effect when `MicroInputs.environment_minmax` is
`nothing` — if a monthly/daily min-max forcing is also supplied, it takes priority for
`reference_temperature`/`reference_humidity`/`reference_wind_speed`/`cloud_cover`
regardless of what's in `environment_hourly` (see [`example_hourly_environment`](@ref)).

## A synthetic two-day hourly record

In place of a real weather-station download, this builds a plausible two-day hourly
record by hand — a simple diel cycle for each variable — to stand in for observed
data:

```@example daily_hourly
using Microclimate, Unitful, FluidProperties, CairoMakie

site = example_site()
depths = Microclimate.DEFAULT_DEPTHS
n = length(depths)
days = collect(180:181)   # two days, late June
hours_per_day = collect(0.0:1.0:23.0)
hours = repeat(hours_per_day, length(days))

diel(hour, minval, maxval) = minval + (maxval - minval) * max(0.0, sin((hour - 6) / 12 * π))

reference_temperature = [diel(h, 12.0, 28.0) for h in hours]u"°C"
reference_humidity = [1.0 - 0.5 * max(0.0, sin((h - 6) / 12 * π)) for h in hours]
reference_wind_speed = fill(2.0, length(hours))u"m/s"
observed_cloud_cover = fill(0.3, length(hours))   # the real cloud cover, if it were known
nothing # hide
```

Two back-to-back daily cycles for each variable feeding the run:

```@example daily_hourly
fig0 = Figure(size = (700, 500))
ax0a = Axis(fig0[1, 1]; ylabel = "Air temperature (°C)")
lines!(ax0a, ustrip.(u"°C", reference_temperature))
ax0b = Axis(fig0[2, 1]; ylabel = "Relative humidity")
lines!(ax0b, reference_humidity)
ax0c = Axis(fig0[3, 1]; xlabel = "Hour (two days)", ylabel = "Wind speed (m/s)")
lines!(ax0c, ustrip.(u"m/s", reference_wind_speed))
fig0
```

## Estimating cloud cover from solar radiation

Weather stations rarely measure cloud cover directly, but it's needed for the
longwave sky-temperature calculation (see [Radiation](../manual/radiation.md)).
Kearney and Porter (2017) describe the technique used to fill this gap: run a
preliminary clear-sky simulation for the site, then use the ratio of observed to
clear-sky solar radiation as an estimate of the proportion of cloud cover for each
hour:

```@example daily_hourly
clear_sky_hourly = HourlyTimeseries(;
    pressure = fill(atmospheric_pressure(site.elevation), length(hours)),
    reference_temperature, reference_humidity, reference_wind_speed,
    global_radiation = nothing, longwave_radiation = nothing,
    cloud_cover = zeros(length(hours)), rainfall = nothing, zenith_angle = nothing,
)
clear_sky_problem = MicroProblem(
    MicroModel(; soil_properties_model = example_soil_properties_model(),
                 soil_hydraulic_model = example_soil_hydraulic_model(), depths),
    MicroInputs(; site, soil_profile = example_soil_profile(depths),
        environment_minmax = nothing, environment_daily = example_daily_environment(days; rainfall = zeros(length(days))u"kg/m^2"),
        environment_hourly = clear_sky_hourly,
        initial_soil_temperature = fill(u"K"(20.0u"°C"), n), initial_soil_moisture = fill(0.15, n)),
    ; days, time_mode = ConsecutiveDayMode(),
)
clear_sky_radiation = solve(clear_sky_problem).global_radiation

# stand-in for a real observed record: clear-sky scaled down by observed_cloud_cover
observed_radiation = clear_sky_radiation .* (1 .- 0.7 .* observed_cloud_cover)
estimated_cloud_cover = clamp.(1 .- observed_radiation ./ max.(clear_sky_radiation, 1e-6u"W/m^2"), 0.0, 1.0)
nothing # hide
```

```@example daily_hourly
fig1 = Figure(size = (700, 350))
ax1a = Axis(fig1[1, 1]; ylabel = "Radiation (W/m²)")
lines!(ax1a, ustrip.(u"W/m^2", clear_sky_radiation); label = "clear-sky")
lines!(ax1a, ustrip.(u"W/m^2", observed_radiation); label = "observed")
axislegend(ax1a)
ax1b = Axis(fig1[2, 1]; xlabel = "Hour (two days)", ylabel = "Cloud cover")
lines!(ax1b, estimated_cloud_cover; label = "estimated")
lines!(ax1b, observed_cloud_cover; label = "true")
axislegend(ax1b)
fig1
```

## Driving the model from hourly data

With cloud cover (real, or estimated as above) in hand, build the real
[`HourlyTimeseries`](@ref) and solve, again with `environment_minmax = nothing` so the
hourly record drives the run:

```@example daily_hourly
observed_hourly = HourlyTimeseries(;
    pressure = fill(atmospheric_pressure(site.elevation), length(hours)),
    reference_temperature, reference_humidity, reference_wind_speed,
    global_radiation = nothing, longwave_radiation = nothing,
    cloud_cover = estimated_cloud_cover, rainfall = nothing, zenith_angle = nothing,
)
problem = MicroProblem(
    MicroModel(; soil_properties_model = example_soil_properties_model(),
                 soil_hydraulic_model = example_soil_hydraulic_model(), depths),
    MicroInputs(; site, soil_profile = example_soil_profile(depths),
        environment_minmax = nothing, environment_daily = example_daily_environment(days; rainfall = zeros(length(days))u"kg/m^2"),
        environment_hourly = observed_hourly,
        initial_soil_temperature = fill(u"K"(20.0u"°C"), n), initial_soil_moisture = fill(0.15, n)),
    ; days, time_mode = ConsecutiveDayMode(),
)
out = solve(problem)
size(out.soil_temperature)
```

```@example daily_hourly
fig2 = Figure()
ax2 = Axis(fig2[1, 1]; xlabel = "Hour (two days)", ylabel = "Temperature (°C)")
lines!(ax2, ustrip.(u"°C", reference_temperature); label = "air (reference height)")
lines!(ax2, ustrip.(u"°C", out.soil_temperature[:, 1]); label = "soil surface")
lines!(ax2, ustrip.(u"°C", out.soil_temperature[:, end - 1]); label = "150 cm depth")
axislegend(ax2)
fig2
```

The surface node tracks the air's diel cycle closely, damped and lagged; the 150 cm
node barely moves across just two days — consistent with [Soil thermal
properties](../manual/soil_thermal_properties.md)'s point that deep nodes respond on a
much longer timescale than a single day/night cycle.

[`ConsecutiveDayMode`](@ref) is the right [time mode](../manual/diel_time_handling.md#Time-modes)
for a real, continuous hourly record — each day's initial soil profile carries over
from the previous day's last hour, rather than restarting from a uniform profile as
[`NonConsecutiveDayMode`](@ref) does for independent representative days.

## Next steps

- [Monthly workflow](monthly_workflow.md) — the same model driven by monthly climate
  normals instead.
- [Soil moisture](../manual/soil_moisture.md) for how well this kind of run validates
  against observed soil moisture at continental scale.
