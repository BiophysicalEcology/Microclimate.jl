# Snow demo

This demonstrates [`SnowModel`](@ref) buffering soil temperature through a winter
onset — the mechanism behind Kearney (2020)'s finding that shallow soil warms at
roughly a third the rate of air temperature at snow-affected sites (see
[Snow](../manual/snow.md#Background)).

## A synthetic cold winter

```@example snow_demo
using Microclimate, Unitful, CairoMakie

site = example_site()
depths = Microclimate.DEFAULT_DEPTHS
n = length(depths)
days = collect(335:15:365)   # late Nov through Dec, one representative day per fortnight

reference_temperature_min = (range(-2.0, -18.0; length = length(days)) .+ 273.15)u"K"
reference_temperature_max = (range(4.0, -8.0; length = length(days)) .+ 273.15)u"K"
environment_minmax = MonthlyMinMaxEnvironment(; forcings = minmax_forcings(;
    reference_temperature_min, reference_temperature_max,
    reference_wind_speed_min = fill(0.5, length(days))u"m/s",
    reference_wind_speed_max = fill(4.0, length(days))u"m/s",
    reference_humidity_min = fill(0.6, length(days)),
    reference_humidity_max = fill(0.95, length(days)),
    cloud_cover_min = fill(0.3, length(days)),
    cloud_cover_max = fill(0.6, length(days)),
))
environment_daily = example_daily_environment(days; rainfall = fill(15.0, length(days))u"kg/m^2")
environment_hourly = example_hourly_environment(days; elevation = site.elevation)

snow_model = SnowModel(; snow_temperature_threshold = 1.0u"°C")

function run(; snow_model)
    model = MicroModel(;
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model = example_soil_hydraulic_model(),
        snow_model, depths,
    )
    inputs = MicroInputs(;
        site, soil_profile = example_soil_profile(depths),
        environment_minmax, environment_daily, environment_hourly,
        initial_soil_temperature = fill(u"K"(2.0u"°C"), n),
        initial_soil_moisture = fill(0.2, n),
    )
    solve(MicroProblem(model, inputs; days, time_mode = ConsecutiveDayMode()))
end

with_snow = run(; snow_model)
without_snow = run(; snow_model = NoSnow())
nothing # hide
```

## Snow depth and its effect on soil temperature

```@example snow_demo
fig = Figure(size = (700, 500))
ax1 = Axis(fig[1, 1]; ylabel = "Snow depth (cm)")
lines!(ax1, ustrip.(u"cm", with_snow.snow_depth))
ax2 = Axis(fig[2, 1]; xlabel = "Hour", ylabel = "Surface soil temperature (°C)")
lines!(ax2, ustrip.(u"°C", with_snow.soil_temperature[:, 1]); label = "with snow")
lines!(ax2, ustrip.(u"°C", without_snow.soil_temperature[:, 1]); label = "without snow (rain instead)")
axislegend(ax2; position = :rb)
fig
```

As the pack builds, surface soil temperature decouples from the falling air
temperature and settles close to 0°C, buffered by the latent heat of fusion at the
snow-soil interface — while the no-snow run (precipitation still falls, as rain)
continues to track the cooling air directly.

## Next steps

- [Snow](../manual/snow.md) for the density, albedo, conductivity and melt equations
  behind this behaviour.
- [Soil thermal properties](../manual/soil_thermal_properties.md) for the underlying
  heat equation and phase-transition treatment near 0°C.
