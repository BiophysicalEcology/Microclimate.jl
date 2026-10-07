# Terrain

Slope, aspect and hillshade change the radiation a surface receives, and therefore its
soil temperature — see [Radiation](../manual/radiation.md) for the equations. This
tutorial runs the same Madison, Wisconsin site used elsewhere in these docs
([`example_site`](@ref)) under a few different terrain configurations.

## Slope and aspect

The default is flat ground (`slope = 0°`). Here's the effect of a 45° south-facing
slope, compared to flat ground, on the 0%-shade simulation:

```@example terrain
using Microclimate, Unitful, CairoMakie

site = example_site()
soil_profile = example_soil_profile()
model = MicroModel(;
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model = example_soil_hydraulic_model(),
)
depths = model.depths

function run_at(site)
    inputs = MicroInputs(;
        site, soil_profile,
        environment_minmax = example_monthly_weather(),
        environment_daily = example_daily_environment(),
        environment_hourly = example_hourly_environment(),
        initial_soil_temperature = nothing,  # each representative day resets to its own mean air temperature
        initial_soil_moisture = fill(0.105, length(depths)),
    )
    solve(MicroProblem(model, inputs))
end

flat_out = run_at(site)

sloped_site = example_site(; slope = 45.0u"°", aspect = 180.0u"°")
sloped_out = run_at(sloped_site)
nothing # hide
```

```@example terrain
july = (7 - 1) * 24 + 1:7 * 24
december = (12 - 1) * 24 + 1:12 * 24

fig = Figure(size = (700, 600))
ax1 = Axis(fig[1, 1]; ylabel = "Surface temperature (°C)", title = "July (day 196)")
lines!(ax1, ustrip.(u"°C", flat_out.soil_temperature[july, 1]); label = "flat")
lines!(ax1, ustrip.(u"°C", sloped_out.soil_temperature[july, 1]); label = "45° south-facing")
axislegend(ax1)
ax2 = Axis(fig[2, 1]; xlabel = "Hour", ylabel = "Surface temperature (°C)", title = "December (day 349)")
lines!(ax2, ustrip.(u"°C", flat_out.soil_temperature[december, 1]); label = "flat")
lines!(ax2, ustrip.(u"°C", sloped_out.soil_temperature[december, 1]); label = "45° south-facing")
axislegend(ax2)
fig
```

At this latitude, with the sun low in the sky in winter, the slope's effect is
strongest on the December temperatures — a south-facing slope catches the low winter
sun far more directly than flat ground does.

## Hillshade

Gullies and gorges receive a shorter window of direct sun and have a smaller view of
the sky, so they don't cool down as much at night. `Site.horizon_angles` sets this: a
vector of horizon angles (the angle to the nearest obstruction, starting due north and
running clockwise at equal steps — 32 directions, 11.25° apart, matching
[Geomorphometry.jl](https://github.com/Deltares/Geomorphometry.jl)'s convention),
which also determine `sky_view_fraction` (the fraction of the sky not blocked by
terrain — see [Radiation](../manual/radiation.md)). Here's a steep-sided gully running
north-south — open along the north-south axis, walled in to the east and west:

```@example terrain
gully_horizon_angles = vcat(
    fill(0.0, 2), fill(65.0, 14),
    fill(0.0, 2), fill(65.0, 14),
)u"°"
gully_sky_view_fraction = 1 - sum(sin, gully_horizon_angles) / length(gully_horizon_angles)

gully_site = example_site(;
    horizon_angles = gully_horizon_angles, sky_view_fraction = gully_sky_view_fraction,
)
gully_out = run_at(gully_site)

fig2 = Figure()
ax = Axis(fig2[1, 1]; xlabel = "Hour (July)", ylabel = "Surface temperature (°C)")
lines!(ax, ustrip.(u"°C", flat_out.soil_temperature[july, 1]); label = "open")
lines!(ax, ustrip.(u"°C", gully_out.soil_temperature[july, 1]); label = "north-south gully")
axislegend(ax)
fig2
```

Most of the day the gully floor is shaded and stays well below the open-ground
temperature, but around solar noon the sun briefly aligns with the narrow north-south
opening and direct beam pours in — a sharp, short-lived spike rather than a gradual
warming, since the sun crosses from fully blocked to fully visible within about an
hour at this resolution.

## Next steps

- [Radiation](../manual/radiation.md) for the shortwave slope/aspect and longwave
  hillshade/view-factor equations behind both examples.
- [Configuring the model](../manual/configuring_the_model.md) — every other slot on
  `MicroModel` this tutorial left at its default.
