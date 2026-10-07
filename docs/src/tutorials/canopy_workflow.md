# Canopy workflow

[`MultilayerCanopy`](@ref) resolves leaf temperature, in-canopy air temperature and
wind speed by height, rather than treating vegetation as a single scalar shade
fraction. See [Canopy](../manual/canopy.md) for the sub-models this composes and where
each one comes from.

## Building a canopy over a set of heights

A canopy needs enough entries in `MicroModel.heights` at or below `canopy_height` to
resolve into layers:

```@example canopy_workflow
using Microclimate, Unitful, CairoMakie

heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"   # 10 in-canopy layers up to 1 m, plus the 2 m reference height
canopy_model = example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 2.0)

canopy_out = solve(example_microclimate_problem(; heights, canopy_model))
bare_out = solve(example_microclimate_problem(; heights))   # NoCanopy, same everything else
size(canopy_out.canopy.leaf_temperature)
```

`canopy_out.canopy` (`CanopyOutput`) is zero-columned for a bare-ground
(`NoCanopy`) run and per-layer for `MultilayerCanopy` — the canopy changes the ground
surface temperature too, since it shades and slows the wind reaching the ground:

```@example canopy_workflow
canopy_out.soil_temperature[:, 1] != bare_out.soil_temperature[:, 1]
```

## A daytime canopy profile

```@example canopy_workflow
layer_heights = ustrip.(u"m", sort(heights[heights .<= canopy_model.canopy_height]; rev = true))
step = 13   # early afternoon on the representative day

fig = Figure(size = (700, 350))
ax1 = Axis(fig[1, 1]; xlabel = "Temperature (°C)", ylabel = "Height (m)", title = "Leaf and in-canopy air")
lines!(ax1, ustrip.(u"°C", canopy_out.canopy.leaf_temperature[step, :]), layer_heights; label = "leaf")
lines!(ax1, ustrip.(u"°C", canopy_out.canopy.air_temperature[step, :]), layer_heights; label = "in-canopy air")
axislegend(ax1)
ax2 = Axis(fig[1, 2]; xlabel = "Wind speed (m/s)", title = "In-canopy wind")
lines!(ax2, ustrip.(u"m/s", canopy_out.canopy.wind_speed[step, :]), layer_heights)
fig
```

Wind attenuates going down through the canopy ([`ExponentialCanopyWindAttenuation`](@ref)
or [`MixingLengthCanopyWindAttenuation`](@ref)), and leaf temperature tracks — but
doesn't equal — in-canopy air temperature, since it's solved from each leaf's own
energy balance via HeatExchange.jl (see [Canopy](../manual/canopy.md)).

## Soil moisture stress reaches the leaf

Pairing [`MoistureResponsiveStomatalConductance`](@ref) with
`DynamicSoilMoisture` lets Campbell's soil-water-supply calculation
(see [Soil moisture](../manual/soil_moisture.md)) constrain leaf water potential, and
so leaf temperature, through stomatal closure:

```@example canopy_workflow
config = MicroConfig(; soil_moisture_strategy = DynamicSoilMoisture())
n = length(Microclimate.DEFAULT_DEPTHS)

wet_out = solve(example_microclimate_problem(; heights, canopy_model, config,
    initial_soil_moisture = fill(0.35, n)))
dry_out = solve(example_microclimate_problem(; heights, canopy_model, config,
    initial_soil_moisture = fill(0.03, n)))

wet_out.canopy.leaf_temperature != dry_out.canopy.leaf_temperature
```

With [`PrescribedStomatalConductance`](@ref) (the `MultilayerCanopy` default) this
coupling is absent — stomata only gate day/night, ignoring `leaf_water_potential`
entirely.

## Next steps

- [Canopy](../manual/canopy.md) — which literature and reference implementations each
  sub-model draws on.
- [Configuring the model](../manual/configuring_the_model.md) — swapping any of the
  canopy's own sub-models the same way `MicroModel`'s own slots are swapped.
