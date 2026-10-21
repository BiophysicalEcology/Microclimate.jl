# Canopy

A single scalar `shade` fraction ([`NoCanopy`](@ref), the default) cannot represent
how radiation, wind, temperature and humidity change with height through a
canopy. [`MultilayerCanopy`](@ref) resolves the canopy into layers (one per 
`MicroModel.heights` entry at or below `canopy_height`) and solves shortwave, longwave, 
wind, air-profile, interception and leaf temperature sub-models per layer, each hour, 
converging the canopy against the soil surface temperature with [`PicardCanopyConvergence`](@ref).

The default wind model, [`ExponentialCanopyWindAttenuation`](@ref), requires at least
10 canopy layers (it replaces the bottom tenth of the profile with a local log-law
correction near the ground) — `MicroModel.heights` needs at least 10 entries at or
below `canopy_height` when using it, as in [Canopy workflow](../tutorials/canopy_workflow.md).

## Sources

The literature sources for each sub-model are as follows:

- **Shortwave**: [`TwoStreamRadiation`](@ref) is the classic Dickinson (1983)/
  Sellers (1985) two-stream scheme, with Campbell's (1990) ellipsoidal leaf-angle
  extinction coefficient.
- **Longwave**: three interchangeable schemes.
  [`LayeredRadiosityExchange`](@ref) (the `MultilayerCanopy` default) is an implicit
  reflecting two-stream-slab scheme following Flerchinger's **SHAW** model.
  [`LayeredLongwaveExchange`](@ref) is a from-scratch sequential, non-reflecting
  gap-fraction cascade. [`AllPairsLongwaveExchange`](@ref) is a direct **port** of
  micropoint's `twostreamvegCpp`/microclimlearn's `longwavebelow`.
- **Wind attenuation**: [`ExponentialCanopyWindAttenuation`](@ref) (exponential decay
  with cumulative plant area index from canopy top, after Nikolov and Zeller 2003) and
  [`MixingLengthCanopyWindAttenuation`](@ref) (a shelter-factor alternative).
- **In-canopy air profile**: [`RaupachLTheoryAirProfile`](@ref) is a Lagrangian 
  near-field/far-field scheme (its `:bulk` formulation option matches micropoint's 
  `LangrangianOne` closely); [`KTheoryAirProfile`](@ref) is a simpler K-theory diffusion
  alternative.
- **Rain interception**: [`LayeredRainInterception`](@ref) reuses micropoint's
  `rainintercept` drop-velocity formula (``3.78 \cdot rainfall^{0.067}``) in conjunction
  with the leaf-angle extinction machinery note above for shortwave radiation;
  [`VerticalRainInterception`](@ref) is a simpler bulk alternative;
  [`NoInterception`](@ref) switches it off (the `MultilayerCanopy` default).
- **Leaf temperature**: solved by
  [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) — not by
  Microclimate.jl itself. [`LinearizedLeafTemperature`](@ref) and
  [`RootFindLeafTemperature`](@ref) (`src/canopy/leaf_temperature/`) both call through
  to HeatExchange.jl's `convection`/`evaporation` (shape-specific Nusselt
  correlations, water-potential-driven surface humidity via the Kelvin equation,
  Campbell's (1990) leaf-angle-dependent aerodynamic width); see HeatExchange.jl's own
  docs for the underlying equations.

Related R resources:
[microclimc](https://github.com/ilyamaclean/microclimc)/
[microclimlearn tutorial](https://rpubs.com/ilyamaclean/microclimlearn), and TrenchR's
["Estimating microclimates"](https://cran.r-project.org/web/packages/TrenchR/vignettes/MicroclimateTutorial.html).

## Leaf traits and stomatal conductance

[`LeafParameters`](@ref) holds per-leaf structural/physiological traits shared across
sub-models — `leaf_length`/`leaf_width` (for HeatExchange.jl's boundary-layer
convection), `leaf_emissivity`, `canopy_projection_ratio` (Campbell's ellipsoidal
leaf-angle `x`, shared by the shortwave and rain-interception extinction
calculations), and `leaf_water_potential` (a fallback default, normally overridden
each hour once a moisture-coupled stomatal model is in use).

Stomatal conductance gates how much of a leaf's evaporative demand is actually met.
[`PrescribedStomatalConductance`](@ref) (the `MultilayerCanopy` default) is a simple
day/night switch — full conductance whenever the sun is up, cuticular-only at night —
with no coupling to soil moisture or photosynthesis.
[`MoistureResponsiveStomatalConductance`](@ref) instead closes the stomata smoothly as
leaf water potential falls, reusing `CampbellSoilHydraulics`'s own
`stomatal_closure_potential`/`stomatal_stability_parameter` (Campbell 1985) so the two
models can't drift out of sync — the same empirical closure relation used for ground
vegetation in [Soil moisture](soil_moisture.md#Stomatal-closure-and-actual-transpiration).
It only pairs with `CampbellSoilHydraulics` as the soil hydraulic model (a
`MethodError` otherwise, by construction), since it reads live leaf water potential
from that model's own solve.

## Rain interception

[`LayeredRainInterception`](@ref) tracks water storage on each canopy layer's leaves
separately: rain is extinguished layer by layer using the same leaf-angle machinery as
direct-beam shortwave (`ellipsoidal_extinction_coefficient`), at a rain-fall
angle `atan(wind_speed / raindrop_fall_velocity)` that assumes drops track the local
wind instantaneously. Each layer's storage capacity scales with its own plant area
index (`leaf_water_storage_capacity`, 0.1 kg/m² by default — roughly a 0.1 mm film);
water beyond capacity drips through to the next layer as ordinary incident rain.
[`VerticalRainInterception`](@ref) is a simpler bulk (whole-canopy, not per-layer)
alternative; [`NoInterception`](@ref) (the `MultilayerCanopy` default) switches
interception off entirely, so all rain reaches the ground unimpeded.

## Leaf dew and frost

[`MonteithLeafCondensation`](@ref) is the canopy's per-leaf analogue of
[ground dew/frost formation](evaporation_condensation.md#Dew-and-frost): a
single-term combination equation (net radiative loss weighted by the
saturation-vapour-pressure-curve slope, plus an aerodynamic term using the leaf's own
heat transfer coefficient directly) rather than the ground's two-term
energy/aerodynamic split — a thin leaf has negligible heat storage, so there's no
ground-heat-flux analogue term, matching Garratt and Segal's (1988) own omission of
that term for canopy dew.

## Assembling a canopy

```@example canopy
using Microclimate, Unitful

canopy_model = example_multilayer_canopy()
canopy_model.canopy_height, canopy_model.shortwave_model, canopy_model.longwave_model
```

[`plant_area_index_from_density`](@ref) builds a per-layer plant area index from a
vertical density profile, for canopies where leaf area isn't uniform with height. See
[Canopy workflow](../tutorials/canopy_workflow.md) for a full worked example, and
[Configuring the model](configuring_the_model.md) for how `canopy_model` slots into
`MicroModel` alongside the ground-surface physics.

## References

See the [References](references.md) page for the full bibliography.
