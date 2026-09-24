# Phase transition and frozen soil

Water freezing and thawing is handled by two separate mechanisms in two separate
parts of the model — an hourly correction for the *soil* column, and a continuous
one inside the *snow* column's own heat ODE — plus a downstream effect on soil
hydraulics: once a soil layer is frozen, water is treated as unable to flow through
it, either through the soil matrix or through roots.

## The hourly latent-heat correction (soil)

[`PhaseTransitionLatentHeat`](@ref) is soil's only freeze/thaw mechanism: the
[soil heat ODE](soil_thermal_properties.md) and its
[`CampbelldeVriesSoilProperties`](@ref) thermal properties have no phase-transition
term of their own. Instead, `PhaseTransitionLatentHeat` applies a discrete correction
once per hour, after the ODE has advanced: it detects whether
the mean temperature of each pair of adjacent nodes has crossed 0°C between the
previous and current hour, and if so, accumulates (on freezing) or depletes (on
thawing) a per-layer latent-heat budget against `mass × L_fusion`. While that budget
is non-empty, the layer's temperature is clamped to exactly 0°C; once it's exhausted,
the clamp releases and the layer cools or warms freely again. This closely follows
NicheMapR Fortran's `OSUB.f:571-602`, and fires only on an actual crossing — a layer
that's already clamped to 0°C stays there without re-triggering the branch.

`frozen_water_content!` converts each layer's accumulated latent-heat budget
into a frozen fraction of that layer's `soil_moisture` (`accumulated_latent_heat /
(mass × L_fusion)`, clamped to `[0, 1]`) — the quantity the soil hydraulics below
actually consumes.

## Apparent heat capacity, inside the ODE (snow)

Snow, unlike soil, is solved with an *apparent* specific heat capacity, boosted near
0°C to absorb the latent heat of fusion as [`SnowModel`](@ref)'s own heat ODE
integrates through the freezing point, rather than treating it as a separate hourly
event. This is a distinct mechanism from soil's hourly correction above, plugged in
via `SnowModel`'s `apparent_heat_capacity` field (default [`BonacinaStep`](@ref)).
Four interchangeable forms are available:

- [`BonacinaStep`](@ref) — a step function over a fixed band (-0.45°C to 0.4°C by
  default), reproducing NicheMapR Fortran's `SOILPROPS.f` exactly. Its discontinuity
  is hard on adaptive ODE solvers: `Tsit5` tends to pin the temperature at the band
  edge rather than cross it.
- [`TanhSmoothed`](@ref) — the same step, replaced by a pair of `tanh` ramps at the
  band edges, so the solver sees a C∞-smooth `cp` and can cross the band in normal
  steps, while matching `BonacinaStep`'s integrated latent budget closely.
- [`Gaussian`](@ref) — a normal-distribution bump centred at the melting point, whose
  integral equals the latent heat of fusion exactly.
- [`WestermannSigmoid`](@ref) — derives `cp` from the derivative of an assumed
  unfrozen-water-fraction sigmoid (the CryoGrid formulation; Westermann et al. 2011,
  2016), so the shape emerges from a physical liquid-fraction model rather than being
  imposed directly.

## Frozen soil blocks water transport

Once a layer has a nonzero frozen water content, two independent things stop water
moving through it:

- **Through the soil matrix**: `ice_impeded_conductivity` scales hydraulic
  conductivity downward following Bloomsburg and Wang's (1969) relation, based on the
  *ice-free* porosity remaining (`porosity - ice_content`). Above a threshold
  (`ICE_IMPEDANCE_MIN_POROSITY = 0.13` m³/m³) conductivity is reduced roughly in
  proportion to how much ice-free pore space is left; at or below it, conductivity is
  floored at a tiny fraction of its unimpeded value (`ICE_CONDUCTIVITY_FLOOR_FACTOR =
  1e-6`, to avoid a literal division/root-finding singularity downstream) bit in practice 
  this results in no meaningful flow.
- **Through roots**: root resistance normally *dwarfs* soil resistance so the the mechanism
  above preventing frozen water moving through soil does not prevent it moving through roots.
  Thus, a layer whose ice-free porosity has crossed the same `ICE_IMPEDANCE_MIN_POROSITY` 
  threshold is additionally treated as **rootless**: its `root_resistance` is set prohibitivley
  high irrespective of that node's `root_density`. See [Soil moisture](soil_moisture.md#Initialisation) 
  for how `root_resistance` normally factors into water uptake.

## References

See the [References](references.md) page for the full bibliography.
