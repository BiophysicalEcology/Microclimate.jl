"""
    LeafParameters(; leaf_length=0.05u"m", leaf_width=0.02u"m", leaf_emissivity=0.97,
                      leaf_water_potential=0.0u"J/kg", leaf_angle_distribution_parameter=1.0)

Per-leaf structural/physiological data — not a swappable model choice, just
the leaf's own traits (analogous to `SoilProfile` holding per-depth soil
data rather than a soil-process choice).

- `leaf_length`, `leaf_width` — leaf dimensions, for the `HeatExchange.jl`
  boundary-layer convection model (`leaf_body`)
- `leaf_emissivity` — leaf longwave emissivity
- `leaf_water_potential` — a *fallback default* for standalone/no-soil-coupling
  use. The natural, intended source is dynamic: `CampbellSoilHydraulics`'s
  own leaf-water-potential solve (`soil_hydraulics/campbell.jl`), read fresh
  each hour by whichever orchestration loop calls `leaf_temperature` — this
  field is not that loop, just what's used when nothing supplies a live value.
- `leaf_angle_distribution_parameter` — Campbell's ellipsoidal leaf-angle
  distribution parameter `x` (1.0 = spherical, the default; 0.0 = horizontal
  leaves; `Inf` = vertical leaves). Lives here rather than on any one
  radiative-transfer model because it's a geometric leaf trait more than one
  sub-model reads: [`TwoStreamRadiation`](@ref) (shortwave extinction) and
  the rain-interception model (rain-angle extinction) both use it.
"""
@kwdef struct LeafParameters{LL,LW,LE,WP,X}
    leaf_length::LL = 0.05u"m"
    leaf_width::LW = 0.02u"m"
    leaf_emissivity::LE = 0.97
    leaf_water_potential::WP = 0.0u"J/kg"
    leaf_angle_distribution_parameter::X = 1.0
end
