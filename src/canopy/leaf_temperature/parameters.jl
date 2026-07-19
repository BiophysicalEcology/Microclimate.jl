"""
    LeafParameters(; leaf_length=0.05u"m", leaf_width=0.02u"m", leaf_emissivity=0.97,
                      leaf_water_potential=0.0u"J/kg", leaf_angle_distribution_parameter=1.0)

Per-leaf structural/physiological data — leaf traits, not a swappable
model choice (analogous to `SoilProfile`'s per-depth soil data).

- `leaf_length`, `leaf_width` — for `HeatExchange.jl`'s boundary-layer
  convection model (`leaf_body`)
- `leaf_emissivity` — leaf longwave emissivity
- `leaf_water_potential` — fallback default for standalone/no-soil-coupling
  use; normally overridden each hour by `CampbellSoilHydraulics`'s live
  leaf-water-potential solve
- `leaf_angle_distribution_parameter` — Campbell's ellipsoidal leaf-angle `x`
  (1.0 = spherical, default; 0.0 = horizontal; `Inf` = vertical). Lives here
  rather than on one radiative-transfer model since both
  [`TwoStreamRadiation`](@ref) and the rain-interception model use it.
"""
@kwdef struct LeafParameters{LL,LW,LE,WP,X}
    leaf_length::LL = 0.05u"m"
    leaf_width::LW = 0.02u"m"
    leaf_emissivity::LE = 0.97
    leaf_water_potential::WP = 0.0u"J/kg"
    leaf_angle_distribution_parameter::X = 1.0
end
