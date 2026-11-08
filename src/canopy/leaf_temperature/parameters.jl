"""
    LeafParameters(; leaf_length=0.05u"m", leaf_width=0.02u"m", leaf_emissivity=0.97,
                      leaf_water_potential=0.0u"J/kg", canopy_projection_ratio=1.0)

Leaf traits.

- `leaf_length`, `leaf_width` — leaf size, for convection (`leaf_body`).
- `leaf_emissivity` — leaf longwave emissivity.
- `leaf_water_potential` — used when the caller gives none.
- `canopy_projection_ratio` — Campbell's leaf-angle ratio (1.0 = spherical, default;
  0.0 = vertical; `Inf` = horizontal), used by radiation and rain interception.
"""
@kwdef struct LeafParameters{LL,LW,LE,WP,X}
    leaf_length::LL = 0.05u"m"
    leaf_width::LW = 0.02u"m"
    leaf_emissivity::LE = 0.97
    leaf_water_potential::WP = 0.0u"J/kg"
    canopy_projection_ratio::X = 1.0
end
