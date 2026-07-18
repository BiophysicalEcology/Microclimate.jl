abstract type AbstractCanopyModel end

"""
    allocate_canopy(model, heights)

Pre-allocate the per-layer buffers a canopy model needs to evaluate its
physics each hour, sized to the sub-canopy portion of `heights` (the same
`MicroModel.heights` vector the boundary-layer model already uses for its
above-canopy profile). Structural quantities that depend only on `model`
(the two-stream radiative-transfer optics, view weights, and similar
canopy-structure constants) are computed once here, not per hour.
"""
function allocate_canopy end

"""
    canopy_radiation!(buffers, model; zenith_angle, direct_horizontal_irradiance,
                       diffuse_horizontal_irradiance, ground_reflectance)

Compute layer-resolved absorbed shortwave radiation through the canopy for
the current hour, writing per-layer results into `buffers` and returning
`(; ground_absorbed_shortwave, canopy_absorbed_shortwave)`. `ground_reflectance`
is supplied by the caller (e.g. from `Site.albedo` or the current
`environment_instant`) rather than stored on the canopy model, matching how
other radiation models here read albedo/emissivity from their caller.
"""
function canopy_radiation! end
