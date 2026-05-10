"""
    MoninObukhov(; karman_constant=0.4, dyer_constant=16.0)

Monin–Obukhov similarity theory boundary-layer formulation. Holds the
empirical constants of the Φ_m relation. Previously these constants were
buried as fields on `MicroTerrain`, where they masqueraded as terrain
properties.

# References
Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971).
Flux–profile relationships in the atmospheric surface layer.
*Journal of the Atmospheric Sciences*, 28(2), 181–189.

Dyer, A. J. (1974). A review of flux–profile relationships.
*Boundary-Layer Meteorology*, 7(3), 363–372.
"""
@kwdef struct MoninObukhov{KC,DC} <: AbstractBoundaryLayerModel
    karman_constant::KC = 0.4
    dyer_constant::DC = 16.0
end
