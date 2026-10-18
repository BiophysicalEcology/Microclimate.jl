"""
    BloomsburgWangFrozenConductivity(; min_ice_free_porosity=0.13, conductivity_floor=1e-6)

Hydraulic conductivity falls linearly with ice-free porosity and is blocked
once that drops to `min_ice_free_porosity` (Bloomsburg & Wang 1969, as used in
SHAW's `SOILHK`). `conductivity_floor` is the factor kept when blocked, so
conductivity stays nonzero and root uptake resistance stays finite.

# References
Bloomsburg, G. L. & Wang, S. J. (1969). Effect of moisture content on
permeability of frozen soils. Soil Sci. Soc. Am. J. 33:686-691.
"""
@kwdef struct BloomsburgWangFrozenConductivity{P,F} <: AbstractFrozenConductivity
    min_ice_free_porosity::P = 0.13
    conductivity_floor::F = 1e-6
end

function conductivity_factor(model::BloomsburgWangFrozenConductivity, ice_content, porosity)
    ice_content <= zero(ice_content) && return 1.0
    (; min_ice_free_porosity, conductivity_floor) = model
    ice_free_porosity = porosity - ice_content
    ice_free_porosity <= min_ice_free_porosity && return conductivity_floor
    return max(conductivity_floor, (ice_free_porosity - min_ice_free_porosity) / (porosity - min_ice_free_porosity))
end

conductivity_blocked(model::BloomsburgWangFrozenConductivity, ice_content, porosity) =
    porosity - ice_content <= model.min_ice_free_porosity
