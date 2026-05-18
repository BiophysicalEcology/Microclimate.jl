abstract type AbstractApparentHeatCapacity end

"""
    apparent_heat_capacity(method, T, cp_base, L_fusion)

Return the effective specific heat at temperature `T` for an ice/water phase-change
region. `cp_base` is the snow specific heat away from the phase transition;
`L_fusion` is the latent heat of fusion. The augmentation models energy absorbed
during melting without tracking phase state explicitly.
"""
function apparent_heat_capacity end
