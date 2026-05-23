abstract type SoilPhaseTransitionScheme end

"""
    phase_transition!(scheme, buffers; temperatures, temperatures_past, accumulated_latent_heat, soil_moisture, depths)

Soil ice/water phase-change correction: detects 0°C crossings per layer,
accumulates latent heat against the layer's freezing budget, clamps the
temperature to 0°C while the budget is non-empty. Mutates `buffers` and
`accumulated_latent_heat` in place; returns the corrected temperature
SVector.
"""
function phase_transition! end

"""
    allocate_phase_transition(scheme, num_nodes)

Allocate per-layer scratch buffers for the freezing scheme.
"""
function allocate_phase_transition end
