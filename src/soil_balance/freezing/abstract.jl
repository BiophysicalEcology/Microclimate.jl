abstract type SoilPhaseTransitionModel end

"""
    phase_transition!(model, buffers; temperatures, temperatures_past, accumulated_latent_heat, soil_moisture, depths)

Soil ice/water phase-change correction: detects 0°C crossings per layer,
accumulates latent heat against the layer's freezing budget, clamps the
temperature to 0°C while the budget is non-empty. Mutates `buffers` and
`accumulated_latent_heat` in place; returns the corrected temperature
SVector.
"""
function phase_transition! end

"""
    allocate_phase_transition(model, num_nodes)

Allocate per-layer scratch buffers for the freezing model.
"""
function allocate_phase_transition end

"""
    frozen_water_content!(model, buffers, accumulated_latent_heat, soil_moisture)

Per-layer frozen water content (m³/m³): `accumulated_latent_heat / (mass × L_fusion)`
times `soil_moisture`. Mutates and returns `buffers`' own scratch vector.
"""
function frozen_water_content! end
