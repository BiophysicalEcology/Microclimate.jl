"""
    AbstractCanopyInterceptionModel

How a canopy resolves rain interception: per-layer catch vs. through-fall,
and how a wetted leaf surface enhances evaporation. Held on
`MultilayerCanopy.interception_model`.

- [`NoInterception`](@ref) (default): rain passes straight through,
  matching `NoCanopy`/`NoSnow`'s off-by-default convention.
- [`LayeredRainInterception`](@ref): per-layer extinction/storage/drip.
"""
abstract type AbstractCanopyInterceptionModel end

"""
    allocate_interception(model, canopy_height, plant_area_index, heights, n_layers, boundary_layer_model)

Pre-allocate and precompute the structural buffers a model needs.
"""
function allocate_interception end

"""
    canopy_interception!(buffers, model, canopy_projection_ratio; rainfall, wind_speed)

Cascade `rainfall` (kg/m², arriving at canopy top) down through the layers,
updating each layer's stored leaf-surface water in place, and return
`(; ground_throughfall)` — the rain reaching the ground this hour.
"""
function canopy_interception! end

"""
    deplete_canopy_water!(model, buffers, layer, evaporated_mass)

Subtract `evaporated_mass` (kg/m², ground area) from layer `layer`'s stored
leaf-surface water, clamped at zero. Called once per hour after Picard
convergence (not itself Picard-iterated). No-op for `NoInterception`.
"""
function deplete_canopy_water! end

"""
    reset_interception!(model, buffers)

Zero accumulated leaf-surface water. Day-boundary reset (like
`pool`/`soil_moisture`), not per-Picard-iteration. No-op for `NoInterception`.
"""
function reset_interception! end

"""
    wet_canopy_fraction(model, buffers, layer) -> Float64 ∈ [0, 1]

Fractional saturation of layer `layer`'s stored leaf-surface water (0 =
dry, 1 = at capacity), for [`blend_stomatal_conductance`](@ref). Always 0
for `NoInterception`.
"""
function wet_canopy_fraction end

"""
    available_canopy_water(model, buffers, layer) -> kg/m^2 (ground area)

Leaf-surface water currently stored in layer `layer`, the hard ceiling on
how much [`deplete_canopy_water!`](@ref) can remove this hour. Used to cap
the wet-surface evaporation flux so it can't claim more water than is
physically stored (see `canopy_energy_balance!`'s depletion loop). Always 0
for `NoInterception`.
"""
function available_canopy_water end

"""
    snapshot_interception!(model, buffers)

Save day-start leaf-surface water into `buffers`' own scratch copy. No-op
for `NoInterception`.
"""
function snapshot_interception! end

"""
    restore_interception!(model, buffers)

Restore leaf-surface water from the last [`snapshot_interception!`](@ref).
No-op for `NoInterception`.
"""
function restore_interception! end

"""
    WET_SURFACE_CONDUCTANCE

A stomatal/cuticular conductance large enough that `HeatExchange.evaporation`'s
series combination with the boundary-layer mass-transfer coefficient
collapses to the boundary layer alone (a free water film has no
stomatal/cuticular resistance). Verified at ~1 m/s in-canopy wind: a further
10x increase changes evaporation by <0.1% — see
`test/canopy_interception_test.jl`.
"""
const WET_SURFACE_CONDUCTANCE = 1000.0u"mol/m^2/s"

"""
    blend_stomatal_conductance(dry::LeafEvaporationParameters, wet_fraction)

Blend `dry` toward [`WET_SURFACE_CONDUCTANCE`](@ref) as `wet_fraction` (from
[`wet_canopy_fraction`](@ref)) rises 0→1. Cuticular conductance is
unaffected — already surface-independent.
"""
function blend_stomatal_conductance(dry::LeafEvaporationParameters, wet_fraction)
    wet_fraction <= 0.0 && return dry
    return LeafEvaporationParameters(;
        abaxial_vapour_conductance = dry.abaxial_vapour_conductance + wet_fraction * (WET_SURFACE_CONDUCTANCE - dry.abaxial_vapour_conductance),
        adaxial_vapour_conductance = dry.adaxial_vapour_conductance + wet_fraction * (WET_SURFACE_CONDUCTANCE - dry.adaxial_vapour_conductance),
        cuticular_conductance = dry.cuticular_conductance,
    )
end
