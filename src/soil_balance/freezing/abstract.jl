"""
    SoilPhaseTransitionModel

Supertype for soil ice/water phase-transition correction models:
[`PhaseTransitionLatentHeat`](@ref) and [`IsothermalPhaseTransition`](@ref).
"""
abstract type SoilPhaseTransitionModel end

"""
    phase_transition!(model, buffers; temperatures, temperatures_past, accumulated_latent_heat, soil_moisture, depths)

Soil ice/water phase-change correction: holds layers at 0°C as `model` says,
accumulating latent heat against each layer's freezing budget (water mass × latent
heat of fusion). Mutates `buffers` and `accumulated_latent_heat` in place; returns the
corrected temperature SVector.
"""
function phase_transition! end

"""
    allocate_phase_transition(model, num_nodes)

Allocate per-layer scratch buffers for the freezing model.
"""
function allocate_phase_transition(::SoilPhaseTransitionModel, num_nodes::Int)
    layer_mass = zeros(Float64, num_nodes)u"kg"
    temperature_scratch = MVector{num_nodes, typeof(0.0u"K")}(undef)
    frozen_water_content = zeros(Float64, num_nodes)
    return (; layer_mass, temperature_scratch, frozen_water_content)
end

# Each layer runs from its node to the node below, and the deepest to 1 m below the deepest
# node. A layer's temperature is the mean of its two nodes, and `layer_phase_transition!`
# corrects those two nodes.
function phase_transition!(model::SoilPhaseTransitionModel, buffers::NamedTuple;
    temperatures::AbstractVector,
    temperatures_past::AbstractVector,
    accumulated_latent_heat::AbstractVector,
    soil_moisture::AbstractVector,
    depths::AbstractVector,
)
    (; layer_mass, temperature_scratch) = buffers
    temperature = temperature_scratch
    temperature .= temperatures
    for i in eachindex(depths)
        layer_mass[i] = layer_water_mass(depths, soil_moisture, i)
        if layer_mass[i] <= 0.0u"kg"
            accumulated_latent_heat[i] = 0.0u"J"
            continue
        end
        layer_phase_transition!(model, temperature, accumulated_latent_heat, i;
            layer_mass=layer_mass[i],
            mean_temperature=layer_mean_temperature(temperature, i),
            mean_temperature_past=layer_mean_temperature(temperatures_past, i))
    end
    return SVector(temperature)
end

"""
    layer_phase_transition!(model, temperature, accumulated_latent_heat, i;
        layer_mass, mean_temperature, mean_temperature_past)

Phase change in layer `i` with water mass `layer_mass`, now at `mean_temperature` and at
`mean_temperature_past` an hour ago. Updates `accumulated_latent_heat[i]` and the layer's
nodes in `temperature`.
"""
function layer_phase_transition! end

function layer_water_mass(depths, soil_moisture, i)
    soil_moisture[i] > 0.0 || return 0.0u"kg"
    thickness = i < length(depths) ? depths[i+1] - depths[i] : 100.0u"cm"
    return u"m"(thickness) * 1000.0u"kg/m" * soil_moisture[i]
end

layer_mean_temperature(temperature, i) =
    i < length(temperature) ? 0.5 * (temperature[i] + temperature[i+1]) : temperature[i]

# Applies `f` to the temperatures of layer `i`'s nodes
function update_layer_nodes!(f, temperature, i)
    temperature[i] = f(temperature[i])
    if i < length(temperature)
        temperature[i+1] = f(temperature[i+1])
    end
    return temperature
end

"""
    frozen_water_content!(model, buffers, accumulated_latent_heat, soil_moisture)

Per-layer frozen water content (m³/m³): `accumulated_latent_heat / (mass × L_fusion)`
times `soil_moisture`. Mutates and returns `buffers`' own scratch vector.
"""
function frozen_water_content! end

# Fixed for the whole hydrology sub-step loop -- not re-derived from
# water_content_new mid-solve, which would freeze/thaw water with no
# matching latent-heat update.
function frozen_water_content!(::SoilPhaseTransitionModel, buffers::NamedTuple, accumulated_latent_heat::AbstractVector, soil_moisture::AbstractVector)
    (; layer_mass, frozen_water_content) = buffers
    for i in eachindex(frozen_water_content)
        max_latent_heat = LATENT_HEAT_FUSION * layer_mass[i]
        fraction = max_latent_heat > zero(max_latent_heat) ?
            clamp(accumulated_latent_heat[i] / max_latent_heat, 0.0, 1.0) : 0.0
        frozen_water_content[i] = fraction * soil_moisture[i]
    end
    return frozen_water_content
end

