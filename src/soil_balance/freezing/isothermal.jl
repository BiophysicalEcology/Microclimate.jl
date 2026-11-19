"""
    IsothermalPhaseTransition()

Soil ice/water phase change that holds each layer at 0°C while it freezes or thaws. Heat
lost below 0°C freezes the layer's water, and heat gained above 0°C melts its ice, until
the layer's latent heat budget (water mass × latent heat of fusion) is used up. Only the
remainder changes the layer's temperature. As in [`PhaseTransitionLatentHeat`](@ref), a
layer's heat per kelvin is that of its water.
"""
struct IsothermalPhaseTransition <: SoilPhaseTransitionModel end

function layer_phase_transition!(::IsothermalPhaseTransition, temperature, accumulated_latent_heat, i;
    layer_mass, mean_temperature, mean_temperature_past,
)
    freezing_point = u"K"(0.0u"°C")
    max_latent_heat = LATENT_HEAT_FUSION * layer_mass
    mean_temperature < freezing_point && accumulated_latent_heat[i] >= max_latent_heat && return temperature
    mean_temperature > freezing_point && accumulated_latent_heat[i] <= 0.0u"J" && return temperature
    heat_per_kelvin = layer_mass * water_heat_capacity
    # Heat below 0°C freezes water, heat above it melts ice
    latent_heat = accumulated_latent_heat[i] + (freezing_point - mean_temperature) * heat_per_kelvin
    accumulated_latent_heat[i] = clamp(latent_heat, 0.0u"J", max_latent_heat)
    new_mean_temperature = freezing_point - (latent_heat - accumulated_latent_heat[i]) / heat_per_kelvin
    shift = new_mean_temperature - mean_temperature
    update_layer_nodes!(node_temperature -> node_temperature + shift, temperature, i)
    return temperature
end
