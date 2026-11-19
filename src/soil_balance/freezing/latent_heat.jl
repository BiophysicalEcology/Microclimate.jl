"""
    PhaseTransitionLatentHeat()

Soil ice/water phase-change accounting (as in NicheMapR). Each
layer tracks accumulated latent heat against `mass × L_fusion`. While the
budget is non-empty, layer temperatures are clamped to 0°C; once exhausted,
the clamp is released and the layer cools/warms freely.
"""
struct PhaseTransitionLatentHeat <: SoilPhaseTransitionModel end

# Fortran OSUB.f:571-602 fires these branches ONLY on an actual crossing. Once a layer has
# been clamped to 0°C, its past mean is ≈ 0 and neither branch fires: the layer then
# cools/warms freely on subsequent hours.
function layer_phase_transition!(::PhaseTransitionLatentHeat, temperature, accumulated_latent_heat, i;
    layer_mass, mean_temperature, mean_temperature_past,
)
    tolerance = 1.0e-4u"K"
    freezing_point = u"K"(0.0u"°C")
    max_latent_heat = LATENT_HEAT_FUSION * layer_mass
    to_freezing_point(_) = freezing_point
    # Freezing: above → below 0°C
    if (mean_temperature_past - freezing_point > tolerance) && (mean_temperature - freezing_point <= -tolerance)
        accumulated_latent_heat[i] += (mean_temperature_past - mean_temperature) * layer_mass * water_heat_capacity
        if accumulated_latent_heat[i] >= max_latent_heat
            accumulated_latent_heat[i] = max_latent_heat
        end
        update_layer_nodes!(to_freezing_point, temperature, i)
    # Thawing: below → above 0°C
    elseif (mean_temperature_past - freezing_point < -tolerance) && (mean_temperature - freezing_point >= tolerance)
        accumulated_latent_heat[i] -= (mean_temperature - mean_temperature_past) * layer_mass * water_heat_capacity
        if accumulated_latent_heat[i] <= 0.0u"J"
            accumulated_latent_heat[i] = 0.0u"J"
        else
            update_layer_nodes!(to_freezing_point, temperature, i)
        end
    end
    return temperature
end
