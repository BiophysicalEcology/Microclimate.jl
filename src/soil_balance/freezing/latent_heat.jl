"""
    PhaseTransitionLatentHeat()

Soil ice/water phase-change accounting (as in NicheMapR). Each
layer tracks accumulated latent heat against `mass × L_fusion`. While the
budget is non-empty, layer temperatures are clamped to 0°C; once exhausted,
the clamp is released and the layer cools/warms freely.
"""
struct PhaseTransitionLatentHeat <: SoilPhaseTransitionModel end

function allocate_phase_transition(::PhaseTransitionLatentHeat, num_nodes::Int)
    layer_mass = zeros(Float64, num_nodes)u"kg"
    phase_change_heat = zeros(Float64, num_nodes)u"J"
    mean_temperature = zeros(typeof(0.0u"K"), num_nodes)
    mean_temperature_past = zeros(typeof(0.0u"K"), num_nodes)
    # Mutable scratch for `phase_transition!`; reused across hourly calls.
    temperature_scratch = MVector{num_nodes, typeof(0.0u"K")}(undef)
    return (; layer_mass, phase_change_heat, mean_temperature, mean_temperature_past, temperature_scratch)
end

function phase_transition!(::PhaseTransitionLatentHeat,
    buffers::NamedTuple;
    temperatures::AbstractVector,
    temperatures_past::AbstractVector,
    accumulated_latent_heat::AbstractVector,
    soil_moisture::AbstractVector,
    depths::AbstractVector,
)
    (; layer_mass, phase_change_heat, mean_temperature, mean_temperature_past, temperature_scratch) = buffers
    latent_heat_fusion = 333550.0u"J/kg"
    specific_heat_water = 4184.0u"J/kg/K"
    num_nodes = length(depths)
    temperature = temperature_scratch
    @inbounds for i in eachindex(temperatures)
        temperature[i] = temperatures[i]
    end
    tolerance = 1.0e-4u"K"  # tolerance is a temperature *difference*, so use K directly

    for i in 1:num_nodes
        # --- Always compute mean layer temperatures ---
        if i < num_nodes
            mean_temperature[i] = 0.5 * (temperature[i] + temperature[i+1])
            mean_temperature_past[i] = 0.5 * (temperatures_past[i] + temperatures_past[i+1])
        else
            mean_temperature[i] = temperature[i]
            mean_temperature_past[i] = temperatures_past[i]
        end
        # --- Compute layer mass (kg), handle dry layers ---
        if soil_moisture[i] > 0
            if i < num_nodes
                layer_mass[i] = (u"m"(depths[i+1] - depths[i])) * 1000.0u"kg/m" * soil_moisture[i]
            else
                layer_mass[i] = (u"m"(depths[i] + 100.0u"cm" - depths[i])) * 1000.0u"kg/m" * soil_moisture[i]
            end
        else
            layer_mass[i] = 0.0u"kg"
        end
        max_latent_heat = latent_heat_fusion * layer_mass[i]

        # --- If no water, reset and skip ---
        if soil_moisture[i] <= 0.0 || layer_mass[i] <= 0.0u"kg"
            accumulated_latent_heat[i] = 0.0u"J"
            phase_change_heat[i] = 0.0u"J"
            continue
        end

        # Zero-point used for signed comparison of mean temperatures.
        # Expressed in the same absolute units as the mean (K) so `Δ > tolerance`
        # means "the mean is more than 0.0001°C above 0°C".
        zero_c_abs = u"K"(0.0u"°C")

        # Fortran OSUB.f:571-602 fires these branches ONLY on an actual crossing.
        # Once a node has been clamped to 0°C, meanTpast becomes ≈ 0 and neither
        # branch fires — the node then cools/warms freely on subsequent hours.
        # --- FREEZING (above → below 0°C) ---
        if (mean_temperature_past[i] - zero_c_abs > tolerance) && (mean_temperature[i] - zero_c_abs <= -tolerance)
            phase_change_heat[i] = (mean_temperature_past[i] - mean_temperature[i]) * layer_mass[i] * specific_heat_water
            accumulated_latent_heat[i] += phase_change_heat[i]
            if accumulated_latent_heat[i] >= max_latent_heat
                accumulated_latent_heat[i] = max_latent_heat
                phase_change_heat[i] = 0.0u"J"
            end

            temperature[i] = 0.0u"°C"
            if i < num_nodes
                temperature[i+1] = 0.0u"°C"
            end

        # --- THAWING (below → above 0°C) ---
        elseif (mean_temperature_past[i] - zero_c_abs < -tolerance) && (mean_temperature[i] - zero_c_abs >= tolerance)
            phase_change_heat[i] = (mean_temperature[i] - mean_temperature_past[i]) * layer_mass[i] * specific_heat_water
            accumulated_latent_heat[i] -= phase_change_heat[i]

            if accumulated_latent_heat[i] <= 0.0u"J"
                accumulated_latent_heat[i] = 0.0u"J"
                phase_change_heat[i] = 0.0u"J"
            else
                temperature[i] = 0.0u"°C"
                if i < num_nodes
                    temperature[i+1] = 0.0u"°C"
                end
            end

        else
            phase_change_heat[i] = 0.0u"J"
        end
    end
    # `accumulated_latent_heat` and `phase_change_heat` are mutated in place via the
    # passed buffers — caller already has them. We only return the new temperature
    # SVector so we avoid the heap-allocated NamedTuple wrapper.
    return SVector(temperature)
end
