"""
    LayeredLongwaveExchange(; atmospheric_radiation_model=CampbellNormanAtmosphericRadiation())

Multi-layer longwave exchange. Each canopy layer is a semi-transparent,
isothermal slab that is *opaque on interception* — longwave reaching a leaf
is fully absorbed, not partly reflected (canopy-internal scattering
neglected, standard for ε=0.95-0.98 leaves, Norman 1979); emission still
scales with true `leaf_emissivity` via Stefan-Boltzmann.

Layer transmission (gap fraction) uses Goudriaan's (1977) result that
diffuse extinction is ≈1 per unit plant area index, almost independent of
leaf angle — unlike the direct beam (see
[`ellipsoidal_extinction_coefficient`](@ref)), so no leaf-angle parameter
is needed here.

Sky longwave reuses [`precompute_longwave_sky`](@ref) with `shade=0` — the
canopy resolves vegetation shading layer by layer, so the bulk `shade`
adjustment is switched off. `atmospheric_radiation_model` is the same
clear-sky block [`ViewFactorLongwave`](@ref) uses for bare ground.

Ground temperature/emissivity are supplied per call, not stored — the
lagged (previous-hour) soil-canopy coupling.

# References
- Goudriaan, J. (1977). *Crop micrometeorology: a simulation study*. PUDOC,
  Wageningen.
- Norman, J. M. (1979). Modeling the complete crop canopy. In *Modification
  of the Aerial Environment of Crops*, ASAE, 249-277.
"""
@kwdef struct LayeredLongwaveExchange{ARM<:AbstractAtmosphericRadiationModel} <: AbstractCanopyLongwaveModel
    atmospheric_radiation_model::ARM = CampbellNormanAtmosphericRadiation()
end

function allocate_longwave(::LayeredLongwaveExchange, plant_area_index, n_layers)
    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    layer_transmission = exp.(-layer_plant_area_index)  # Goudriaan (1977): k_diffuse ≈ 1
    return (;
        layer_transmission,
        boundary_downward_longwave = zeros(typeof(0.0u"W/m^2"), n_layers + 1),
        boundary_upward_longwave = zeros(typeof(0.0u"W/m^2"), n_layers + 1),
        absorbed_longwave = zeros(typeof(0.0u"W/m^2"), n_layers),
    )
end

"""
    canopy_longwave!(buffers, model::LayeredLongwaveExchange, leaf_emissivity;
                      leaf_temperature, ground_temperature, ground_emissivity,
                      site, environment_instant, vapour_pressure_equation=GoffGratch())

Sequential top-down (sky → ground) then bottom-up (ground → sky) gap-fraction
exchange, each pass `O(n_layers)` with no matrix solve (non-reflecting
layers mean the two streams don't couple). Writes `absorbed_longwave`
(gross per layer, per unit ground area — `leaf_heat_balance` subtracts the
leaf's own emission separately) and returns `(; ground_absorbed_longwave,
canopy_absorbed_longwave)`.
"""
function canopy_longwave!(buffers, model::LayeredLongwaveExchange, leaf_emissivity;
    leaf_temperature, ground_temperature, ground_emissivity,
    site, environment_instant, vapour_pressure_equation=GoffGratch(),
)
    (; layer_transmission, boundary_downward_longwave, boundary_upward_longwave, absorbed_longwave) = buffers
    n_layers = length(layer_transmission)

    sky = precompute_longwave_sky(model.atmospheric_radiation_model;
        site, environment_instant, vapour_pressure_equation, shade=0.0)
    boundary_downward_longwave[1] = sky.incoming_longwave
    @inbounds for i in 1:n_layers
        τ = layer_transmission[i]
        emission = leaf_emissivity * σ * u"K"(leaf_temperature[i])^4
        boundary_downward_longwave[i + 1] = boundary_downward_longwave[i] * τ + emission * (1.0 - τ)
    end

    ground_emission = ground_emissivity * σ * u"K"(ground_temperature)^4
    boundary_upward_longwave[n_layers + 1] = ground_emission + (1.0 - ground_emissivity) * boundary_downward_longwave[n_layers + 1]
    @inbounds for i in n_layers:-1:1
        τ = layer_transmission[i]
        emission = leaf_emissivity * σ * u"K"(leaf_temperature[i])^4
        boundary_upward_longwave[i] = boundary_upward_longwave[i + 1] * τ + emission * (1.0 - τ)
    end

    @inbounds for i in 1:n_layers
        τ = layer_transmission[i]
        absorbed_longwave[i] = (1.0 - τ) * (boundary_downward_longwave[i] + boundary_upward_longwave[i + 1])
    end

    ground_absorbed_longwave = ground_emissivity * boundary_downward_longwave[n_layers + 1]
    canopy_absorbed_longwave = boundary_downward_longwave[1] - boundary_upward_longwave[1]

    return (; ground_absorbed_longwave, canopy_absorbed_longwave)
end
