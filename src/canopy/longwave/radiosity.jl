"""
    LayeredRadiosityExchange(; atmospheric_radiation_model=CampbellNormanAtmosphericRadiation())

Multi-layer longwave exchange as a stack of reflecting, transmitting,
emitting slabs, solved implicitly in one pass. Each layer's gap
transmission `τᵢ` (probability of no interception) leaves a scattered
remainder `(1-leaf_emissivity)(1-τᵢ)`, split into a backscattered part
(returns the way it came) and a forward-scattered part (folded into an
effective transmission alongside `τᵢ`) by `backscatter_fraction` — Shaw's
leaf-angle-dependent `FDDU` (Flerchinger & Yu 2007), computed from
`canopy_projection_ratio`. `leaf_emissivity=1` (no scattering) reduces
this to [`LayeredLongwaveExchange`](@ref)'s cascade. Non-adjacent-layer
coupling (repeated reflection between layers/ground) falls out of the
implicit solve rather than needing an explicit all-pairs sum, unlike
[`AllPairsLongwaveExchange`](@ref).

Solved via the two-stream "adding" recursion: positing
`upward[i] = reflectance_response[i] * downward[i] + effective_source[i]`
and sweeping ground-to-sky gives `reflectance_response`/`response_gain` as
pure layer geometry and `effective_source` as the only temperature-
dependent piece; a single sky-to-ground pass then recovers every boundary
flux. `O(n_layers)` per call, matching [`LayeredLongwaveExchange`](@ref)'s
cost while adding reflection. Backscatter/forward-scatter split by
construction (`ρᵢ+tᵢ+εᵢ(1-τᵢ)=1`, both hemispheres) keeps this energy-
conserving and Kirchhoff-consistent for any `backscatter_fraction`, not
just full backscatter.

Loosely follows Flerchinger's SHAW canopy longwave scheme (implicit banded
solve of reflecting/transmitting/emitting layers).

# References
- Flerchinger, G. N. and Saxton, K. E. (1989). Simultaneous heat and water
  model of a freezing snow-residue-soil system I: theory and development.
  *Transactions of the ASAE*, 32, 565-571.
- Flerchinger, G. N. and Pierson, F. B. (1991). Modeling plant canopy
  effects on variability of soil temperature and water. *Agricultural and
  Forest Meteorology*, 56, 227-246.
- Flerchinger, G. N. (2000). *The Simultaneous Heat and Water (SHAW) Model:
  Technical Documentation*. USDA-ARS Northwest Watershed Research Center.
- Flerchinger, G. N. and Yu, Q. (2007). Simplified expressions for
  diffuse and scattered radiation transfer in discontinuous vegetation
  canopies. *Agronomy Journal*, 99(1), 243-250.
- Norman, J. M. (1979). Modeling the complete crop canopy. In
  *Modification of the Aerial Environment of Crops*, ASAE, 249-277.
"""
@kwdef struct LayeredRadiosityExchange{ARM<:AbstractAtmosphericRadiationModel} <: AbstractCanopyLongwaveModel
    atmospheric_radiation_model::ARM = CampbellNormanAtmosphericRadiation()
end

function allocate_longwave(::LayeredRadiosityExchange, plant_area_index, n_layers, canopy_projection_ratio)
    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    layer_transmission = exp.(-layer_plant_area_index)

    # Shaw's FDDU (Flerchinger & Yu 2007): diffuse upward-backscatter fraction
    x = canopy_projection_ratio
    backscatter_fraction = 0.5 + 0.5 * (atan(x) / (π / 2))^1.585

    return (;
        layer_transmission, backscatter_fraction,
        boundary_downward_longwave = zeros(typeof(0.0u"W/m^2"), n_layers + 1),
        boundary_upward_longwave = zeros(typeof(0.0u"W/m^2"), n_layers + 1),
        absorbed_longwave = zeros(typeof(0.0u"W/m^2"), n_layers),
        layer_reflectance = zeros(n_layers),
        layer_effective_transmission = zeros(n_layers),
        layer_emission_source = zeros(typeof(0.0u"W/m^2"), n_layers),
        reflectance_response = zeros(n_layers + 1),
        response_gain = zeros(n_layers),
        effective_source = zeros(typeof(0.0u"W/m^2"), n_layers + 1),
    )
end

"""
    canopy_longwave!(buffers, model::LayeredRadiosityExchange, leaf_emissivity;
                      leaf_temperature, ground_temperature, ground_emissivity,
                      site, environment_instant, vapour_pressure_equation=GoffGratch())

Layer index 1 is canopy top. Returns `(; ground_absorbed_longwave,
landscape_absorbed_longwave)`.
"""
function canopy_longwave!(buffers, model::LayeredRadiosityExchange, leaf_emissivity;
    leaf_temperature, ground_temperature, ground_emissivity,
    site, environment_instant, vapour_pressure_equation=GoffGratch(),
)
    (; layer_transmission, backscatter_fraction, boundary_downward_longwave, boundary_upward_longwave,
       absorbed_longwave, layer_reflectance, layer_effective_transmission, layer_emission_source,
       reflectance_response, response_gain, effective_source) = buffers
    n_layers = length(layer_transmission)

    sky = precompute_longwave_sky(model.atmospheric_radiation_model;
        site, environment_instant, vapour_pressure_equation, shade=0.0)
    incoming_longwave = sky.incoming_longwave
    ground_emission = ground_emissivity * σ * u"K"(ground_temperature)^4

    @inbounds for i in 1:n_layers
        τ_gap = layer_transmission[i]
        scattered = (1.0 - leaf_emissivity) * (1.0 - τ_gap)
        layer_reflectance[i] = backscatter_fraction * scattered
        layer_effective_transmission[i] = τ_gap + (1.0 - backscatter_fraction) * scattered
        layer_emission_source[i] = leaf_emissivity * (1.0 - τ_gap) * σ * u"K"(leaf_temperature[i])^4
    end

    # ground-to-sky sweep: upward[i] = reflectance_response[i]*downward[i] + effective_source[i]
    reflectance_response[n_layers + 1] = 1.0 - ground_emissivity
    effective_source[n_layers + 1] = ground_emission
    @inbounds for i in n_layers:-1:1
        t = layer_effective_transmission[i]
        ρ = layer_reflectance[i]
        R_next = reflectance_response[i + 1]
        S_next = effective_source[i + 1]
        γ = 1.0 / (1.0 - ρ * R_next)
        response_gain[i] = γ
        reflectance_response[i] = ρ + t^2 * γ * R_next
        effective_source[i] = t * γ * R_next * (ρ * S_next + layer_emission_source[i]) +
            t * S_next + layer_emission_source[i]
    end

    # sky-to-ground sweep: recover every boundary flux from the precomputed response/gain/source
    boundary_downward_longwave[1] = incoming_longwave
    @inbounds for i in 1:n_layers
        D = boundary_downward_longwave[i]
        boundary_upward_longwave[i] = reflectance_response[i] * D + effective_source[i]
        boundary_downward_longwave[i + 1] = response_gain[i] *
            (layer_effective_transmission[i] * D + layer_reflectance[i] * effective_source[i + 1] + layer_emission_source[i])
    end
    boundary_upward_longwave[n_layers + 1] = reflectance_response[n_layers + 1] * boundary_downward_longwave[n_layers + 1] +
        effective_source[n_layers + 1]

    @inbounds for i in 1:n_layers
        absorbed_longwave[i] = leaf_emissivity * (1.0 - layer_transmission[i]) *
            (boundary_downward_longwave[i] + boundary_upward_longwave[i + 1])
    end

    ground_absorbed_longwave = ground_emissivity * boundary_downward_longwave[n_layers + 1]
    landscape_absorbed_longwave = boundary_downward_longwave[1] - boundary_upward_longwave[1]

    return (; ground_absorbed_longwave, landscape_absorbed_longwave)
end
