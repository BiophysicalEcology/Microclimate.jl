"""
    AllPairsLongwaveExchange(; atmospheric_radiation_model=CampbellNormanAtmosphericRadiation())

Multi-layer longwave exchange via a pairwise `exp(-ΔPAI)` gap-transmission
weight matrix between every layer, sky, and ground, rescaled per layer so
each layer's total view sums to a full sphere. Ported from micropoint's
`twostreamvegCpp`/microclimlearn's `longwavebelow` for direct comparison
against those references — self-view is omitted (zero diagonal) and ground
does not reflect, matching that source exactly rather than a from-scratch
design; see [`LayeredRadiosityExchange`](@ref) for a reflection- and
self-view-consistent alternative. `O(n_layers^2)` per call; the weight
matrix is precomputed once in [`allocate_longwave`](@ref).

`boundary_downward_longwave`/`boundary_upward_longwave` (for
[`LayeredLongwaveExchange`](@ref)-style output/diagnostics) are `[1:n_layers]
== downward_longwave/upward_longwave` exactly — same layer positions, not a
separate boundary grid — with only `[n_layers+1]` a genuinely distinct
point (the ground surface itself).

# References
- Maclean, I. M. D. micropoint/microclimlearn below-canopy longwave scheme.
"""
@kwdef struct AllPairsLongwaveExchange{ARM<:AbstractAtmosphericRadiationModel} <: AbstractCanopyLongwaveModel
    atmospheric_radiation_model::ARM = CampbellNormanAtmosphericRadiation()
end

# Row-normalized exp(-ΔPAI) weights from each of `receiver_pai_above` (n_receivers) to every
# `source_pai_above`/`source_pai` layer, rescaled per row so canopy+sky+ground view sums to 2.
# `self_view` toggles the zero-diagonal (only meaningful when receivers == layers themselves).
function _view_factor_weights(receiver_pai_above, source_pai_above, source_pai, pai_total, self_view::Bool)
    n_receivers = length(receiver_pai_above)
    n_layers = length(source_pai_above)
    τ_sky = zeros(n_receivers)
    τ_ground = zeros(n_receivers)
    @inbounds for i in 1:n_receivers
        τ_sky[i] = exp(-receiver_pai_above[i])
        τ_ground[i] = exp(-(pai_total - receiver_pai_above[i]))
    end
    w = zeros(n_receivers, n_layers)
    @inbounds for j in 1:n_layers, i in 1:n_receivers
        w[i, j] = (self_view && i == j) ? 0.0 :
            exp(-abs(receiver_pai_above[i] - source_pai_above[j])) * source_pai[j]
    end
    @inbounds for i in 1:n_receivers
        Σw = zero(eltype(w))
        for j in 1:n_layers
            Σw += w[i, j]
        end
        Σv = Σw + τ_sky[i] + τ_ground[i] # total view: canopy + sky + ground
        μ = 1.0 + (2.0 - Σv) / _safe_nonzero(Σw)
        for j in 1:n_layers
            w[i, j] *= μ
        end
    end
    return w, τ_sky, τ_ground
end

function allocate_longwave(::AllPairsLongwaveExchange, plant_area_index, n_layers, canopy_projection_ratio)
    (; layer_plant_area_index, layer_plant_area_index_boundaries) = canopy_layer_geometry(plant_area_index, n_layers)
    pai_total = sum(layer_plant_area_index)

    # layer position = cumulative PAI at the BOTTOM of the layer (boundaries[i+1]),
    # matching R's paia = rev(cumsum(rev(paii))) exactly (not the layer midpoint)
    layer_pai_above = layer_plant_area_index_boundaries[2:end]

    view_weights, sky_transmission, ground_transmission =
        _view_factor_weights(layer_pai_above, layer_pai_above, layer_plant_area_index, pai_total, true)
    # single extra receiver point at the ground surface itself (pai_above = pai_total)
    ground_view_weights, ground_sky_transmission, ground_ground_transmission =
        _view_factor_weights([pai_total], layer_pai_above, layer_plant_area_index, pai_total, false)

    layer_transmission = exp.(-layer_plant_area_index)

    return (;
        view_weights, sky_transmission, ground_transmission,
        ground_view_weights=vec(ground_view_weights), ground_sky_transmission=ground_sky_transmission[1],
        ground_ground_transmission=ground_ground_transmission[1],
        layer_transmission,
        downward_longwave = zeros(typeof(0.0u"W/m^2"), n_layers),
        upward_longwave = zeros(typeof(0.0u"W/m^2"), n_layers),
        boundary_downward_longwave = zeros(typeof(0.0u"W/m^2"), n_layers + 1),
        boundary_upward_longwave = zeros(typeof(0.0u"W/m^2"), n_layers + 1),
        absorbed_longwave = zeros(typeof(0.0u"W/m^2"), n_layers),
        leaf_emission = zeros(typeof(0.0u"W/m^2"), n_layers),
    )
end

"""
    canopy_longwave!(buffers, model::AllPairsLongwaveExchange, leaf_emissivity;
                      leaf_temperature, ground_temperature, ground_emissivity,
                      site, environment_instant, vapour_pressure_equation=GoffGratch())

Layer index 1 is canopy top. Downward flux sums the sky term and every
layer strictly above the receiver; upward flux sums the ground term (no
ground reflection, matching the ported source) and every layer at or
below. `boundary_downward_longwave[1:n_layers]`/`boundary_upward_longwave
[1:n_layers]` equal `downward_longwave`/`upward_longwave` exactly; index
`n_layers+1` is the ground surface. Returns `(; ground_absorbed_longwave,
net_absorbed_longwave)`.
"""
function canopy_longwave!(buffers, model::AllPairsLongwaveExchange, leaf_emissivity;
    leaf_temperature, ground_temperature, ground_emissivity,
    site, environment_instant, vapour_pressure_equation=GoffGratch(),
)
    (; view_weights, sky_transmission, ground_transmission,
       ground_view_weights, ground_sky_transmission, ground_ground_transmission,
       layer_transmission, downward_longwave, upward_longwave,
       boundary_downward_longwave, boundary_upward_longwave, absorbed_longwave, leaf_emission) = buffers
    n_layers = length(layer_transmission)

    sky = precompute_longwave_sky(model.atmospheric_radiation_model;
        site, environment_instant, vapour_pressure_equation, shade=0.0)
    incoming_longwave = sky.incoming_longwave
    ground_emission = ground_emissivity * σ * u"K"(ground_temperature)^4

    @inbounds for i in 1:n_layers
        leaf_emission[i] = leaf_emissivity * σ * u"K"(leaf_temperature[i])^4
    end

    @inbounds for i in 1:n_layers
        down = incoming_longwave * sky_transmission[i]
        up = ground_emission * ground_transmission[i]
        for j in 1:n_layers
            wij = view_weights[i, j]
            j < i && (down += wij * leaf_emission[j])
            j > i && (up += wij * leaf_emission[j])
        end
        downward_longwave[i] = down
        upward_longwave[i] = up
    end

    @inbounds for i in 1:n_layers
        absorbed_longwave[i] = (1.0 - layer_transmission[i]) * (downward_longwave[i] + upward_longwave[i])
        boundary_downward_longwave[i] = downward_longwave[i]
        boundary_upward_longwave[i] = upward_longwave[i]
    end

    # ground surface itself: every layer is "above" it (all contribute to downward);
    # nothing contributes to upward beyond the ground's own emission (no reflection)
    ground_downward = incoming_longwave * ground_sky_transmission
    @inbounds for j in 1:n_layers
        ground_downward += ground_view_weights[j] * leaf_emission[j]
    end
    boundary_downward_longwave[n_layers + 1] = ground_downward
    boundary_upward_longwave[n_layers + 1] = ground_emission * ground_ground_transmission

    ground_absorbed_longwave = ground_emissivity * boundary_downward_longwave[n_layers + 1]
    net_absorbed_longwave = boundary_downward_longwave[1] - boundary_upward_longwave[1]

    return (; ground_absorbed_longwave, net_absorbed_longwave)
end
