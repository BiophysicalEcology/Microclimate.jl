"""
    BlackBodyLongwaveExchange(; atmospheric_radiation_model=CampbellNormanAtmosphericRadiation())

Multi-layer longwave exchange in which leaves send longwave as black bodies at their own
temperature, as micropoint's `longwavemodelCpp` does. Each layer intercepts
`1 - exp(-plant area index)` and replaces it with black-body leaf longwave, so no layer
reflects, but leaves absorb only `leaf_emissivity` of what they intercept. The leaf energy
balance therefore closes only with `leaf_emissivity = 1`, as in micropoint. The ground
sends up its own surface longwave, emitted plus reflected, so unlike micropoint's its
longwave balance closes.

# References
- Maclean, I. M. D. micropoint: point microclimate model, R package,
  https://github.com/ilyamaclean/micropoint (commit 57b23f3).
"""
@kwdef struct BlackBodyLongwaveExchange{ARM<:AbstractAtmosphericRadiationModel} <: AbstractCanopyLongwaveModel
    atmospheric_radiation_model::ARM = CampbellNormanAtmosphericRadiation()
end

# The same layers as LayeredLongwaveExchange
allocate_longwave(::BlackBodyLongwaveExchange, plant_area_index, n_layers, canopy_projection_ratio) =
    allocate_longwave(LayeredLongwaveExchange(), plant_area_index, n_layers, canopy_projection_ratio)

"""
    canopy_longwave!(buffers, model::BlackBodyLongwaveExchange, leaf_emissivity;
                      leaf_temperature, ground, incoming_longwave)

One pass down, one pass up. Writes each layer's absorbed longwave per ground area and
returns `(; ground_absorbed_longwave, net_absorbed_longwave, ground_downward_longwave)`.
"""
function canopy_longwave!(buffers, model::BlackBodyLongwaveExchange, leaf_emissivity;
    leaf_temperature, ground, incoming_longwave,
)
    (; layer_transmission, boundary_downward_longwave, boundary_upward_longwave, absorbed_longwave) = buffers
    n_layers = length(layer_transmission)

    boundary_downward_longwave[1] = incoming_longwave
    for i in 1:n_layers
        boundary_downward_longwave[i + 1] = black_body_layer_flux(boundary_downward_longwave[i], layer_transmission[i], leaf_temperature[i])
    end

    ground_downward_longwave = boundary_downward_longwave[n_layers + 1]
    boundary_upward_longwave[n_layers + 1] = ground.emitted_longwave + ground.longwave_reflectance * ground_downward_longwave
    for i in n_layers:-1:1
        boundary_upward_longwave[i] = black_body_layer_flux(boundary_upward_longwave[i + 1], layer_transmission[i], leaf_temperature[i])
    end

    for i in 1:n_layers
        absorbed_longwave[i] = leaf_emissivity * (1.0 - layer_transmission[i]) *
            (boundary_downward_longwave[i] + boundary_upward_longwave[i + 1])
    end

    ground_absorbed_longwave = (1.0 - ground.longwave_reflectance) * ground_downward_longwave
    net_absorbed_longwave = boundary_downward_longwave[1] - boundary_upward_longwave[1]

    return (; ground_absorbed_longwave, net_absorbed_longwave, ground_downward_longwave)
end

# A layer passes on what it transmits and replaces what it intercepts with black-body leaf longwave
black_body_layer_flux(incoming, transmission, leaf_temperature) =
    transmission * incoming + (1.0 - transmission) * σ * u"K"(leaf_temperature)^4

# Layers don't reflect, so nothing the ground sends up comes back down
function ground_longwave_conditions(buffers, model::BlackBodyLongwaveExchange, leaf_emissivity;
    leaf_temperature, incoming_longwave,
)
    (; layer_transmission) = buffers
    downward_longwave = incoming_longwave
    for i in eachindex(layer_transmission)
        downward_longwave = black_body_layer_flux(downward_longwave, layer_transmission[i], leaf_temperature[i])
    end
    return (; downward_longwave, longwave_return_fraction=0.0)
end
