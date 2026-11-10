"""
    LayeredLongwaveExchange(; atmospheric_radiation_model=CampbellNormanAtmosphericRadiation())

Multi-layer longwave exchange. Each canopy layer is a semi-transparent,
isothermal slab with no reflection, following Norman (1979). Leaves absorb and
emit in proportion to `leaf_emissivity`, and longwave they intercept but don't
absorb carries on in the same direction. Layer transmission assumes that
diffuse extinction is almost independent of leaf angle (Goudriaan 1977).

Sky longwave comes from `atmospheric_radiation_model` via
[`precompute_longwave_sky`](@ref) with no shade.

# References
- Goudriaan, J. (1977). *Crop micrometeorology: a simulation study*. PUDOC,
  Wageningen.
- Norman, J. M. (1979). Modeling the complete crop canopy. In *Modification
  of the Aerial Environment of Crops*, ASAE, 249-277.
"""
@kwdef struct LayeredLongwaveExchange{ARM<:AbstractAtmosphericRadiationModel} <: AbstractCanopyLongwaveModel
    atmospheric_radiation_model::ARM = CampbellNormanAtmosphericRadiation()
end

function allocate_longwave(::LayeredLongwaveExchange, plant_area_index, n_layers, canopy_projection_ratio)
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
                      incoming_longwave)

One pass down, one pass up. Writes each layer's absorbed longwave per ground area and
returns `(; ground_absorbed_longwave, net_absorbed_longwave)`.
"""
function canopy_longwave!(buffers, model::LayeredLongwaveExchange, leaf_emissivity;
    leaf_temperature, ground_temperature, ground_emissivity,
    incoming_longwave,
)
    (; layer_transmission, boundary_downward_longwave, boundary_upward_longwave, absorbed_longwave) = buffers
    n_layers = length(layer_transmission)
    L_d = boundary_downward_longwave
    L_u = boundary_upward_longwave

    L_d[1] = incoming_longwave
    for i in 1:n_layers
        absorptance = leaf_emissivity * (1.0 - layer_transmission[i])
        L_d[i + 1] = L_d[i] * (1.0 - absorptance) + absorptance * σ * u"K"(leaf_temperature[i])^4
    end

    ground_emission = ground_emissivity * σ * u"K"(ground_temperature)^4
    L_u[n_layers + 1] = ground_emission + (1.0 - ground_emissivity) * L_d[n_layers + 1]
    for i in n_layers:-1:1
        absorptance = leaf_emissivity * (1.0 - layer_transmission[i])
        L_u[i] = L_u[i + 1] * (1.0 - absorptance) + absorptance * σ * u"K"(leaf_temperature[i])^4
    end

    for i in 1:n_layers
        absorptance = leaf_emissivity * (1.0 - layer_transmission[i])
        absorbed_longwave[i] = absorptance * (L_d[i] + L_u[i + 1])
    end

    ground_absorbed_longwave = ground_emissivity * L_d[n_layers + 1]
    net_absorbed_longwave = L_d[1] - L_u[1]

    return (; ground_absorbed_longwave, net_absorbed_longwave)
end
