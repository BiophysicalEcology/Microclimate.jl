"""
    LayeredRadiosityExchange(; atmospheric_radiation_model=CampbellNormanAtmosphericRadiation())

Longwave through layers that absorb, emit, transmit and reflect, solved exactly in one
sweep up and one down. The reflected share comes from `canopy_projection_ratio`
(Flerchinger & Yu 2007). With leaf emissivity 1 it equals [`LayeredLongwaveExchange`](@ref).

# References
- Flerchinger, G. N., and Yu, Q. 2007. Simplified expressions for radiation 
  scattering in canopies with ellipsoidal leaf angle distributions.
  *Agricultural and Forest Meteorology*. 144 (3–4): 230–35.
- Flerchinger, G. N. and Saxton, K. E. (1989). Simultaneous heat and water
  model of a freezing snow-residue-soil system I: theory and development.
  *Transactions of the ASAE*, 32, 565-571.
- Flerchinger, G. N. and Pierson, F. B. (1991). Modeling plant canopy
  effects on variability of soil temperature and water. *Agricultural and
  Forest Meteorology*, 56, 227-246.
- Flerchinger, G. N. (2000). *The Simultaneous Heat and Water (SHAW) Model:
  Technical Documentation*. USDA-ARS Northwest Watershed Research Center.
- Norman, J. M. (1979). Modeling the complete crop canopy. In
  *Modification of the Aerial Environment of Crops*, ASAE, 249-277.
"""
@kwdef struct LayeredRadiosityExchange{ARM<:AbstractAtmosphericRadiationModel} <: AbstractCanopyLongwaveModel
    atmospheric_radiation_model::ARM = CampbellNormanAtmosphericRadiation()
end

function allocate_longwave(::LayeredRadiosityExchange, plant_area_index, n_layers, canopy_projection_ratio)
    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    layer_transmission = exp.(-layer_plant_area_index)

    x = canopy_projection_ratio
    backscatter_fraction = 0.5 + 0.5 * (atan(x) / (π / 2))^1.585 # (Flerchinger & Yu 2007, eq. 20)

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
                      leaf_temperature, ground, incoming_longwave)

Layer index 1 is canopy top. Returns `(; ground_absorbed_longwave,
net_absorbed_longwave, ground_downward_longwave)`.
"""
function canopy_longwave!(buffers, model::LayeredRadiosityExchange, leaf_emissivity;
    leaf_temperature, ground, incoming_longwave,
)
    (; layer_transmission, backscatter_fraction, boundary_downward_longwave, boundary_upward_longwave,
       absorbed_longwave, layer_reflectance, layer_effective_transmission, layer_emission_source,
       reflectance_response, response_gain, effective_source) = buffers

    n = length(layer_transmission)

    τ_d = layer_transmission
    εₗ = leaf_emissivity

    ρ = layer_reflectance
    τ = layer_effective_transmission
    E = layer_emission_source

    L_d = boundary_downward_longwave
    L_u = boundary_upward_longwave

    A = absorbed_longwave
    R = reflectance_response
    Γ = response_gain
    S = effective_source

    radiosity_layer_optics!(buffers, leaf_emissivity, leaf_temperature)

    # ground-to-sky sweep
    R[n + 1] = ground.longwave_reflectance
    S[n + 1] = ground.emitted_longwave
    for i in n:-1:1
        Γ[i] = inv(1.0 - ρ[i] * R[i + 1])
        R[i] = ρ[i] + τ[i]^2 * Γ[i] * R[i + 1]
        S[i] = E[i] + τ[i] * S[i + 1] + τ[i] * Γ[i] * R[i + 1] * (ρ[i] * S[i + 1] + E[i])
    end

    # sky-to-ground recovery sweep

    L_d[1] = incoming_longwave
    for i in 1:n
        L_u[i] = R[i] * L_d[i] + S[i]
        L_d[i + 1] = Γ[i] * (τ[i] * L_d[i] + ρ[i] * S[i + 1] + E[i])
    end

    # ground boundary
    L_u[n + 1] = R[n + 1] * L_d[n + 1] + S[n + 1]
    for i in 1:n
        A[i] = εₗ * (1.0 - τ_d[i]) * (L_d[i] + L_u[i + 1])
    end

    ground_absorbed_longwave = (1.0 - ground.longwave_reflectance) * L_d[n + 1]
    net_absorbed_longwave = L_d[1] - L_u[1]

    return (; ground_absorbed_longwave, net_absorbed_longwave, ground_downward_longwave=L_d[n + 1])
end

# The reflectance, transmission and emission of each canopy layer
function radiosity_layer_optics!(buffers, leaf_emissivity, leaf_temperature)
    (; layer_transmission, backscatter_fraction, layer_reflectance, layer_effective_transmission, layer_emission_source) = buffers
    for i in eachindex(layer_transmission)
        intercepted = 1.0 - layer_transmission[i]
        scattered = (1.0 - leaf_emissivity) * intercepted
        layer_reflectance[i] = backscatter_fraction * scattered
        layer_effective_transmission[i] = layer_transmission[i] + (1.0 - backscatter_fraction) * scattered
        layer_emission_source[i] = leaf_emissivity * intercepted * σ * u"K"(leaf_temperature[i])^4
    end
    return buffers
end

# The canopy over a black ground that emits nothing sends down `downward_longwave`. Going down
# from the sky, each layer adds to the reflectance the canopy shows the ground.
function ground_longwave_conditions(buffers, model::LayeredRadiosityExchange, leaf_emissivity;
    leaf_temperature, incoming_longwave,
)
    (; layer_reflectance, layer_effective_transmission, layer_emission_source) = radiosity_layer_optics!(buffers,
        leaf_emissivity, leaf_temperature)
    downward_longwave = incoming_longwave
    # Seen from below, the canopy above layer i reflects `reflectance` of the upward longwave
    reflectance = 0.0
    for i in eachindex(layer_reflectance)
        gain = inv(1.0 - layer_reflectance[i] * reflectance)
        # Downward out of layer i, with no upward longwave from below it
        downward_longwave = layer_effective_transmission[i] * gain * (downward_longwave + reflectance * layer_emission_source[i]) + layer_emission_source[i]
        reflectance = layer_reflectance[i] + layer_effective_transmission[i]^2 * gain * reflectance
    end
    return (; downward_longwave, longwave_return_fraction=reflectance)
end
