"""
    LayeredRadiosityExchange(; atmospheric_radiation_model=CampbellNormanAtmosphericRadiation())

Multi-layer longwave exchange: reflecting/transmitting/emitting slabs,
solved implicitly in one pass. Layer `i`'s gap transmission `τ_d[i]`
leaves a scattered remainder `(1-εₗ)(1-τ_d[i])`, split by
`backscatter_fraction` (Flerchinger & Yu 2007, from
`canopy_projection_ratio`) into reflectance `ρ[i]` and an effective
transmission `τ[i]` (forward-scatter folded in). `εₗ=1` reduces this to
[`LayeredLongwaveExchange`](@ref)'s cascade.

Solved via the two-stream "adding" recursion `L_u[i] = R[i]*L_d[i] + S[i]`:
a ground-to-sky sweep gives `R`/`Γ` as pure geometry and `S` as the only
temperature-dependent piece, then one sky-to-ground pass recovers every
`L_d`/`L_u`. `O(n)` per call. Non-adjacent-layer coupling falls out of the 
implicit solve. `ρ[i]+τ[i]+εₗ(1-τ_d[i])=1` keeps this energy-conserving and 
Kirchhoff-consistent for any `backscatter_fraction`.

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
                      leaf_temperature, ground_temperature, ground_emissivity,
                      site, environment_instant, vapour_pressure_equation=GoffGratch())

Layer index 1 is canopy top. Returns `(; ground_absorbed_longwave,
net_absorbed_longwave)`.
"""
function canopy_longwave!(buffers, model::LayeredRadiosityExchange, leaf_emissivity;
    leaf_temperature, ground_temperature, ground_emissivity,
    site, environment_instant, vapour_pressure_equation=GoffGratch(),
)
    (; layer_transmission, backscatter_fraction, boundary_downward_longwave, boundary_upward_longwave,
       absorbed_longwave, layer_reflectance, layer_effective_transmission, layer_emission_source,
       reflectance_response, response_gain, effective_source) = buffers

    n = length(layer_transmission)

    τ_d = layer_transmission; f_d = backscatter_fraction
    εₗ = leaf_emissivity
    εₛ = ground_emissivity;    Tₛ = u"K"(ground_temperature)

    ρ = layer_reflectance;     τ = layer_effective_transmission
    E = layer_emission_source

    L_d = boundary_downward_longwave
    L_u = boundary_upward_longwave

    A = absorbed_longwave;     R = reflectance_response
    Γ = response_gain;         S = effective_source

    sky = precompute_longwave_sky(model.atmospheric_radiation_model;
    site, environment_instant, vapour_pressure_equation, shade=0.0)   

    # optical coefficients of each canopy layer
    @inbounds for i in 1:n
        interceptedᵢ = 1.0 - τ_d[i]
        scatteredᵢ = (1.0 - εₗ) * interceptedᵢ
        αᵢ = εₗ * interceptedᵢ
        Tₗᵢ = u"K"(leaf_temperature[i]) # not broadcasting this avoids allocation
        ρ[i] = f_d * scatteredᵢ
        τ[i] = τ_d[i] + (1.0 - f_d) * scatteredᵢ
        E[i] = αᵢ * σ * Tₗᵢ^4
    end

    # ground-to-sky sweep
    Eₛ = εₛ * σ * Tₛ^4  # ground emission 
    R[n + 1] = 1.0 - εₛ
    S[n + 1] = Eₛ
    @inbounds for i in n:-1:1
        Γ[i] = inv(1.0 - ρ[i] * R[i + 1])
        R[i] = ρ[i] + τ[i]^2 * Γ[i] * R[i + 1]
        S[i] = E[i] + τ[i] * S[i + 1] + τ[i] * Γ[i] * R[i + 1] * (ρ[i] * S[i + 1] + E[i])
    end

    # sky-to-ground recovery sweep

    L_d[1] = sky.incoming_longwave
    @inbounds for i in 1:n
        L_u[i] = R[i] * L_d[i] + S[i]
        L_d[i + 1] = Γ[i] * (τ[i] * L_d[i] + ρ[i] * S[i + 1] + E[i])
    end

    # ground boundary
    L_u[n + 1] = R[n + 1] * L_d[n + 1] + S[n + 1]
    @inbounds for i in 1:n
        A[i] = εₗ * (1.0 - τ_d[i]) * (L_d[i] + L_u[i + 1])
    end

    ground_absorbed_longwave = εₛ * L_d[n + 1]
    net_absorbed_longwave = L_d[1] - L_u[1]

    return (; ground_absorbed_longwave, net_absorbed_longwave)
end
