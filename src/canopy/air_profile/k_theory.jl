"""
    KTheoryAirProfile()

Steady first-order-closure (K-theory) in-canopy air-temperature and humidity
profile: a tridiagonal diffusion balance (finite-volume discretisation of
`d/dz(K dT/dz) = -S`, and separately `d/dz(K dρv/dz) = -Sv`) between each
layer's leaf-supplied source and its neighbours, re-solved at every canopy
residual evaluation. Assumes in-canopy air reaches quasi-equilibrium within 
an hour, unlike Raupach's (1989) transient Lagrangian near/far-field model.

Layer 1's node sits exactly at `canopy_height` (`canopy_layer_heights`), 
extrapolated there by the above-canopy Monin-Obukhov profile); an explicit
Dirichlet boundary condition. Layers `2:n_layers` are free unknowns,
diffusing from that fixed top value down to `ground_temperature`/
`ground_relative_humidity` at the other end.

Conductance network mirrors `SoilHeatTransport1D`'s `thermal_conductance`
(eddy-diffusivity instead of thermal). The balance is symmetric tridiagonal,
factorized in place (`ldlt!`) in preallocated buffers each call.
The humidity solve reuses the same diffusivity network and buffers as the
temperature solve.

Eddy diffusivity shape follows the wind-attenuation profile
([`MixingLengthCanopyWindAttenuation`](@ref)/`wind_attenuation_profile`), scaled to the
canopy-top value `κ u* (canopy_height - displacement_height)`.

# References
- Raupach, M. R. (1989). A practical Lagrangian method for relating scalar
  concentrations to source distributions in vegetation canopies. *Quarterly
  Journal of the Royal Meteorological Society*, 115(487), 609-632.
"""
struct KTheoryAirProfile <: AbstractCanopyAirProfileModel end

function allocate_air_profile(::KTheoryAirProfile, canopy_height, displacement_height, plant_area_index, heights, n_layers, boundary_layer_model)
    (; layer_heights, layer_thickness) = canopy_layer_heights(heights, canopy_height, n_layers)
    Δz = layer_thickness

    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    K_relative = wind_attenuation_profile(
        layer_plant_area_index, canopy_height, sum(plant_area_index), boundary_layer_model.karman_constant,
    )

    # Layer 1 is a Dirichlet boundary (see docstring) -- only layers 2:n_layers
    # (n_free = n_layers-1 of them) are solved.
    n_free = n_layers - 1
    d, ev, rhs = zeros(n_free), zeros(n_free - 1), zeros(n_free)
    d_v, ev_v, rhs_v = zeros(n_free), zeros(n_free - 1), zeros(n_free)

    return (;
        canopy_height, displacement_height, plant_area_index = sum(plant_area_index), layer_heights,
        layer_spacing = Δz, relative_eddy_diffusivity = K_relative,
        heat_diagonal = d, heat_off_diagonal = ev, heat_right_hand_side = rhs,
        vapor_diagonal = d_v, vapor_off_diagonal = ev_v, vapor_right_hand_side = rhs_v,
        heat_matrix = SymTridiagonal(d, ev), vapor_matrix = SymTridiagonal(d_v, ev_v),
        air_temperature = zeros(typeof(0.0u"K"), n_layers),
        air_temperature_prev = zeros(typeof(0.0u"K"), n_layers),
        vapor_density = zeros(typeof(0.0u"kg/m^3"), n_layers),
        vapor_density_prev = zeros(typeof(0.0u"kg/m^3"), n_layers),
        relative_humidity = zeros(n_layers),
        ground_heat_conductance = Ref(0.0u"W/m^2/K"),
        ground_vapor_conductance = Ref(0.0u"m/s"),
    )
end

# The profile's own ground conductances, or those from the coupling's resistances
_ktheory_ground_conductances(::ProfileGroundResistance, boundary_layer_model, g_ground, gᵥ_ground, ρcp, forcing) =
    (g_ground, gᵥ_ground)
function _ktheory_ground_conductances(coupling::AbstractGroundCoupling, boundary_layer_model, g_ground, gᵥ_ground, ρcp, forcing)
    ground = ground_resistance(coupling, boundary_layer_model, forcing)
    return ρcp / ground.heat, 1 / ground.vapor
end

# Shared conductance-building loop for both the heat and vapor solves, over
# the n_free = n_layers-1 free layers 2:n_layers. `scale` converts the raw
# eddy-diffusivity/spacing ratio into the right flux law (ρ_cp for a heat
# flux, 1 for a vapor-density flux already being a mass flux). `g_top` is the
# conductance between the fixed layer-1 boundary and layer 2; `g_ground` the
# conductance between the last free layer and the ground boundary, kept in `unit`.
function _ktheory_conductances!(ev, relative_eddy_diffusivity, eddy_diffusivity_top, layer_spacing, scale, unit)
    n_free = length(ev) + 1
    g_top = ustrip(unit, scale * eddy_diffusivity_top * 0.5 *
        (relative_eddy_diffusivity[1] + relative_eddy_diffusivity[2]) / layer_spacing[1])
    g_ground = uconvert(unit, scale * eddy_diffusivity_top * relative_eddy_diffusivity[n_free + 1] / layer_spacing[n_free + 1])
    for i in 1:(n_free - 1)
        eddy_diffusivity_boundary = eddy_diffusivity_top * 0.5 *
            (relative_eddy_diffusivity[i + 1] + relative_eddy_diffusivity[i + 2])
        g = ustrip(unit, scale * eddy_diffusivity_boundary / layer_spacing[i + 1])
        ev[i] = g
    end
    return g_top, g_ground
end

function canopy_heat_vapor_transport!(buffers, ::KTheoryAirProfile, boundary_layer_model,
    forcing::CanopyAirProfileForcing; site, vapour_pressure_equation=GoffGratch(), ground_coupling=ProfileGroundResistance(),
)
    (; above_canopy_column, wind_speed, canopy_top_air_temperature, canopy_top_relative_humidity, ground_temperature,
       ground_relative_humidity, sensible_heat_source, evaporation_mass_flow, atmospheric_pressure) = forcing
    (;
        canopy_height, displacement_height, plant_area_index, layer_heights, layer_spacing, relative_eddy_diffusivity,
        heat_diagonal, heat_off_diagonal, heat_right_hand_side,
        vapor_diagonal, vapor_off_diagonal, vapor_right_hand_side,
        heat_matrix, vapor_matrix,
        air_temperature, air_temperature_prev, vapor_density,
        ground_heat_conductance, ground_vapor_conductance,
    ) = buffers

    Δz = layer_spacing
    K̂ = relative_eddy_diffusivity

    d  = heat_diagonal
    ev = heat_off_diagonal
    b  = heat_right_hand_side

    dᵥ  = vapor_diagonal
    evᵥ = vapor_off_diagonal
    bᵥ  = vapor_right_hand_side

    T = air_temperature
    ρᵥ = vapor_density

    h = canopy_height
    d₀ = displacement_height
    uₛ = above_canopy_column.friction_velocity
    κ = boundary_layer_model.karman_constant
    p = atmospheric_pressure

    T_top = canopy_top_air_temperature
    RH_top = canopy_top_relative_humidity
    T_ground = ground_temperature
    RH_ground = ground_relative_humidity

    Sₕ = sensible_heat_source
    E = evaporation_mass_flow

    n = length(T)
    n_free = n - 1

    length(Sₕ) == n || throw(ArgumentError("sensible_heat_source must have one entry per canopy layer"))
    length(E) == n || throw(ArgumentError("evaporation_mass_flow must have one entry per canopy layer"))

    # Fixed canopy-top boundary (Dirichlet boundary)
    T[1] = T_top
    ρᵥ_top = wet_air_properties(T_top, RH_top, p; vapour_pressure_equation).vapour_density
    ρᵥ[1] = ρᵥ_top

    # Canopy-top eddy diffusivity
    ρcp = calc_ρ_cp(T_top)
    K₀ = κ * uₛ * max(h - d₀, 1.0e-3u"m") * above_canopy_column.enhancement # numerical floor, not a physical parameter

    # Boundary conductances, layer 1↔2 to layer n↔ground: heat W/m²/K, vapor m/s
    g_top, g_ground = _ktheory_conductances!(ev, K̂, K₀, Δz, ρcp, u"W/m^2/K")
    gᵥ_top, gᵥ_ground = _ktheory_conductances!(evᵥ, K̂, K₀, Δz, 1.0, u"m/s")
    g_ground, gᵥ_ground = _ktheory_ground_conductances(ground_coupling, boundary_layer_model, g_ground, gᵥ_ground, ρcp,
        GroundExchangeForcing(; profile_heat_resistance=ρcp / g_ground, profile_vapor_resistance=1 / gᵥ_ground,
            layer_height=layer_heights[n], wind_speed=wind_speed[n], layer_air_temperature=air_temperature_prev[n],
            ground_temperature=T_ground, roughness_height=site.roughness_height, friction_velocity=uₛ, plant_area_index,
            atmospheric_pressure=p))
    ground_heat_conductance[] = g_ground
    ground_vapor_conductance[] = gᵥ_ground

    # Heat solve

    # Finite-volume flux balance per free layer: g₋(T₋ - Tᵢ) - g₊(Tᵢ - T₊) + Sₕᵢ = 0
    # Layer 1 and the ground provide the known boundary values.
    T₀ = ustrip(u"K", T[1])
    T_g = ustrip(u"K", T_ground)
    g_g = ustrip(u"W/m^2/K", g_ground)
    for i in 1:n_free
        g₋ = i == 1      ? g_top    : ev[i - 1]
        g₊ = i == n_free ? g_g : ev[i]
        d[i] = -(g₋ + g₊)
        b[i] = -ustrip(u"W/m^2", Sₕ[i + 1])
    end

    b[1] -= g_top * T₀
    b[n_free] -= g_g * T_g

    ldiv!(ldlt!(heat_matrix), b)
    T[2:end] .= b .* u"K"

    # Vapor-density solve, on the same diffusivity network
    ρᵥ_ground = wet_air_properties(T_ground, RH_ground, p; vapour_pressure_equation).vapour_density
    ρᵥ₀ = ustrip(u"kg/m^3", ρᵥ_top)
    ρᵥ_g = ustrip(u"kg/m^3", ρᵥ_ground)
    gᵥ_g = ustrip(u"m/s", gᵥ_ground)
    for i in 1:n_free
        g₋ = i == 1      ? gᵥ_top    : evᵥ[i - 1]
        g₊ = i == n_free ? gᵥ_g : evᵥ[i]
        dᵥ[i] = -(g₋ + g₊)
        bᵥ[i] = -ustrip(u"kg/m^2/s", E[i + 1] / 1.0u"m^2")
    end

    bᵥ[1] -= gᵥ_top * ρᵥ₀
    bᵥ[n_free] -= gᵥ_g * ρᵥ_g

    ldiv!(ldlt!(vapor_matrix), bᵥ)
    ρᵥ[2:end] .= bᵥ .* u"kg/m^3"

    return nothing
end
