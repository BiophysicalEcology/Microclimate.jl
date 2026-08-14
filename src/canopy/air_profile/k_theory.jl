"""
    KTheoryAirProfile()

Steady first-order-closure (K-theory) in-canopy air-temperature and humidity
profile: a tridiagonal diffusion balance (finite-volume discretisation of
`d/dz(K dT/dz) = -S`, and separately `d/dz(K dρv/dz) = -Sv`) between each
layer's leaf-supplied source and its neighbours, re-solved once per
convergence-loop pass. Assumes in-canopy air reaches quasi-equilibrium within 
an hour, unlike Raupach's (1989) transient Lagrangian near/far-field model.

Layer 1's node sits exactly at `canopy_height` (`canopy_layer_heights`), 
extrapolated there by the above-canopy Monin-Obukhov profile); an explicit
Dirichlet boundary condition. Layers `2:n_layers` are free unknowns,
diffusing from that fixed top value down to `ground_temperature`/
`ground_relative_humidity` at the other end.

Conductance network mirrors `SoilHeatTransport1D`'s `thermal_conductance`
(eddy-diffusivity instead of thermal), solved via a `LinearSolve.jl` cache
built once per buffer allocation and reused every call (`cache.A`/`cache.b`
reassigned, `SciMLBase.solve!` re-factorizes only when `A` actually changed)
The humidity solve reuses the same diffusivity network and buffers as the
temperature solve.

Eddy diffusivity shape follows the wind-attenuation profile
([`CanopyWindAttenuation`](@ref)/`wind_attenuation_profile`), scaled to the
canopy-top value `κ u* (canopy_height - displacement_height)`.

# References
- Raupach, M. R. (1989). A practical Lagrangian method for relating scalar
  concentrations to source distributions in vegetation canopies. *Quarterly
  Journal of the Royal Meteorological Society*, 115(487), 609-632.
"""
struct KTheoryAirProfile <: AbstractCanopyAirProfileModel end

function allocate_air_profile(::KTheoryAirProfile, canopy_height, plant_area_index, heights, n_layers, boundary_layer_model)
    (; layer_thickness) = canopy_layer_heights(heights, canopy_height, n_layers)
    Δz = layer_thickness

    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    K_relative = wind_attenuation_profile(
        layer_plant_area_index, canopy_height, sum(plant_area_index), boundary_layer_model.karman_constant,
    )

    # Layer 1 is a Dirichlet boundary (see docstring) -- only layers 2:n_layers
    # (n_free = n_layers-1 of them) are solved.
    n_free = n_layers - 1
    dl, d, du, rhs = zeros(n_free - 1), zeros(n_free), zeros(n_free - 1), zeros(n_free)
    dl_v, d_v, du_v, rhs_v = zeros(n_free - 1), zeros(n_free), zeros(n_free - 1), zeros(n_free)
    heat_cache = SciMLBase.init(LinearProblem(Tridiagonal(dl, d, du), rhs))
    vapor_cache = SciMLBase.init(LinearProblem(Tridiagonal(dl_v, d_v, du_v), rhs_v))

    return (;
        layer_spacing = Δz, relative_eddy_diffusivity = K_relative,
        heat_lower_diagonal = dl, heat_diagonal = d, heat_upper_diagonal = du, heat_right_hand_side = rhs,
        vapor_lower_diagonal = dl_v, vapor_diagonal = d_v, vapor_upper_diagonal = du_v, vapor_right_hand_side = rhs_v,
        heat_solver_cache = heat_cache, vapor_solver_cache = vapor_cache,
        air_temperature = zeros(typeof(0.0u"K"), n_layers),
        vapor_density = zeros(typeof(0.0u"kg/m^3"), n_layers),
        relative_humidity = zeros(n_layers),
        ground_heat_conductance = Ref(0.0u"W/m^2/K"),
        ground_vapor_conductance = Ref(0.0u"m/s"),
    )
end

# Shared conductance-building loop for both the heat and vapor solves, over
# the n_free = n_layers-1 free layers 2:n_layers. `scale` converts the raw
# eddy-diffusivity/spacing ratio into the right flux law (ρ_cp for a heat
# flux, 1 for a vapor-density flux already being a mass flux). `g_top` is the
# conductance between the fixed layer-1 boundary and layer 2; `g_ground` the
# conductance between the last free layer and the ground boundary.
function _ktheory_conductances!(dl, du, relative_eddy_diffusivity, eddy_diffusivity_top, layer_spacing, scale, unit)
    n_free = length(dl) + 1
    g_top = ustrip(unit, scale * eddy_diffusivity_top * 0.5 *
        (relative_eddy_diffusivity[1] + relative_eddy_diffusivity[2]) / layer_spacing[1])
    g_ground = ustrip(unit, scale * eddy_diffusivity_top * relative_eddy_diffusivity[n_free + 1] / layer_spacing[n_free + 1])
    @inbounds for i in 1:(n_free - 1)
        eddy_diffusivity_boundary = eddy_diffusivity_top * 0.5 *
            (relative_eddy_diffusivity[i + 1] + relative_eddy_diffusivity[i + 2])
        g = ustrip(unit, scale * eddy_diffusivity_boundary / layer_spacing[i + 1])
        dl[i] = g
        du[i] = g
    end
    return g_top, g_ground
end

function canopy_air_profile!(buffers, ::KTheoryAirProfile, boundary_layer_model;
    canopy_height, displacement_height, friction_velocity,
    canopy_top_air_temperature, canopy_top_relative_humidity, ground_temperature, ground_relative_humidity,
    sensible_heat_source, evaporation_mass_flow, atmospheric_pressure, vapour_pressure_equation=GoffGratch(),
    obukhov_length=nothing,  # shared dispatch surface with RaupachLTheoryAirProfile; K-theory's steady closure doesn't need it
    wind_attenuation=nothing,  # ditto -- K-theory already derives its own ground diffusivity from relative_eddy_diffusivity
    leaf_temperature=nothing,  # ditto -- only RaupachLTheoryAirProfile's :bulk mode uses it
)
    (;
        layer_spacing, relative_eddy_diffusivity,
        heat_lower_diagonal, heat_diagonal, heat_upper_diagonal, heat_right_hand_side,
        vapor_lower_diagonal, vapor_diagonal, vapor_upper_diagonal, vapor_right_hand_side,
        heat_solver_cache, vapor_solver_cache,
        air_temperature, vapor_density, relative_humidity,
        ground_heat_conductance, ground_vapor_conductance,
    ) = buffers

    Δz = layer_spacing
    K̂ = relative_eddy_diffusivity

    dl = heat_lower_diagonal
    d  = heat_diagonal
    du = heat_upper_diagonal
    b  = heat_right_hand_side

    dlᵥ = vapor_lower_diagonal
    dᵥ  = vapor_diagonal
    duᵥ = vapor_upper_diagonal
    bᵥ  = vapor_right_hand_side

    heat_cache = heat_solver_cache
    vapor_cache = vapor_solver_cache

    T = air_temperature
    ρᵥ = vapor_density
    RH = relative_humidity

    h = canopy_height
    d₀ = displacement_height
    uₛ = friction_velocity
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
    K₀ = κ * uₛ * max(h - d₀, 1.0e-3u"m") # numerical floor, not a physical parameter

    # Heat solve — boundary conductances g[0] (layer 1 <-> layer 2) .. g[n] (layer n <-> ground), W/m^2/K.
    g_top, g_ground = _ktheory_conductances!(dl, du, K̂, K₀, Δz, ρcp, u"W/m^2/K")
    ground_heat_conductance[] = g_ground * u"W/m^2/K"

    # Finite-volume flux balance per free layer: g₋(T₋ - Tᵢ) - g₊(Tᵢ - T₊) + Sₕᵢ = 0
    # Layer 1 and the ground provide the known boundary values.
    T₀ = ustrip(u"K", T[1])
    T_g = ustrip(u"K", T_ground)
    @inbounds for i in 1:n_free
        g₋ = i == 1      ? g_top    : dl[i - 1]
        g₊ = i == n_free ? g_ground : dl[i]
        d[i] = -(g₋ + g₊)
        b[i] = -ustrip(u"W/m^2", Sₕ[i + 1])
    end

    b[1] -= g_top * T₀
    b[n_free] -= g_ground * T_g

    heat_cache.A = Tridiagonal(dl, d, du)
    heat_cache.b = b

    T[2:end] .= SciMLBase.solve!(heat_cache).u .* u"K"

    # Vapor-density solve — same diffusivity network, no ρ_cp (A vapor-density gradient multiplied by
    # K already gives a mass flux, so g has units of piston velocity, m/s).
    gᵥ_top, gᵥ_ground = _ktheory_conductances!(dlᵥ, duᵥ, K̂, K₀, Δz, 1.0, u"m/s")
    ground_vapor_conductance[] = gᵥ_ground * u"m/s"

    ρᵥ_ground = wet_air_properties(T_ground, RH_ground, p; vapour_pressure_equation).vapour_density
    ρᵥ₀ = ustrip(u"kg/m^3", ρᵥ_top)
    ρᵥ_g = ustrip(u"kg/m^3", ρᵥ_ground)
    @inbounds for i in 1:n_free
        g₋ = i == 1      ? gᵥ_top    : dlᵥ[i - 1]
        g₊ = i == n_free ? gᵥ_ground : dlᵥ[i]
        dᵥ[i] = -(g₋ + g₊)
        bᵥ[i] = -ustrip(u"kg/m^2/s", E[i + 1] / 1.0u"m^2",)
    end

    bᵥ[1] -= gᵥ_top * ρᵥ₀
    bᵥ[n_free] -= gᵥ_ground * ρᵥ_g

    vapor_cache.A = Tridiagonal(dlᵥ, dᵥ, duᵥ)
    vapor_cache.b = bᵥ

    ρᵥ[2:end] .= SciMLBase.solve!(vapor_cache).u .* u"kg/m^3"

    # Diagnostic relative humidity
    @inbounds for i in 1:n
        ρᵥ_sat = wet_air_properties(T[i], 1.0, p; vapour_pressure_equation).vapour_density
        RH[i] = clamp(ρᵥ[i] / ρᵥ_sat, 0.0, 1.0)
    end

    return (; air_temperature = T, relative_humidity = RH,)
end
