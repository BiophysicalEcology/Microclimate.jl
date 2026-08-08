"""
    KTheoryAirProfile()

Steady first-order-closure (K-theory) in-canopy air-temperature and humidity
profile: a tridiagonal diffusion balance (finite-volume discretisation of
`d/dz(K dT/dz) = -S`, and separately `d/dz(K dρv/dz) = -Sv`) between each
layer's leaf-supplied source and its neighbours, re-solved once per
convergence-loop pass. Assumes in-canopy air reaches quasi-equilibrium within an hour, unlike
Raupach's (1989) transient Lagrangian near/far-field model.

Layer 1's node sits exactly at `canopy_height` (`canopy_layer_heights`), the
same point `canopy_top_air_temperature`/`canopy_top_relative_humidity`
describe (extrapolated there by the above-canopy Monin-Obukhov profile) --
so layer 1's air state is an explicit Dirichlet boundary condition, set
directly rather than solved. Only layers `2:n_layers` are free unknowns,
diffusing from that fixed top value down to `ground_temperature`/
`ground_relative_humidity` at the other end.

Conductance network mirrors `SoilHeatTransport1D`'s `thermal_conductance`
(eddy-diffusivity instead of thermal), solved via a `LinearSolve.jl` cache
built once per buffer allocation and reused every call (`cache.A`/`cache.b`
reassigned, `SciMLBase.solve!` re-factorizes only when `A` actually changed)
-- near-zero allocation, since `A` is fixed within an hour (wind/diffusivity
are hoisted out of the Picard/nonlinear-solve loop) and only `b` changes
across passes. Unitful quantities are stripped/reattached around the solve,
the same boundary `RootFindLeafTemperature` draws around `HeatExchange.zbrent`.
The humidity solve reuses the same diffusivity network and buffers as the
temperature solve, without the `ρ_cp` multiply (a vapor-density gradient is
already a mass flux, not a heat flux).

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
    layer_spacing = layer_thickness

    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    relative_eddy_diffusivity = wind_attenuation_profile(
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
        layer_spacing, relative_eddy_diffusivity,
        dl, d, du, rhs, dl_v, d_v, du_v, rhs_v, heat_cache, vapor_cache,
        air_temperature = zeros(typeof(0.0u"K"), n_layers),
        vapour_density = zeros(typeof(0.0u"kg/m^3"), n_layers),
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
)
    (; layer_spacing, relative_eddy_diffusivity, dl, d, du, rhs,
       dl_v, d_v, du_v, rhs_v, heat_cache, vapor_cache, air_temperature, vapour_density, relative_humidity,
       ground_heat_conductance, ground_vapor_conductance) = buffers
    n = length(air_temperature)
    n_free = n - 1
    length(sensible_heat_source) == n || throw(ArgumentError("sensible_heat_source must have one entry per canopy layer"))

    # Layer 1: Dirichlet boundary (see docstring), set directly. Its own
    # source isn't fed into the diffusion system -- the fixed boundary
    # absorbs it without changing its prescribed state.
    air_temperature[1] = canopy_top_air_temperature
    ρv_top = wet_air_properties(canopy_top_air_temperature, canopy_top_relative_humidity, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    vapour_density[1] = ρv_top

    ρ_cp = calc_ρ_cp(canopy_top_air_temperature)
    eddy_diffusivity_top = boundary_layer_model.karman_constant * friction_velocity *
        max(canopy_height - displacement_height, 1.0e-3u"m")  # numerical floor, not a physical parameter

    # Heat solve — boundary conductances g[0] (layer 1 <-> layer 2) .. g[n] (layer n <-> ground), W/m^2/K.
    g_top, g_ground = _ktheory_conductances!(dl, du, relative_eddy_diffusivity, eddy_diffusivity_top,
        layer_spacing, ρ_cp, u"W/m^2/K")
    ground_heat_conductance[] = g_ground * u"W/m^2/K"

    # Finite-volume flux balance per free layer: g_prev*(T_prev - T_i) - g_next*(T_i - T_next) + S_i = 0
    # (T_prev/T_next are the known boundary values at layer 1/ground for the first/last free layers).
    T_top = ustrip(u"K", air_temperature[1])
    T_ground = ustrip(u"K", ground_temperature)
    @inbounds for i in 1:n_free
        g_prev = i == 1 ? g_top : dl[i - 1]
        g_next = i == n_free ? g_ground : dl[i]
        d[i] = -(g_prev + g_next)
        rhs[i] = -ustrip(u"W/m^2", sensible_heat_source[i + 1])
    end
    rhs[1] -= g_top * T_top
    rhs[n_free] -= g_ground * T_ground

    heat_cache.A = Tridiagonal(dl, d, du)
    heat_cache.b = rhs
    air_temperature[2:end] .= SciMLBase.solve!(heat_cache).u .* u"K"

    # Vapor-density solve — same diffusivity network, no ρ_cp (a vapor-density
    # gradient is already a mass flux, so the conductance is a plain piston
    # velocity, m/s).
    g_top_v, g_ground_v = _ktheory_conductances!(dl_v, du_v, relative_eddy_diffusivity, eddy_diffusivity_top,
        layer_spacing, 1.0, u"m/s")
    ground_vapor_conductance[] = g_ground_v * u"m/s"

    ρv_ground = wet_air_properties(ground_temperature, ground_relative_humidity, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    ρv_top_stripped = ustrip(u"kg/m^3", ρv_top)
    ρv_ground_stripped = ustrip(u"kg/m^3", ρv_ground)
    @inbounds for i in 1:n_free
        g_prev = i == 1 ? g_top_v : dl_v[i - 1]
        g_next = i == n_free ? g_ground_v : dl_v[i]
        d_v[i] = -(g_prev + g_next)
        rhs_v[i] = -ustrip(u"kg/m^2/s", evaporation_mass_flow[i + 1] / 1.0u"m^2")
    end
    rhs_v[1] -= g_top_v * ρv_top_stripped
    rhs_v[n_free] -= g_ground_v * ρv_ground_stripped

    vapor_cache.A = Tridiagonal(dl_v, d_v, du_v)
    vapor_cache.b = rhs_v
    vapour_density[2:end] .= SciMLBase.solve!(vapor_cache).u .* u"kg/m^3"

    @inbounds for i in 1:n
        ρv_sat = wet_air_properties(air_temperature[i], 1.0, atmospheric_pressure; vapour_pressure_equation).vapour_density
        relative_humidity[i] = clamp(vapour_density[i] / ρv_sat, 0.0, 1.0)  # same-unit quantities cancel to a bare Float64
    end

    return (; air_temperature, relative_humidity)
end
