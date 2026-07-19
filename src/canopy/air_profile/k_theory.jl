"""
    KTheoryAirProfile()

Steady first-order-closure (K-theory) in-canopy air-temperature profile: a
tridiagonal diffusion balance (finite-volume discretisation of
`d/dz(K dT/dz) = -S`) between each layer's leaf-supplied sensible-heat
source and its neighbours, re-solved once per Picard pass. Assumes in-canopy
air reaches quasi-equilibrium within an hour, unlike Raupach's (1989)
transient Lagrangian near/far-field model.

Conductance network mirrors `SoilHeatTransport1D`'s `thermal_conductance`
(eddy-diffusivity instead of thermal), solved via `LinearAlgebra`'s
tridiagonal solver. Unitful quantities are stripped/reattached around the
solve (generic tridiagonal LU needs a dimensionless `oneunit`), the same
boundary `RootFindLeafTemperature` draws around `HeatExchange.zbrent`.

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
    layer_heights = sort(heights[heights .<= canopy_height]; rev=true)  # top-to-bottom, matching radiation/wind convention
    length(layer_heights) == n_layers ||
        throw(ArgumentError("number of `heights` at or below `canopy_height` must equal n_layers"))

    # Floored: node heights can legitimately coincide with a boundary (e.g. the
    # existing shortwave/wind tests put a node exactly at canopy_height — fine
    # for PAI-based physics, but a literal zero gap here would divide by zero.
    min_spacing = 1.0e-3u"m"
    layer_spacing = similar(layer_heights)
    @inbounds for i in 1:(n_layers - 1)
        layer_spacing[i] = max(layer_heights[i] - layer_heights[i + 1], min_spacing)
    end
    layer_spacing[n_layers] = max(layer_heights[n_layers], min_spacing)  # bottom layer down to the ground (z=0)
    top_spacing = max(canopy_height - layer_heights[1], min_spacing)

    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    relative_eddy_diffusivity = wind_attenuation_profile(
        layer_plant_area_index, canopy_height, plant_area_index, boundary_layer_model.karman_constant,
    )

    return (;
        layer_spacing, top_spacing, relative_eddy_diffusivity,
        dl = zeros(n_layers - 1), d = zeros(n_layers), du = zeros(n_layers - 1), rhs = zeros(n_layers),
        air_temperature = zeros(typeof(0.0u"K"), n_layers),
    )
end

function canopy_air_profile!(buffers, ::KTheoryAirProfile, boundary_layer_model;
    canopy_height, displacement_height, friction_velocity,
    canopy_top_air_temperature, ground_temperature, sensible_heat_source,
)
    (; layer_spacing, top_spacing, relative_eddy_diffusivity, dl, d, du, rhs, air_temperature) = buffers
    n = length(air_temperature)
    length(sensible_heat_source) == n || throw(ArgumentError("sensible_heat_source must have one entry per canopy layer"))

    ρ_cp = calc_ρ_cp(canopy_top_air_temperature)
    eddy_diffusivity_top = boundary_layer_model.karman_constant * friction_velocity *
        max(canopy_height - displacement_height, 1.0e-3u"m")

    # Boundary conductances g[0] (canopy top <-> layer 1) .. g[n] (layer n <-> ground), W/m^2/K.
    g_top = ustrip(u"W/m^2/K", ρ_cp * eddy_diffusivity_top * relative_eddy_diffusivity[1] / top_spacing)
    g_ground = ustrip(u"W/m^2/K", ρ_cp * eddy_diffusivity_top * relative_eddy_diffusivity[n] / layer_spacing[n])
    @inbounds for i in 1:(n - 1)
        eddy_diffusivity_boundary = eddy_diffusivity_top * 0.5 * (relative_eddy_diffusivity[i] + relative_eddy_diffusivity[i + 1])
        g = ustrip(u"W/m^2/K", ρ_cp * eddy_diffusivity_boundary / layer_spacing[i])
        dl[i] = g
        du[i] = g
    end

    # Finite-volume flux balance per layer: g_prev*(T_prev - T_i) - g_next*(T_i - T_next) + S_i = 0
    # (T_prev/T_next are the known boundary values at the canopy top/ground for layers 1/n).
    T_top = ustrip(u"K", canopy_top_air_temperature)
    T_ground = ustrip(u"K", ground_temperature)
    @inbounds for i in 1:n
        g_prev = i == 1 ? g_top : dl[i - 1]
        g_next = i == n ? g_ground : dl[i]
        d[i] = -(g_prev + g_next)
        rhs[i] = -ustrip(u"W/m^2", sensible_heat_source[i])
    end
    rhs[1] -= g_top * T_top
    rhs[n] -= g_ground * T_ground

    air_temperature .= (Tridiagonal(dl, d, du) \ rhs) .* u"K"

    return (; air_temperature)
end
