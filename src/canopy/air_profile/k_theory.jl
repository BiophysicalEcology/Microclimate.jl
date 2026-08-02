"""
    KTheoryAirProfile()

Steady first-order-closure (K-theory) in-canopy air-temperature and humidity
profile: a tridiagonal diffusion balance (finite-volume discretisation of
`d/dz(K dT/dz) = -S`, and separately `d/dz(K dρv/dz) = -Sv`) between each
layer's leaf-supplied source and its neighbours, re-solved once per Picard
pass. Assumes in-canopy air reaches quasi-equilibrium within an hour, unlike
Raupach's (1989) transient Lagrangian near/far-field model.

Conductance network mirrors `SoilHeatTransport1D`'s `thermal_conductance`
(eddy-diffusivity instead of thermal), solved via `LinearAlgebra`'s
tridiagonal solver. Unitful quantities are stripped/reattached around the
solve (generic tridiagonal LU needs a dimensionless `oneunit`), the same
boundary `RootFindLeafTemperature` draws around `HeatExchange.zbrent`. The
humidity solve reuses the same diffusivity network and buffers as the
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
    (; layer_heights, layer_thickness) = canopy_layer_heights(heights, canopy_height, n_layers)
    layer_spacing = layer_thickness
    top_spacing = max(canopy_height - layer_heights[1], 1.0e-3u"m")  # numerical floor, not a physical parameter

    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    relative_eddy_diffusivity = wind_attenuation_profile(
        layer_plant_area_index, canopy_height, sum(plant_area_index), boundary_layer_model.karman_constant,
    )

    return (;
        layer_spacing, top_spacing, relative_eddy_diffusivity,
        dl = zeros(n_layers - 1), d = zeros(n_layers), du = zeros(n_layers - 1), rhs = zeros(n_layers),
        dl_v = zeros(n_layers - 1), d_v = zeros(n_layers), du_v = zeros(n_layers - 1), rhs_v = zeros(n_layers),
        air_temperature = zeros(typeof(0.0u"K"), n_layers),
        vapour_density = zeros(typeof(0.0u"kg/m^3"), n_layers),
        relative_humidity = zeros(n_layers),
    )
end

# Shared conductance-building loop for both the heat and vapor solves — `scale`
# converts the raw eddy-diffusivity/spacing ratio into the right flux law
# (ρ_cp for a heat flux, 1 for a vapor-density flux already being a mass flux).
function _ktheory_conductances!(dl, du, relative_eddy_diffusivity, eddy_diffusivity_top, layer_spacing, top_spacing, scale, unit)
    n = length(dl) + 1
    g_top = ustrip(unit, scale * eddy_diffusivity_top * relative_eddy_diffusivity[1] / top_spacing)
    g_ground = ustrip(unit, scale * eddy_diffusivity_top * relative_eddy_diffusivity[n] / layer_spacing[n])
    @inbounds for i in 1:(n - 1)
        eddy_diffusivity_boundary = eddy_diffusivity_top * 0.5 * (relative_eddy_diffusivity[i] + relative_eddy_diffusivity[i + 1])
        g = ustrip(unit, scale * eddy_diffusivity_boundary / layer_spacing[i])
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
)
    (; layer_spacing, top_spacing, relative_eddy_diffusivity, dl, d, du, rhs,
       dl_v, d_v, du_v, rhs_v, air_temperature, vapour_density, relative_humidity) = buffers
    n = length(air_temperature)
    length(sensible_heat_source) == n || throw(ArgumentError("sensible_heat_source must have one entry per canopy layer"))

    ρ_cp = calc_ρ_cp(canopy_top_air_temperature)
    eddy_diffusivity_top = boundary_layer_model.karman_constant * friction_velocity *
        max(canopy_height - displacement_height, 1.0e-3u"m")  # numerical floor, not a physical parameter

    # Heat solve — boundary conductances g[0] (canopy top <-> layer 1) .. g[n] (layer n <-> ground), W/m^2/K.
    g_top, g_ground = _ktheory_conductances!(dl, du, relative_eddy_diffusivity, eddy_diffusivity_top,
        layer_spacing, top_spacing, ρ_cp, u"W/m^2/K")

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

    # Vapor-density solve — same diffusivity network, no ρ_cp (a vapor-density
    # gradient is already a mass flux, so the conductance is a plain piston
    # velocity, m/s).
    g_top_v, g_ground_v = _ktheory_conductances!(dl_v, du_v, relative_eddy_diffusivity, eddy_diffusivity_top,
        layer_spacing, top_spacing, 1.0, u"m/s")

    ρv_top = wet_air_properties(canopy_top_air_temperature, canopy_top_relative_humidity, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    ρv_ground = wet_air_properties(ground_temperature, ground_relative_humidity, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    ρv_top_stripped = ustrip(u"kg/m^3", ρv_top)
    ρv_ground_stripped = ustrip(u"kg/m^3", ρv_ground)
    @inbounds for i in 1:n
        g_prev = i == 1 ? g_top_v : dl_v[i - 1]
        g_next = i == n ? g_ground_v : dl_v[i]
        d_v[i] = -(g_prev + g_next)
        rhs_v[i] = -ustrip(u"kg/m^2/s", evaporation_mass_flow[i] / 1.0u"m^2")
    end
    rhs_v[1] -= g_top_v * ρv_top_stripped
    rhs_v[n] -= g_ground_v * ρv_ground_stripped

    vapour_density .= (Tridiagonal(dl_v, d_v, du_v) \ rhs_v) .* u"kg/m^3"

    @inbounds for i in 1:n
        ρv_sat = wet_air_properties(air_temperature[i], 1.0, atmospheric_pressure; vapour_pressure_equation).vapour_density
        relative_humidity[i] = clamp(ustrip(u"kg/m^3", vapour_density[i]) / ustrip(u"kg/m^3", ρv_sat), 0.0, 1.0)
    end

    return (; air_temperature, relative_humidity)
end
