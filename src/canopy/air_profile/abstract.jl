"""
    AbstractCanopyAirProfileModel

How a canopy resolves the in-canopy air-temperature profile from per-layer
leaf sensible-heat sources plus canopy-top and ground boundary values.
Concrete variants implement [`allocate_air_profile`](@ref) and
[`canopy_air_profile!`](@ref). Currently: [`KTheoryAirProfile`](@ref)
(steady first-order-closure diffusion).
"""
abstract type AbstractCanopyAirProfileModel end

"""
    allocate_air_profile(model, canopy_height, displacement_height, plant_area_index, heights, n_layers, boundary_layer_model)

Pre-allocate and precompute the structural buffers a model needs.
"""
function allocate_air_profile end

"""
    CanopyAirProfileForcing(; friction_velocity, obukhov_length, wind_speed,
        canopy_top_air_temperature, canopy_top_relative_humidity, ground_temperature,
        ground_relative_humidity, sensible_heat_source, evaporation_mass_flow,
        leaf_temperature, atmospheric_pressure)

Forcing and state for one [`canopy_air_profile!`](@ref) solve: the
above-canopy turbulence, the canopy-top and ground boundary values, and the
per-layer leaf sources and temperatures of the current Picard pass. Every
air-profile model receives the same forcing and reads the fields it needs.
"""
@kwdef struct CanopyAirProfileForcing{FV,OL,WS,TT,TRH,GT,GRH,SH,EM,LT,P}
    friction_velocity::FV
    obukhov_length::OL
    wind_speed::WS
    canopy_top_air_temperature::TT
    canopy_top_relative_humidity::TRH
    ground_temperature::GT
    ground_relative_humidity::GRH
    sensible_heat_source::SH
    evaporation_mass_flow::EM
    leaf_temperature::LT
    atmospheric_pressure::P
end

"""
    canopy_air_profile!(buffers, model, boundary_layer_model, forcing::CanopyAirProfileForcing;
                        site, vapour_pressure_equation=GoffGratch())

Solve the in-canopy air-temperature and humidity profile for this Picard
pass; writes into `buffers`, returns `(; air_temperature, relative_humidity)`.
"""
function canopy_air_profile! end

"""
    canopy_heat_vapor_transport!(buffers, model, boundary_layer_model, forcing::CanopyAirProfileForcing;
                              site, vapour_pressure_equation=GoffGratch())

The model's physical map from leaf sources and boundary values to the
in-canopy air temperature and vapour density, written into
`buffers.air_temperature` and `buffers.vapor_density`. A model whose map
depends on the current air state reads it from `buffers.air_temperature_prev`
and `buffers.vapor_density_prev`. No relaxation, clamping or relative
humidity: [`canopy_air_profile!`](@ref) adds those for the Picard solve.
"""
function canopy_heat_vapor_transport! end

"""
    clamp_canopy_top_air_temperature(model, canopy_top_air_temperature, reference_temperature, ground_temperature)

Any model-specific bound on the canopy-top boundary temperature. No bound by default.
"""
clamp_canopy_top_air_temperature(::AbstractCanopyAirProfileModel, canopy_top_air_temperature, reference_temperature, ground_temperature) =
    canopy_top_air_temperature

"""
    reset_air_profile_scratch!(model, buffers)

Clear any warm-start/relaxation state that shouldn't carry across a
day-boundary reset (see `reset_canopy_scratch!`) — e.g. Aitken Δ²
residual history, which would otherwise seed the new day's first Picard
pass from an unrelated previous day. No-op by default.
"""
reset_air_profile_scratch!(::AbstractCanopyAirProfileModel, buffers) = nothing
