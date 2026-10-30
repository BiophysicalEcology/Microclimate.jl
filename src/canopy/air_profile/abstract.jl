"""
    AbstractCanopyAirProfileModel

How a canopy resolves the in-canopy air-temperature profile from per-layer
leaf sensible-heat sources plus canopy-top and ground boundary values.
Concrete variants implement [`allocate_air_profile`](@ref) and
[`canopy_heat_vapor_transport!`](@ref). Currently: [`KTheoryAirProfile`](@ref)
(steady first-order-closure diffusion).
"""
abstract type AbstractCanopyAirProfileModel end

"""
    allocate_air_profile(model, canopy_height, displacement_height, plant_area_index, heights, n_layers, boundary_layer_model)

Pre-allocate and precompute the structural buffers a model needs.
"""
function allocate_air_profile end

"""
    CanopyAirProfileForcing(; above_canopy_column, wind_speed,
        canopy_top_air_temperature, canopy_top_vapor_density, ground_temperature,
        ground_relative_humidity, sensible_heat_source, evaporation_mass_flow,
        leaf_temperature, atmospheric_pressure)

Forcing and state for one [`canopy_heat_vapor_transport!`](@ref) call: the
[`AboveCanopyColumn`](@ref), the canopy-top and ground boundary values, and the
per-layer leaf sources and temperatures of the current canopy state. Every
air-profile model receives the same forcing and reads the fields it needs.
"""
@kwdef struct CanopyAirProfileForcing{AC,WS,TT,TVD,GT,GRH,SH,EM,LT,P}
    above_canopy_column::AC
    wind_speed::WS
    canopy_top_air_temperature::TT
    canopy_top_vapor_density::TVD
    ground_temperature::GT
    ground_relative_humidity::GRH
    sensible_heat_source::SH
    evaporation_mass_flow::EM
    leaf_temperature::LT
    atmospheric_pressure::P
end

"""
    canopy_heat_vapor_transport!(buffers, model, boundary_layer_model, forcing::CanopyAirProfileForcing;
                              site, vapour_pressure_equation=GoffGratch(), ground_coupling)

The model's physical map from leaf sources and boundary values to the
in-canopy air temperature and vapour density, written into
`buffers.air_temperature` and `buffers.vapor_density`. A model whose map
depends on the current air state reads it from `buffers.air_temperature_prev`
and `buffers.vapor_density_prev`. [`evaluate_canopy!`](@ref) takes its
difference from the current state as the air-profile residual.
`ground_coupling` sets the resistance between the ground and the lowest layer.
"""
function canopy_heat_vapor_transport! end
