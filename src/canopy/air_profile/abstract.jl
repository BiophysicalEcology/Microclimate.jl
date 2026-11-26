"""
    AbstractCanopyAirProfileModel

How a canopy turns each layer's leaf heat and vapor sources, plus the values at the canopy top
and the ground, into air temperature and vapor density through the canopy. Concrete models
implement [`allocate_air_profile`](@ref) and [`canopy_heat_vapor_transport!`](@ref):
[`RaupachLTheoryAirProfile`](@ref) and [`KTheoryAirProfile`](@ref).
"""
abstract type AbstractCanopyAirProfileModel end

"""
    allocate_air_profile(model, canopy_height, displacement_height, plant_area_index, heights, n_layers, boundary_layer_model)

Allocate the buffers a model needs, and compute what depends only on the canopy's structure.
"""
function allocate_air_profile end

"""
    CanopyAirProfileForcing(; friction_velocity, inverse_obukhov_length, wind_speed,
        canopy_top_air_temperature, canopy_top_vapor_density, ground_temperature,
        ground_relative_humidity, sensible_heat_flux, evaporation_mass_flux, atmospheric_pressure)

Everything [`canopy_heat_vapor_transport!`](@ref) takes from the current canopy state: the
turbulence above the canopy, the values at the canopy top and the ground, and each layer's
leaf sources. Every air-profile model receives the same forcing and reads the fields it needs.
"""
@kwdef struct CanopyAirProfileForcing{FV,OL,WS,TT,TVD,GT,GRH,SH,EM,P}
    friction_velocity::FV
    inverse_obukhov_length::OL
    wind_speed::WS
    canopy_top_air_temperature::TT
    canopy_top_vapor_density::TVD
    ground_temperature::GT
    ground_relative_humidity::GRH
    sensible_heat_flux::SH
    evaporation_mass_flux::EM
    atmospheric_pressure::P
end

"""
    canopy_heat_vapor_transport!(buffers, model, boundary_layer_model, forcing::CanopyAirProfileForcing;
                              site, vapour_pressure_equation=GoffGratch())

Writes the air temperature and vapor density that the leaf sources and boundary values
produce into `buffers.air_temperature` and `buffers.vapor_density`. A model that also depends
on the air's current state reads it from `buffers.current_air_temperature` and
`buffers.current_vapor_density`. [`evaluate_canopy!`](@ref) takes the difference between the
current state and this result as the air-profile residual.
"""
function canopy_heat_vapor_transport! end
