@kwdef struct SoilEnergyInputs{F,B,SP,D<:Vector{<:Number},H<:Vector{<:Number},ST,MT,EI,SW}
    forcing::F
    buffers::B
    soil_thermal_model::SP
    depths::D
    heights::H
    nodes::Vector{Float64}
    solar_terrain::ST #TODO make just one terrain
    micro_terrain::MT
    environment_instant::EI
    runmoist::Bool
    soil_wetness::SW
end

@kwdef struct MicroForcing{
    S<:AbstractInterpolation,ZE<:AbstractInterpolation,ZS<:AbstractInterpolation,T<:AbstractInterpolation,
    V<:AbstractInterpolation,RH<:AbstractInterpolation,CL<:AbstractInterpolation,P<:AbstractInterpolation,
}
    interpolate_solar::S
    interpolate_zenith::ZE
    interpolate_slope_zenith::ZS
    interpolate_temperature::T
    interpolate_wind::V
    interpolate_humidity::RH
    interpolate_cloud::CL
    interpolate_pressure::P

end

abstract type AbstractEnvironment end

@kwdef struct MicroResult{P,AT,WS,RH,CC,GS,DrS,DfS,ZA,SkT,SoT,SM,SWP,SH,STC,SPH,SBD,SW,SR,Pr} <: AbstractEnvironment
    pressure::P 
    reference_temperature::AT 
    reference_wind_speed::WS
    reference_humidity::RH
    cloud_cover::CC
    solar_radiation::GS
    direct_solar::DrS
    diffuse_solar::DfS
    zenith_angle::ZA 
    sky_temperature::SkT
    # TODO: should things like soil_temperature be sub-components? soil.temperature ?
    soil_temperature::SoT
    soil_moisture::SM
    soil_water_potential::SWP
    soil_humidity::SH 
    soil_thermal_conductivity::STC
    soil_heat_capacity::SPH 
    soil_bulk_density::SBD 
    surface_water::SW
    solrad::SR
    profile::Pr
end
function MicroResult(nsteps::Int, numnodes_a::Int)
    return MicroResult(;
        pressure = Array{typeof(1.0u"Pa")}(undef, nsteps),
        reference_temperature = Array{typeof(1.0u"K")}(undef, nsteps),
        reference_wind_speed = Array{typeof(1.0u"m/s")}(undef, nsteps),
        reference_humidity = Array{Float64}(undef, nsteps),
        cloud_cover = Array{Float64}(undef, nsteps),
        solar_radiation = Array{typeof(1.0u"W/m^2")}(undef, nsteps),
        direct_solar = Array{typeof(1.0u"W/m^2")}(undef, nsteps),
        diffuse_solar = Array{typeof(1.0u"W/m^2")}(undef, nsteps),
        zenith_angle = Array{typeof(1.0u"°")}(undef, nsteps),
        sky_temperature = Array{typeof(1.0u"K")}(undef, nsteps),
        soil_temperature = Array{typeof(1.0u"K")}(undef, nsteps, numnodes_a),
        soil_moisture = Array{Float64}(undef, nsteps, numnodes_a),
        soil_water_potential = Array{typeof(1.0u"J/kg")}(undef, nsteps, numnodes_a),
        soil_humidity = Array{Float64}(undef, nsteps, numnodes_a),
        soil_thermal_conductivity = Array{typeof(1.0u"W/m/K")}(undef, nsteps, numnodes_a),
        soil_heat_capacity = Array{typeof(1.0u"J/kg/K")}(undef, nsteps, numnodes_a),
        soil_bulk_density = Array{typeof(1.0u"kg/m^3")}(undef, nsteps, numnodes_a),
        surface_water = Array{typeof(1.0u"kg/m^2")}(undef, nsteps),
        solrad = nothing,
        profile = Array{Any}(undef, nsteps),
    )
end

Base.show(io::IO, mr::MicroResult) = print(io, "MicroResult")


abstract type AbstractSoilThermalModel end

# TODO are these parameters for a specific named model
@kwdef struct CampbelldeVriesSoilThermal <: AbstractSoilThermalModel 
    deVries_shape_factor
    mineral_conductivity
    mineral_density
    mineral_heat_capacity
    bulk_density
    saturation_moisture
    recirculation_power
    return_flow_threshold
end

abstract type AbstractSoilMoistureModel end

# TODO whos model is this what is it called
@kwdef struct SoilMoistureModel <: AbstractSoilMoistureModel
    air_entry_water_potential
    saturated_hydraulic_conductivity
    Campbells_b_parameter
    soil_bulk_density2
    soil_mineral_density2
    root_density
    root_resistance
    stomatal_closure_potential
    leaf_resistance
    stomatal_stability_parameter
    root_radius
    moist_error
    moist_count
    moist_step
    maxpool
end

abstract type AbstractTerrain end

# TODO is there a more specific name for this collection of terrain variables
@kwdef struct MicroTerrain <: AbstractTerrain
    elevation
    roughness_height = nothing
    karman_constant = nothing
    dyer_constant = nothing
    viewfactor = nothing
end

# TODO: this should be more generic.
# We could possible make a field type that is either interpolated or indexed
# so we just mix min-max fields with e.g. daily fields in a single environment object
@kwdef struct MonthlyMinMaxEnvironment{AT,W,H,C,M}# <: AbstractEnvironment
    reference_temperature_min::AT
    reference_temperature_max::AT
    reference_wind_min::W
    reference_wind_max::W
    reference_humidity_min::H
    reference_humidity_max::H
    cloud_min::C
    cloud_max::C
    minima_times::M
    maxima_times::M
end
@kwdef struct DailyTimeseries <: AbstractEnvironment
    shade
    soil_wetness
    surface_emissivity
    cloud_emissivity
    rainfall
    deep_soil_temperature
    leaf_area_index
end
@kwdef struct HourlyTimeseries <: AbstractEnvironment
    pressure
    reference_temperature
    reference_humidity
    reference_wind_speed
    solar_radiation
    longwave_radiation
    cloud_cover
    rainfall
    zenith_angle
end