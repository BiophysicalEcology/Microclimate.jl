@kwdef struct SoilEnergyInputs{F,B,SP,D<:Vector{<:Number},H<:Vector{<:Number},ST,MT,EI,SW,VP,LW}
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
    vapour_pressure_equation::VP = GoffGratch()
    longwave_sky::LW
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

struct MicroProfile{AT,WS,RH,CHF,FV}
    air_temperature::AT        # Matrix (nsteps × nheights)
    wind_speed::WS             # Matrix (nsteps × nheights)
    relative_humidity::RH      # Matrix (nsteps × nheights)
    convective_heat_flux::CHF  # Vector (nsteps)
    friction_velocity::FV      # Vector (nsteps)
end
function MicroProfile(nsteps::Int, nheights::Int)
    MicroProfile(
        Matrix{typeof(1.0u"K")}(undef, nsteps, nheights),
        Matrix{typeof(1.0u"m/s")}(undef, nsteps, nheights),
        Matrix{Float64}(undef, nsteps, nheights),
        Vector{typeof(1.0u"W/m^2")}(undef, nsteps),
        Vector{typeof(1.0u"m/s")}(undef, nsteps),
    )
end

@kwdef struct MicroResult{P,AT,WS,RH,CC,GS,DF,SkT,SoT,SM,SWP,SH,STC,SPH,SBD,SW,SR,Pr} <: AbstractEnvironment
    pressure::P
    reference_temperature::AT
    reference_wind_speed::WS
    reference_humidity::RH
    cloud_cover::CC
    global_radiation::GS
    diffuse_fraction::DF
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
    solar_radiation::SR
    profile::Pr
end
function MicroResult(nsteps::Int, num_nodes::Int, nheights::Int, solar_radiation::NamedTuple)

    return MicroResult(;
        pressure = Array{typeof(1.0u"Pa")}(undef, nsteps),
        reference_temperature = Array{typeof(1.0u"K")}(undef, nsteps),
        reference_wind_speed = Array{typeof(1.0u"m/s")}(undef, nsteps),
        reference_humidity = Array{Float64}(undef, nsteps),
        cloud_cover = Array{Float64}(undef, nsteps),
        global_radiation = Array{typeof(1.0u"W/m^2")}(undef, nsteps),
        diffuse_fraction = Array{Float64}(undef, nsteps),
        sky_temperature = Array{typeof(1.0u"K")}(undef, nsteps),
        soil_temperature = Array{typeof(1.0u"K")}(undef, nsteps, num_nodes),
        soil_moisture = Array{Float64}(undef, nsteps, num_nodes),
        soil_water_potential = Array{typeof(1.0u"J/kg")}(undef, nsteps, num_nodes),
        soil_humidity = Array{Float64}(undef, nsteps, num_nodes),
        soil_thermal_conductivity = Array{typeof(1.0u"W/m/K")}(undef, nsteps, num_nodes),
        soil_heat_capacity = Array{typeof(1.0u"J/kg/K")}(undef, nsteps, num_nodes),
        soil_bulk_density = Array{typeof(1.0u"kg/m^3")}(undef, nsteps, num_nodes),
        surface_water = Array{typeof(1.0u"kg/m^2")}(undef, nsteps),
        solar_radiation = solar_radiation,
        profile = MicroProfile(nsteps, nheights),
    )
end

Base.show(io::IO, mr::MicroResult) = print(io, "MicroResult")


abstract type AbstractSoilThermalModel end

# TODO are these parameters for a specific named model
    @kwdef struct CampbelldeVriesSoilThermal{SF,MC,MD,MHC,BD,SM,RP,RFT} <: AbstractSoilThermalModel
    de_vries_shape_factor::SF
    mineral_conductivity::MC
    mineral_density::MD
    mineral_heat_capacity::MHC
    bulk_density::BD
    saturation_moisture::SM
    recirculation_power::RP
    return_flow_threshold::RFT
end

abstract type AbstractSoilMoistureModel end

# TODO whos model is this what is it called
@kwdef struct SoilMoistureModel{AEWP,SHC,CBP,SBD,SMD,RDen,RRes,SCP,LRes,SSP,RRad,ME,MC,MS,MP} <: AbstractSoilMoistureModel
    air_entry_water_potential::AEWP
    saturated_hydraulic_conductivity::SHC
    campbell_b_parameter::CBP
    soil_bulk_density2::SBD
    soil_mineral_density2::SMD
    root_density::RDen
    root_resistance::RRes
    stomatal_closure_potential::SCP
    leaf_resistance::LRes
    stomatal_stability_parameter::SSP
    root_radius::RRad
    moist_error::ME
    moist_count::MC
    moist_step::MS
    maxpool::MP
end

abstract type AbstractTerrain end

# TODO is there a more specific name for this collection of terrain variables
@kwdef struct MicroTerrain{E,RH,KC,DC,VF} <: AbstractTerrain
    elevation::E
    roughness_height::RH = nothing
    karman_constant::KC = nothing
    dyer_constant::DC = nothing
    viewfactor::VF = nothing
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
"""
    DailyMinMaxEnvironment

Per-day analogue of `MonthlyMinMaxEnvironment` for consecutive-day simulations
(ERA5, station data, etc.).  Each entry corresponds to one actual calendar day.
Passing this to `simulate_microclimate` automatically sets `daily=true` so that
consecutive days inherit soil state and iterate once.
"""
@kwdef struct DailyMinMaxEnvironment{AT,W,H,C,M}
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
@kwdef struct DailyTimeseries{Sh,SW,SE,CE,R,DST,LAI} <: AbstractEnvironment
    shade::Sh
    soil_wetness::SW
    surface_emissivity::SE
    cloud_emissivity::CE
    rainfall::R
    deep_soil_temperature::DST
    leaf_area_index::LAI
end
@kwdef struct HourlyTimeseries{P,RT,RH,RWS,GR,LW,CC,R,ZA} <: AbstractEnvironment
    pressure::P
    reference_temperature::RT
    reference_humidity::RH
    reference_wind_speed::RWS
    global_radiation::GR
    longwave_radiation::LW
    cloud_cover::CC
    rainfall::R
    zenith_angle::ZA
end
