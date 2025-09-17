Base.@kwdef struct SoilLayers
    depp::Vector{typeof(0.0u"m")}
    wc::Vector{typeof(1.0u"J/K/m^2")}
    c::Vector{typeof(1.0u"W/K/m^2")}
end


Base.@kwdef struct MicroParams{SP,D<:Vector{<:Number},ReH<:Number,RoH<:Number,D0,Z<:Number,Sl<:Number,E<:Number,TD}
    soilprops::SP
    depths::D
    reference_height::ReH
    roughness_height::RoH
    d0::D0
    zh::Z
    slope::Sl
    shade::Float64
    viewfactor::Float64
    elevation::E
    albedo::Float64
    sle::Float64
    slep::Float64
    pctwet::Float64
    nodes::Vector{Float64}
    tdeep::TD
    θ_soil::Vector{Float64}
    runmoist::Bool
end

Base.@kwdef struct MicroForcing{
    S<:AbstractInterpolation,ZE<:AbstractInterpolation,ZS<:AbstractInterpolation,T<:AbstractInterpolation,
    V<:AbstractInterpolation,RH<:AbstractInterpolation,CL<:AbstractInterpolation,
}
    SOLRt::S
    ZENRt::ZE
    ZSLt::ZS
    TAIRt::T
    VELt::V
    RHt::RH
    CLDt::CL
end

Base.@kwdef struct MicroInputs{MP<:MicroParams,MF<:MicroForcing,SL<:SoilLayers,B}
    params::MP
    forcing::MF
    soillayers::SL
    buffers::B
end

function init_soillayers(N)
    depp = fill(0.0u"cm", N + 1)
    wc = fill(1.0u"J/K/m^2", N)
    c = fill(1.0u"W/K/m^2", N)
    return SoilLayers(depp, wc, c)
end

abstract type AbstractEnvironment end

@kwdef struct MicroResult{AT,WS,RH,CC,GS,DrS,DfS,ZA,SkT,SoT,SM,SWP,SH,STC,SPH,SBD,SW,SR,Pr} <: AbstractEnvironment
    air_temperature::AT 
    wind_speed::WS
    relative_humidity::RH
    cloud_cover::CC
    global_solar::GS
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
    soil_specific_heat::SPH 
    soil_bulk_density::SBD 
    surface_water::SW
    solrad::SR
    profile::Pr
end

Base.show(io::IO, mr::MicroResult) = print(io, "MicroResult")
