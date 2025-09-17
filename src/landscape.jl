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
