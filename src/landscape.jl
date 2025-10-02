Base.@kwdef struct SoilLayers
    depp::Vector{typeof(0.0u"m")}
    wc::Vector{typeof(1.0u"J/K/m^2")}
    c::Vector{typeof(1.0u"W/K/m^2")}
end

Base.@kwdef struct MoistLayers
    P::Vector{typeof(1.0u"J/kg")} 
    Z::Vector{typeof(1.0u"m")}
    V::Vector{typeof(1.0u"kg/m^2")}
    W::Vector{typeof(1.0u"m^3/m^3")}
    WN::Vector{typeof(1.0u"m^3/m^3")}
    K::Vector{typeof(1.0u"kg*s/m^3")}
    H::Vector{Float64}
    T::Vector{typeof(1.0u"K")}
    rh_soil::Vector{Float64}
    ψ_soil::Vector{typeof(1.0u"J/kg")}
    ψ_root::Vector{typeof(1.0u"J/kg")}
    PR::Vector{typeof(1.0u"J/kg")}
    PP::Vector{typeof(1.0u"J/kg")}
    B1::Vector{Float64}
    N::Vector{Float64}
    N1::Vector{Float64}
    WS::Vector{Float64}
    RR::Vector{typeof(1.0u"m^4/kg/s")}
    BZ::Vector{typeof(1.0u"m")}
    JV::Vector{typeof(1.0u"kg/m^2/s")}
    DJ::Vector{typeof(1.0u"kg*s/m^4")}
    CP::Vector{typeof(1.0u"kg*s/m^4")}
    A::Vector{typeof(1.0u"kg*s/m^4")}
    B::Vector{typeof(1.0u"kg*s/m^4")}
    C::Vector{typeof(1.0u"kg*s/m^4")}
    C2::Vector{Float64}
    F::Vector{typeof(1.0u"kg/m^2/s")}
    F2::Vector{typeof(1.0u"J/kg")}
    DP::Vector{typeof(1.0u"J/kg")}
end


Base.@kwdef struct MicroParams{SP,D<:Vector{<:Number},H<:Vector{<:Number},ReH,RoH<:Number,Sl<:Number,E<:Number,P<:Number,TD<:Number}
    soilprops::SP
    depths::D
    heights::H
    reference_height::ReH
    roughness_height::RoH
    κ::Float64
    slope::Sl
    shade::Float64
    viewfactor::Float64
    elevation::E
    P_atmos::P
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
    buffers::B = (;)
end

#Base.@kwdef struct MicroInputs{MP::MicroParams,MF<:MicroForcing,SL<:SoilLayers}
#    params::MP
#    forcing::MF
#    soillayers::SL
#end

function init_soillayers(N)
    depp = fill(0.0u"cm", N + 1)
    wc = fill(1.0u"J/K/m^2", N)
    c = fill(1.0u"W/K/m^2", N)
    return SoilLayers(depp, wc, c)
end

function init_moistlayers(M)
    MoistLayers(
        P = zeros(typeof(1.0u"J/kg"), M+1),
        Z = zeros(typeof(1.0u"m"), M+1),
        V = zeros(typeof(1.0u"kg/m^2"), M+1),
        W = zeros(typeof(1.0u"m^3/m^3"), M+1),
        WN = zeros(typeof(1.0u"m^3/m^3"), M+1),
        K = zeros(typeof(1.0u"kg*s/m^3"), M+1),
        H = zeros(M+1),
        T = zeros(typeof(1.0u"K"), M+1),
        rh_soil = zeros(M),
        ψ_soil = zeros(typeof(1.0u"J/kg"), M),
        ψ_root = zeros(typeof(1.0u"J/kg"), M),
        PR = zeros(typeof(1.0u"J/kg"), M+1),
        PP = zeros(typeof(1.0u"J/kg"), M+1),
        B1 = zeros(M+1),
        N = zeros(M+1),
        N1 = zeros(M+1),
        WS = zeros(M+1),
        RR = zeros(typeof(1.0u"m^4/kg/s"), M+1),
        BZ = zeros(typeof(1.0u"m"), M+1),
        JV = zeros(typeof(1.0u"kg/m^2/s"), M+1),
        DJ = zeros(typeof(1.0u"kg*s/m^4"), M+1),
        CP = zeros(typeof(1.0u"kg*s/m^4"), M+1),
        A = zeros(typeof(1.0u"kg*s/m^4"), M+1),
        B = zeros(typeof(1.0u"kg*s/m^4"), M+1),
        C = zeros(typeof(1.0u"kg*s/m^4"), M+1),
        C2 = zeros(M+1),
        F = zeros(typeof(1.0u"kg/m^2/s"), M+1),
        F2 = zeros(typeof(1.0u"J/kg"), M+1),
        DP = zeros(typeof(1.0u"J/kg"), M+1)
    )
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
