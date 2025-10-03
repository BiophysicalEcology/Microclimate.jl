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

# TODO: these are really buffers
function init_soillayers(N)
    depp = fill(0.0u"cm", N + 1)
    wc = fill(1.0u"J/K/m^2", N)
    c = fill(1.0u"W/K/m^2", N)
    return SoilLayers(depp, wc, c)
end

function init_moistlayers(M)
    MoistLayers(
        P = zeros(M+1) .* u"J/kg",
        Z = zeros(M+1) .* u"m",
        V = zeros(M+1) .* u"kg/m^2",
        W = zeros(M+1) .* u"m^3/m^3",
        WN = zeros(M+1) .* u"m^3/m^3",
        K = zeros(M+1) .* u"kg*s/m^3",
        H = zeros(M+1),
        T = zeros(M+1) .* u"K",
        rh_soil = zeros(M),
        ψ_soil = zeros(M) .* u"J/kg",
        ψ_root = zeros(M) .* u"J/kg",
        PR = zeros(M+1) .* u"J/kg",
        PP = zeros(M+1) .* u"J/kg",
        B1 = zeros(M+1),
        N = zeros(M+1),
        N1 = zeros(M+1),
        WS = zeros(M+1),
        RR = zeros(M+1) .* u"m^4/kg/s",
        BZ = zeros(M+1) .* u"m",
        JV = zeros(M+1) .* u"kg/m^2/s",
        DJ = zeros(M+1) .* u"kg*s/m^4",
        CP = zeros(M+1) .* u"kg*s/m^4",
        A = zeros(M+1) .* u"kg*s/m^4",
        B = zeros(M+1) .* u"kg*s/m^4",
        C = zeros(M+1) .* u"kg*s/m^4",
        C2 = zeros(M+1),
        F = zeros(M+1) .* u"kg/m^2/s",
        F2 = zeros(M+1) .* u"J/kg",
        DP = zeros(M+1) .* u"J/kg"
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
function MicroResult(nsteps::Int, numodes_a::Int)
    return MicroResult(;
        # TODO ... the rest
        air_temperatures,
        wind_speeds,
        relative_humidity,
        cloud_cover = cloud_covers,
        global_solar,
        direct_solar,
        diffuse_solar,
        zenith_angles,
        sky_temperature = Array{typeof(1.0u"K")}(undef, nsteps),
        soil_temperature = Array{SVector{numnodes_a,typeof(1.0u"K")}}(undef, nsteps),
        soil_moisture = Array{Float64}(undef, nsteps, numnodes_a),
        soil_water_potential = Array{typeof(1.0u"J/kg")}(undef, nsteps, numnodes_a),
        soil_humidity = Array{Float64}(undef, nsteps, numnodes_a),
        soil_thermal_conductivity = Array{typeof(1.0u"W/m/K")}(undef, nsteps, numnodes_a),
        soil_specific_heat = Array{typeof(1.0u"J/kg/K")}(undef, nsteps, numnodes_a),
        soil_bulk_density = Array{typeof(1.0u"kg/m^3")}(undef, nsteps, numnodes_a),
        surface_water = Array{typeof(1.0u"kg/m^2")}(undef, nsteps),
    )
    # return MicroResult(;
    #     air_temperatures,
    #     wind_speeds,
    #     relative_humidity,
    #     cloud_cover = cloud_covers,
    #     global_solar,
    #     direct_solar,
    #     diffuse_solar,
    #     zenith_angles,
    #     sky_temperature = T_skys,
    #     soil_temperature = reduce(vcat, transpose.(T_soils)),
    #     soil_moisture = θ_soils,
    #     soil_water_potential = ψ_soils,
    #     soil_humidity = rh_soils,
    #     soil_thermal_conductivity = λ_bulk,
    #     soil_specific_heat = c_p_bulk,
    #     soil_bulk_density = ρ_bulk,
    #     surface_water = pools,
    #     solrad=solrad_out,
    #     profile=profile_out,
    # )
end

Base.show(io::IO, mr::MicroResult) = print(io, "MicroResult")
