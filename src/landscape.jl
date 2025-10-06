Base.@kwdef struct SoilEnergyInputs{F,B,SP,D<:Vector{<:Number},H<:Vector{<:Number},ReH,RoH<:Number,Sl<:Number,E<:Number,P<:Number,TD<:Number}
    forcing::F
    buffers::B
    soilprops::SP
    depths::D
    heights::H
    nodes::Vector{Float64}
    terrain::T
    environment_instant::EI
    θ_soil::Vector{Float64}
    runmoist::Bool
end

Base.@kwdef struct MicroForcing{
    S<:AbstractInterpolation,ZE<:AbstractInterpolation,ZS<:AbstractInterpolation,T<:AbstractInterpolation,
    V<:AbstractInterpolation,RH<:AbstractInterpolation,CL<:AbstractInterpolation,
}
    # TODO readable names
    SOLRt::S
    ZENRt::ZE
    ZSLt::ZS
    TAIRt::T
    VELt::V
    RHt::RH
    CLDt::CL
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
