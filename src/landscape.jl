@kwdef struct SoilEnergyInputs{F,B,SP,D<:Vector{<:Number},H<:Vector{<:Number},T,EI,SW}
    forcing::F
    buffers::B
    soil_thermal_model::SP
    depths::D
    heights::H
    nodes::Vector{Float64}
    terrain::T
    environment_instant::EI
    runmoist::Bool
    soil_wetness::SW
end

@kwdef struct MicroForcing{
    S<:AbstractInterpolation,ZE<:AbstractInterpolation,ZS<:AbstractInterpolation,T<:AbstractInterpolation,
    V<:AbstractInterpolation,RH<:AbstractInterpolation,CL<:AbstractInterpolation,
}
    interpolate_solar::S
    interpolate_zenith::ZE
    interpolate_slope_zenith::ZS
    interpolate_temperature::T
    interpolate_wind::V
    interpolate_humidity::RH
    interpolate_cloud::CL
end

abstract type AbstractEnvironment end

@kwdef struct MicroResult{AT,WS,RH,CC,GS,DrS,DfS,ZA,SkT,SoT,SM,SWP,SH,STC,SPH,SBD,SW,SR,Pr} <: AbstractEnvironment
    reference_temperature::AT 
    reference_wind_speed::WS
    reference_humidity::RH
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
    soil_heat_capacity::SPH 
    soil_bulk_density::SBD 
    surface_water::SW
    solrad::SR
    profile::Pr
end
function MicroResult(nsteps::Int, numnodes_a::Int)
    return MicroResult(;
        reference_temperature = Array{typeof(1.0u"K")}(undef, nsteps),
        reference_wind_speed = Array{typeof(1.0u"m/s")}(undef, nsteps),
        reference_humidity = Array{Float64}(undef, nsteps),
        cloud_cover = Array{Float64}(undef, nsteps),
        global_solar = Array{typeof(1.0u"W/m^2")}(undef, nsteps),
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
@kwdef struct Terrain <: AbstractTerrain
    elevation
    horizon_angles
    slope
    aspect
    # TODO these are not needed in solrad
    roughness_height = nothing
    karman_constant = nothing
    P_atmos = atmospheric_pressure(elevation)
    viewfactor = 1 - sum(sin.(horizon_angles)) / length(horizon_angles) # convert horizon angles to radians and calc view factor(s)
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
    albedo
    shade
    soil_wetness
    surface_emissivity
    cloud_emissivity
    rainfall
    deep_soil_temperature
    leaf_area_index
end
@kwdef struct HourlyTimeseries <: AbstractEnvironment
    reference_temperature
    reference_humidity
    reference_wind_speed
    solar_radiation
    longwave_radiation
    cloud_cover
    rainfall
    zenith_angle
end


abstract type AbstractSolarRadiation end

"""
    SolarRadiation

# TODO who wrote this model what is it called

# Keyword Arguments

- `cmH2O::Real=1`: Precipitable water in cm for atmospheric column (e.g. 0.1: dry, 1.0: moist, 2.0: humid).
- `ϵ::Real=0.0167238`: Orbital eccentricity of Earth.
- `ω::Real=2π/365`: Mean angular orbital velocity of Earth (radians/day).
- `se::Real=0.39779`: Precomputed solar elevation constant.
- `d0::Real=80`: Reference day for declination calculations.
- `iuv::Bool=false`: If `true`, uses the full gamma-function model for diffuse radiation (expensive).
- `scattered::Bool=true`: If `true`, disables scattered light computations (faster).
- `amr::Quantity=25.0u"km"`: Mixing ratio height of the atmosphere.
- `nmax::Integer=111`: Maximum number of wavelength intervals.
- `Iλ::Vector{Quantity}`: Vector of wavelength bins (e.g. in `nm`).
- `OZ::Matrix{Float64}`: Ozone column depth table indexed by latitude band and month (size 19×12).
- `τR`, `τO`, `τA`, `τW`: Vectors of optical depths per wavelength for Rayleigh scattering, ozone, aerosols, and water vapor.
- `Sλ::Vector{Quantity}`: Solar spectral irradiance per wavelength bin (e.g. in `mW * cm^-2 * nm^-1`).
- `FD`, `FDQ`: Radiation scattered from the direct solar beam and reflected radiation
    rescattered downward as a function of wavelength, from tables in Dave & Furukawa (1966).
- `s̄`: a function of τR linked to molecular scattering in the UV range (< 360 nm)
"""
@kwdef struct SolarRadiation <: AbstractSolarRadiation
    solar_geometry_model = McCulloughPorterSolarGeometry()
    cmH2O = 1 # precipitable cm H2O in air column 0.1 = very dry; 1 = moist air conditions; 2 = humid tropical conditions (note this is for the whole atmospheric profile not just near the ground)
    iuv = false # if `true` uses the full gamma-function model for diffuse radiation (expensive)
    scattered = true # if `false` disables scattered light computations (faster)
    amr = 25.0u"km" # mixing ratio height of the atmosphere
    nmax = 111 # Maximum number of wavelength intervals
    # TODO better field names
    Iλ = DEFAULT_Iλ # cector of wavelength bins (e.g. in `nm`)
    OZ = DEFAULT_OZ # ozone column depth table indexed by latitude band and month (size 19×12)
    τR = DEFAULT_τR # vector of optical depths per wavelength for Rayleigh scattering
    τO = DEFAULT_τO # vector of optical depths per wavelength for ozone
    τA = DEFAULT_τA # vector of optical depths per wavelength for aerosols
    τW = DEFAULT_τW # vector of optical depths per wavelength for water vapor
    Sλ = DEFAULT_Sλ # solar spectral irradiance per wavelength bin (e.g. in `mW * cm^-2 * nm^-1`)
    FD = DEFAULT_FD # interpolated function of radiation scattered from the direct solar beam
    FDQ = DEFAULT_FDQ # interpolated function of radiation scattered from ground-reflected radiation
    s̄ = DEFAULT_s̄ # a function of τR linked to molecular scattering in the UV range (< 360 nm)
end
