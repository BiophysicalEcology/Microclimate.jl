struct AtmosphericProfile{AT,WS,RH,CHF,FV}
    air_temperature::AT        # Matrix (nsteps × nheights)
    wind_speed::WS             # Matrix (nsteps × nheights)
    relative_humidity::RH      # Matrix (nsteps × nheights)
    convective_heat_flux::CHF  # Vector (nsteps)
    friction_velocity::FV      # Vector (nsteps)
end
function AtmosphericProfile(nsteps::Int, nheights::Int)
    AtmosphericProfile(
        Matrix{typeof(1.0u"K")}(undef, nsteps, nheights),
        Matrix{typeof(1.0u"m/s")}(undef, nsteps, nheights),
        Matrix{Float64}(undef, nsteps, nheights),
        Vector{typeof(1.0u"W/m^2")}(undef, nsteps),
        Vector{typeof(1.0u"m/s")}(undef, nsteps),
    )
end

"""
    CanopyOutput

Per-layer and per-hour canopy diagnostics. Zero-column/all-zero for
`NoCanopy` (mirrors `snow_temperature` being `(nsteps, 0)` for `NoSnow`).
"""
struct CanopyOutput{LT,AT,WS,RH,GAS,GAL,CAS,CAL,GTF,IT}
    leaf_temperature::LT              # nsteps × n_canopy_layers Matrix
    air_temperature::AT               # nsteps × n_canopy_layers Matrix (in-canopy, distinct from profile.air_temperature)
    wind_speed::WS                    # nsteps × n_canopy_layers Matrix (in-canopy)
    relative_humidity::RH             # nsteps × n_canopy_layers Matrix (in-canopy)
    ground_absorbed_shortwave::GAS    # nsteps Vector
    ground_absorbed_longwave::GAL     # nsteps Vector
    canopy_absorbed_shortwave::CAS    # nsteps Vector
    canopy_absorbed_longwave::CAL     # nsteps Vector
    ground_throughfall::GTF           # nsteps Vector
    iterations::IT                    # nsteps Vector{Int}, Picard iteration count
end
function CanopyOutput(nsteps::Int, n_canopy_layers::Int)
    CanopyOutput(
        zeros(typeof(1.0u"K"), nsteps, n_canopy_layers),
        zeros(typeof(1.0u"K"), nsteps, n_canopy_layers),
        zeros(typeof(1.0u"m/s"), nsteps, n_canopy_layers),
        zeros(Float64, nsteps, n_canopy_layers),
        zeros(typeof(1.0u"W/m^2"), nsteps),
        zeros(typeof(1.0u"W/m^2"), nsteps),
        zeros(typeof(1.0u"W/m^2"), nsteps),
        zeros(typeof(1.0u"W/m^2"), nsteps),
        zeros(typeof(1.0u"kg/m^2"), nsteps),
        zeros(Int, nsteps),
    )
end

@kwdef struct MicroResult{P,AT,WS,RH,CC,GS,DF,SkT,SoT,SM,SWP,SH,STC,SPH,SBD,SW,SR,Pr,SF,SD,SDN,SNT,CB}
    pressure::P
    reference_temperature::AT
    reference_wind_speed::WS
    reference_humidity::RH
    cloud_cover::CC
    global_radiation::GS
    diffuse_fraction::DF
    sky_temperature::SkT
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
    snow_fall::SF
    snow_depth::SD
    snow_density::SDN
    snow_temperature::SNT  # nsteps × n_snow matrix; size (nsteps, 0) when NoSnow
    canopy::CB
end
function MicroResult(nsteps::Int, num_nodes::Int, nheights::Int, solar_radiation::NamedTuple, n_snow::Int=0, n_canopy_layers::Int=0)
    return MicroResult(;
        pressure = Array{typeof(1.0u"Pa")}(undef, nsteps),
        reference_temperature = Array{typeof(1.0u"K")}(undef, nsteps),
        reference_wind_speed = Array{typeof(1.0u"m/s")}(undef, nsteps),
        reference_humidity = Array{Float64}(undef, nsteps),
        cloud_cover = Array{Float64}(undef, nsteps),
        global_radiation = Array{typeof(1.0u"W/m^2")}(undef, nsteps),
        diffuse_fraction = Array{Float64}(undef, nsteps),
        sky_temperature = Array{typeof(1.0u"K")}(undef, nsteps),
        soil_temperature = zeros(typeof(1.0u"K"), nsteps, num_nodes),
        soil_moisture = zeros(Float64, nsteps, num_nodes),
        soil_water_potential = zeros(typeof(1.0u"J/kg"), nsteps, num_nodes),
        soil_humidity = zeros(Float64, nsteps, num_nodes),
        soil_thermal_conductivity = Array{typeof(1.0u"W/m/K")}(undef, nsteps, num_nodes),
        soil_heat_capacity = Array{typeof(1.0u"J/kg/K")}(undef, nsteps, num_nodes),
        soil_bulk_density = Array{typeof(1.0u"kg/m^3")}(undef, nsteps, num_nodes),
        surface_water = Array{typeof(1.0u"kg/m^2")}(undef, nsteps),
        solar_radiation = solar_radiation,
        profile = AtmosphericProfile(nsteps, nheights),
        snow_fall = zeros(typeof(1.0u"cm/hr"), nsteps),
        snow_depth = zeros(typeof(1.0u"cm"), nsteps),
        snow_density = zeros(typeof(1.0u"g/cm^3"), nsteps),
        snow_temperature = zeros(typeof(1.0u"K"), nsteps, n_snow),
        canopy = CanopyOutput(nsteps, n_canopy_layers),
    )
end

Base.show(io::IO, ::MicroResult) = print(io, "MicroResult")
