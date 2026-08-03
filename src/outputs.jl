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
struct CanopyOutput{LT,AT,WS,RH,GAS,GAL,CAS,CAL,GTF,IT,DSW,USW,DLW,ULW,AR,NB,CSH,CLH}
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
    boundary_downward_shortwave::DSW  # nsteps × (n_canopy_layers+1) Matrix, canopy top -> ground
    boundary_upward_shortwave::USW    # nsteps × (n_canopy_layers+1) Matrix, canopy top -> ground
    boundary_downward_longwave::DLW   # nsteps × (n_canopy_layers+1) Matrix, canopy top -> ground
    boundary_upward_longwave::ULW     # nsteps × (n_canopy_layers+1) Matrix, canopy top -> ground
    absorbed_radiation::AR            # nsteps × n_canopy_layers Matrix, W/m^2 leaf area (leaf_heat_balance's own input)
    net_balance::NB                   # nsteps × n_canopy_layers Matrix, W, leaf_heat_balance's residual at the converged temperature
    canopy_sensible_heat_flux::CSH    # nsteps Vector, W/m^2 ground area, layers summed
    canopy_latent_heat_flux::CLH      # nsteps Vector, W/m^2 ground area, layers summed
end
function CanopyOutput(nsteps::Int, n_canopy_layers::Int)
    n_boundaries = n_canopy_layers == 0 ? 0 : n_canopy_layers + 1
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
        zeros(typeof(1.0u"W/m^2"), nsteps, n_boundaries),
        zeros(typeof(1.0u"W/m^2"), nsteps, n_boundaries),
        zeros(typeof(1.0u"W/m^2"), nsteps, n_boundaries),
        zeros(typeof(1.0u"W/m^2"), nsteps, n_boundaries),
        zeros(typeof(1.0u"W/m^2"), nsteps, n_canopy_layers),
        zeros(typeof(1.0u"W"), nsteps, n_canopy_layers),
        zeros(typeof(1.0u"W/m^2"), nsteps),
        zeros(typeof(1.0u"W/m^2"), nsteps),
    )
end

@kwdef struct MicroResult{P,AT,WS,RH,CC,GS,DF,SkT,SoT,SM,SWP,SH,STC,SPH,SBD,SW,SR,Pr,SF,SD,SDN,SNT,GHF,CB}
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
    ground_heat_flux::GHF  # nsteps Vector, W/m^2, Fourier's law between the top two soil nodes, downward-positive
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
        ground_heat_flux = zeros(typeof(1.0u"W/m^2"), nsteps),
        canopy = CanopyOutput(nsteps, n_canopy_layers),
    )
end

Base.show(io::IO, ::MicroResult) = print(io, "MicroResult")
