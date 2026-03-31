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

# ── Soil moisture mode (defined early so CampbellSoilHydraulics can reference it) ──

abstract type AbstractSoilMoistureMode end

"""
    PrescribedSoilMoisture(; precomputed_soil_moisture=nothing)

Use prescribed soil moisture values rather than running a dynamic soil moisture solver.
Soil wetness is taken from the daily environment timeseries.

When `precomputed_soil_moisture` is provided as a matrix (ndepths × ndays), those values
override `initial_soil_moisture` at the start of each day (monthly-representative mode only).
"""
struct PrescribedSoilMoisture{PSM} <: AbstractSoilMoistureMode
    precomputed_soil_moisture::PSM
end
PrescribedSoilMoisture(; precomputed_soil_moisture=nothing) =
    PrescribedSoilMoisture(precomputed_soil_moisture)

"""
    DynamicSoilMoisture()

Run the Campbell soil water balance solver to dynamically compute soil moisture
at each hourly timestep. Tracks the evolving surface soil wetness fraction
internally.
"""
mutable struct DynamicSoilMoisture <: AbstractSoilMoistureMode
    soil_wetness::Float64
end
DynamicSoilMoisture() = DynamicSoilMoisture(0.0)

"""
    CampbellSoilHydraulics(; ..., mode=PrescribedSoilMoisture())

Soil hydraulic parameters for Campbell's (1985) soil water balance model.
The `mode` field controls how soil moisture is handled during simulation:

- `PrescribedSoilMoisture()` — use prescribed wetness from environment data (default)
- `DynamicSoilMoisture()` — run the full Campbell soil water balance solver

# References
Campbell, G. S. (1985). Soil Physics with BASIC. Elsevier.
"""
@kwdef struct CampbellSoilHydraulics{AEWP,SHC,CBP,SBD,SMD,RDen,RRes,SCP,LRes,SSP,RRad,ME,MC,MS,MP,Mode<:AbstractSoilMoistureMode} <: AbstractSoilMoistureModel
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
    mode::Mode = PrescribedSoilMoisture()
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

# Deprecation alias for renamed SoilMoistureModel
const SoilMoistureModel = CampbellSoilHydraulics

# ── Soil moisture mode dispatch functions ─────────────────────────────────

get_soil_wetness(::PrescribedSoilMoisture, environment_instant) = environment_instant.soil_wetness
get_soil_wetness(mode::DynamicSoilMoisture, _) = mode.soil_wetness

init_soil_wetness!(::PrescribedSoilMoisture) = nothing
init_soil_wetness!(mode::DynamicSoilMoisture) = (mode.soil_wetness = 0.0; nothing)

function initialise_soil_humidity!(::DynamicSoilMoisture, output, soil_water_potential, T0)
    water_molar_mass = 0.01801528u"kg/mol"
    output.soil_humidity[1, :] = clamp.(exp.(water_molar_mass .* soil_water_potential ./ (R .* T0)), 0, 1)
end
initialise_soil_humidity!(::PrescribedSoilMoisture, output, soil_water_potential, T0) = nothing

function reset_day_soil_moisture!(mode::PrescribedSoilMoisture, soil_moisture, initial_soil_moisture, day_index)
    if !isnothing(mode.precomputed_soil_moisture)
        soil_moisture .= mode.precomputed_soil_moisture[:, day_index]
    else
        soil_moisture .= initial_soil_moisture
    end
end
function reset_day_soil_moisture!(::DynamicSoilMoisture, soil_moisture, initial_soil_moisture, day_index)
    soil_moisture .= initial_soil_moisture
end

# ── Soil temperature convergence ────────────────────────────────────────────

abstract type AbstractSoilTemperatureConvergence end

"""
    FixedSoilTemperatureIterations(iterations_per_day::Int)

Run a fixed number of iterations of the soil temperature solver per simulated day.
"""
struct FixedSoilTemperatureIterations <: AbstractSoilTemperatureConvergence
    iterations_per_day::Int
end

"""
    SoilTemperatureConvergenceTolerance(; tolerance, max_iterations_per_day)

Iterate the soil temperature solver until the maximum nodal temperature change
between successive full-day passes falls below `tolerance`, or until
`max_iterations_per_day` passes have been completed.
"""
@kwdef struct SoilTemperatureConvergenceTolerance{T} <: AbstractSoilTemperatureConvergence
    tolerance::T
    max_iterations_per_day::Int
end

max_iterations(c::FixedSoilTemperatureIterations) = c.iterations_per_day
max_iterations(c::SoilTemperatureConvergenceTolerance) = c.max_iterations_per_day

may_iterate(c::FixedSoilTemperatureIterations) = c.iterations_per_day > 1
may_iterate(::SoilTemperatureConvergenceTolerance) = true

function is_converged(::FixedSoilTemperatureIterations, iter, niter, T0, T0_prev)
    return iter >= niter
end
function is_converged(c::SoilTemperatureConvergenceTolerance, iter, niter, T0, T0_prev)
    iter >= niter && return true
    iter <= 1 && return niter <= 1
    tol = ustrip(u"K", c.tolerance)
    max_change = maximum(abs.(ustrip.(u"K", T0) .- ustrip.(u"K", T0_prev)))
    return max_change < tol
end

# ── Time mode ───────────────────────────────────────────────────────────────

abstract type AbstractTimeMode end

"""
    MonthlyRepresentativeMode()

Each simulated day is treated independently as a representative day for its month.
Soil state is reset at the start of each day, and the temperature solver iterates
multiple times per day (controlled by the convergence strategy).
"""
struct MonthlyRepresentativeMode <: AbstractTimeMode end

"""
    HourlyMode(; spinup_first_day=false)

Consecutive-day simulation where soil state carries over between days.
Only one iteration per day is performed (the soil state from the previous day
provides the initial condition). When `spinup_first_day=true`, the first day
is iterated multiple times to establish a quasi-equilibrium initial state.
"""
struct HourlyMode <: AbstractTimeMode
    spinup_first_day::Bool
end
HourlyMode(; spinup_first_day=false) = HourlyMode(spinup_first_day)

independent_days(::MonthlyRepresentativeMode) = true
independent_days(::HourlyMode) = false

function iterations_for_day(::MonthlyRepresentativeMode, convergence::AbstractSoilTemperatureConvergence, day_index)
    return max_iterations(convergence)
end
function iterations_for_day(mode::HourlyMode, convergence::AbstractSoilTemperatureConvergence, day_index)
    return (mode.spinup_first_day && day_index == 1) ? max_iterations(convergence) : 1
end

# ── Diffuse fraction model ──────────────────────────────────────────────────

abstract type AbstractDiffuseFractionModel end

"""
    ErbsDiffuseFraction()

Estimate the diffuse fraction of global solar radiation from the clearness index
using the piecewise model of Erbs et al. (1982).

# References
Erbs, D. G., Klein, S. A., & Duffie, J. A. (1982). Estimation of the diffuse
radiation fraction for hourly, daily and monthly-average global radiation.
*Solar Energy*, 28(4), 293–302.
"""
struct ErbsDiffuseFraction <: AbstractDiffuseFractionModel end

function calc_diffuse_fraction(::ErbsDiffuseFraction, clearness_index)
    if clearness_index <= 0.22
        return 1 - 0.09 * clearness_index
    elseif clearness_index <= 0.80
        return 0.9511 - 0.1604 * clearness_index + 4.388 * clearness_index^2 -
               16.638 * clearness_index^3 + 12.336 * clearness_index^4
    else
        return 0.165
    end
end
