"""
    SoilProfile(; bulk_density, mineral_density, mineral_conductivity, mineral_heat_capacity, hydraulics)

Per-depth soil column profile read by the soil properties model (thermal)
and the soil hydraulic model. Sized to match `MicroModel.soil.depths`.

- `bulk_density` — dry soil bulk density (kg/m³) at each depth node
- `mineral_density` — soil mineral density (kg/m³) at each depth node
- `mineral_conductivity` — soil mineral thermal conductivity (W/m/K) at each depth node
- `mineral_heat_capacity` — soil mineral specific heat (J/kg/K) at each depth node
- `hydraulics` — per-depth hydraulic parameters (e.g. `CampbellHydraulicProfile`)
  matched to the chosen `soil_hydraulic_model`
"""
@kwdef struct SoilProfile{BD,MD,MC,MHC,H}
    bulk_density::BD
    mineral_density::MD
    mineral_conductivity::MC
    mineral_heat_capacity::MHC
    hydraulics::H
end

function example_soil_profile(depths=DEFAULT_DEPTHS;
    bulk_density = 1.3u"Mg/m^3",
    mineral_density = 2.560u"Mg/m^3",
    mineral_conductivity = 1.25u"W/m/K",
    mineral_heat_capacity = 870.0u"J/kg/K",
    hydraulics = example_campbell_hydraulic_profile(depths),
)
    SoilProfile(;
        bulk_density = bulk_density isa AbstractVector ? bulk_density : fill(bulk_density, length(depths)),
        mineral_density = mineral_density isa AbstractVector ? mineral_density : fill(mineral_density, length(depths)),
        mineral_conductivity = mineral_conductivity isa AbstractVector ? mineral_conductivity : fill(mineral_conductivity, length(depths)),
        mineral_heat_capacity = mineral_heat_capacity isa AbstractVector ? mineral_heat_capacity : fill(mineral_heat_capacity, length(depths)),
        hydraulics,
    )
end

"""
    RadiationModel(; solar_radiation_model=SolarProblem(),
                     longwave_model=ViewFactorLongwave(),
                     shortwave_model=AngstromMaxwellShortwave())

Bundle of radiation models used by the simulation:

- `solar_radiation_model::SolarProblem` — clear-sky direct/diffuse/global
  irradiance and solar geometry, computed once up front for every (day, hour)
- `shortwave_model::AbstractShortwaveModel` — cloud-adjusted shortwave
  budget that consumes the clear-sky output
- `longwave_model::AbstractLongwaveModel` — surface longwave budget
  evaluated inside the soil energy ODE
"""
@kwdef struct RadiationModel{SRM,LWM,SWM}
    solar_radiation_model::SRM = SolarProblem()
    longwave_model::LWM = ViewFactorLongwave()
    shortwave_model::SWM = AngstromMaxwellShortwave()
end

"""
    SoilModel(; depths=DEFAULT_DEPTHS, profile, properties_model, hydraulic_model,
                energy_model=SoilHeatTransport1D(), moisture_strategy=PrescribedSoilMoisture())

The soil column:

- `depths` — soil nodes; keep spacing close near the surface
- `profile::SoilProfile` — per-depth `bulk_density`, `mineral_density`,
  `mineral_conductivity`, `mineral_heat_capacity` and hydraulic profile
- `properties_model::AbstractSoilProperties` (e.g. `CampbelldeVriesSoilProperties(...)`)
- `hydraulic_model::AbstractSoilHydraulicsModel` (e.g. `CampbellSoilHydraulics(...)`)
- `energy_model::SoilHeatTransportModel` — soil column energy ODE (carries the
  phase-transition `freezing_model` and ODE solver settings)
- `moisture_strategy` — `PrescribedSoilMoisture()` or `DynamicSoilMoisture(; ...)`
"""
@kwdef struct SoilModel{Dep,SP<:SoilProfile,SPM,SHM,SEM,SMS}
    depths::Dep = DEFAULT_DEPTHS
    profile::SP
    properties_model::SPM
    hydraulic_model::SHM
    energy_model::SEM = SoilHeatTransport1D()
    moisture_strategy::SMS = PrescribedSoilMoisture()
end

"""
    SurfaceModel(; evaporation_model=BulkTransferEvaporation(),
                   condensation_model=GarrattSegalCondensation(),
                   max_surface_pool=1.0e4u"kg/m^2")

Exchange of water at the ground surface:

- `evaporation_model` — surface latent flux
- `condensation_model::AbstractCondensationModel` — ground dew/frost formation
  (`GarrattSegalCondensation()` default, `BulkTransferCondensation()`, or
  `NoCondensation()` to switch it off)
- `max_surface_pool` — maximum water held on the soil surface; water above it is removed as `runoff`
"""
@kwdef struct SurfaceModel{EVM,CDM,MSP}
    evaporation_model::EVM = BulkTransferEvaporation()
    condensation_model::CDM = GarrattSegalCondensation()
    max_surface_pool::MSP = 1.0e4u"kg/m^2"
end

"""
    AtmosphereModel(; heights=[0.01, 2]u"m", radiation=RadiationModel(),
                      boundary_layer_model=MoninObukhov())

The air above the ground:

- `heights` — air nodes for the temperature, wind speed and humidity profile;
  the last is the reference height of the weather data
- `radiation::RadiationModel` — solar, shortwave and longwave models
- `boundary_layer_model` — surface-layer profile and fluxes (`MoninObukhov()`)
"""
@kwdef struct AtmosphereModel{Ht,RAD,BLM}
    heights::Ht = [0.01, 2]u"m"
    radiation::RAD = RadiationModel()
    boundary_layer_model::BLM = MoninObukhov()
end

"""
    MicroConfig(; ...)

Solver iteration and data-delivery strategy. Lives on `MicroModel.config`.

- `convergence`: `FixedIterationConvergence(3)` or `IterationToleranceConvergence(; tolerance, max_iterations_per_day)`
- `rainfall_schedule`: `DailyRainfall()` (default) or `HourlyRainfall()`
- `canopy_soil_convergence`: per-hour convergence between the canopy solve
  and the soil-heat ODE, iterating the two to a jointly self-consistent state
  within the hour by default (`IterationToleranceConvergence(; tolerance=0.05u"K",
  max_iterations_per_day=80)`); `FixedIterationConvergence(1)` for a single
  pass using the previous hour's soil temperature instead.
- `canopy_soil_relaxation`: under-relaxation on each `canopy_soil_convergence`
  pass's `ground_temperature` update (`x_new = relaxation*x_solved +
  (1-relaxation)*x_prev`), `1.0` (no relaxation) by default.
"""
@kwdef struct MicroConfig{CV,RFS,CSC,CSR}
    convergence::CV = FixedIterationConvergence(3)
    rainfall_schedule::RFS = DailyRainfall()
    canopy_soil_convergence::CSC = IterationToleranceConvergence(; tolerance=0.05u"K", max_iterations_per_day=80)
    canopy_soil_relaxation::CSR = 1.0
end

"""
    MicroModel(; hours=DEFAULT_HOURS, site, soil, surface=SurfaceModel(),
                 atmosphere=AtmosphereModel(), snow=NoSnow(), canopy=NoCanopy(),
                 vapour_pressure_equation=GoffGratch(), config=MicroConfig())

Every parameter of the simulation, as one tree:

- `hours` — hours of the day to solve # TODO generalise to sub-hourly
- `site::Site` — properties of the place
- `soil::SoilModel` — the soil column
- `surface::SurfaceModel` — water exchange at the ground surface
- `atmosphere::AtmosphereModel` — air nodes, radiation and boundary layer
- `snow::AbstractSnowModel` — `NoSnow()` (default) or `SnowModel(...)`
- `canopy::AbstractCanopyModel` — `NoCanopy()` (default, scalar `shade`
  vegetation handling) or `MultilayerCanopy(...)` (layer-resolved two-stream
  radiative transfer through the canopy)
- `vapour_pressure_equation` — `GoffGratch()`, `Teten()` or `Huang()`
- `config::MicroConfig` — iteration and data-delivery strategy

Combine with a `MicroInputs` via `MicroProblem(model, inputs; days)` to run.
"""
@kwdef struct MicroModel{H,S,SO<:SoilModel,SU<:SurfaceModel,A<:AtmosphereModel,SNM,CAN,VPE,C}
    hours::H = DEFAULT_HOURS
    site::S
    soil::SO
    surface::SU = SurfaceModel()
    atmosphere::A = AtmosphereModel()
    snow::SNM = NoSnow()
    canopy::CAN = NoCanopy()
    vapour_pressure_equation::VPE = GoffGratch()
    config::C = MicroConfig()
end

"""
    MicroInputs(; environment_minmax, environment_daily, environment_hourly=nothing,
                  initial_soil_temperature=nothing, initial_soil_moisture,
                  initial_snow_depth=0.0u"cm", initial_snow_temperature=u"K"(0.0u"°C"),
                  initial_snow_density=nothing)

Per-run forcing and initial conditions.

- `environment_minmax`, `environment_daily`, `environment_hourly` — input forcings
- `initial_*` — initial conditions (`initial_soil_temperature=nothing` falls back to
  the day-mean reference air temperature)
- `initial_snow_density=nothing` → use `model.snow.snow_density`

Pass a fresh `MicroInputs` to `reinit!(cache, inputs)` to re-solve the same
model on different data without rebuilding the cache.
"""
@kwdef struct MicroInputs{EMM,EH,ED,IST,ISM,ISD,ISNT,ISND}
    environment_minmax::EMM
    environment_hourly::EH = nothing
    environment_daily::ED
    initial_soil_temperature::IST = nothing
    initial_soil_moisture::ISM
    initial_snow_depth::ISD = 0.0u"cm"
    initial_snow_temperature::ISNT = u"K"(0.0u"°C")
    initial_snow_density::ISND = nothing
end

"""
    MicroProblem(model::MicroModel, inputs::MicroInputs;
                 days=DEFAULT_DAYS, time_mode=NonConsecutiveDayMode())

Pairing of a model description with a set of input data, the days of
year to simulate, and how to iterate through them. Pass to
`init(problem)` to allocate a `MicroCache` or `solve(problem)` to run.

- `days` — day-of-year integers to simulate (e.g. `DEFAULT_DAYS` for monthly
  mid-month days, or `1:365` for a full year of daily forcing). Length must
  match the per-day rows of the supplied `environment_daily`.
- `time_mode` — `NonConsecutiveDayMode()` (independent representative days,
  Fortran monthly behaviour) or `ConsecutiveDayMode(; spinup_first_day=false)`
  (continuous run where state carries day-to-day).
"""
struct MicroProblem{D,TM,M<:MicroModel,I<:MicroInputs}
    days::D
    time_mode::TM
    model::M
    inputs::I
end

MicroProblem(model::MicroModel, inputs::MicroInputs;
    days=DEFAULT_DAYS, time_mode=NonConsecutiveDayMode(),
) = MicroProblem(days, time_mode, model, inputs)

"""
    example_microclimate_problem(; kwargs...)

Build a complete, ready-to-solve [`MicroProblem`](@ref) for Madison, Wisconsin, USA,
from `example_*` defaults throughout. Any field of the model tree, the inputs or
the problem can be overridden by name via `kwargs`.
"""
function example_microclimate_problem(;
    days = DEFAULT_DAYS,
    hours = DEFAULT_HOURS,
    depths = DEFAULT_DEPTHS,
    heights = [0.01, 2]u"m",
    site = example_site(),
    soil_profile = example_soil_profile(depths),
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model = example_soil_hydraulic_model(),
    soil_moisture_strategy = PrescribedSoilMoisture(),
    max_surface_pool = 1.0e4u"kg/m^2",
    snow_model = NoSnow(),
    canopy_model = NoCanopy(),
    environment_minmax = example_monthly_weather(),
    environment_daily = example_daily_environment(days),
    environment_hourly = example_hourly_environment(days, hours; elevation=site.elevation),
    config = MicroConfig(),
    initial_soil_temperature = fill(u"K"(7.741667u"°C"), length(depths)),
    initial_soil_moisture = fill(0.42 * 0.25, length(depths)),
)
    soil = SoilModel(;
        depths, profile=soil_profile, properties_model=soil_properties_model,
        hydraulic_model=soil_hydraulic_model, moisture_strategy=soil_moisture_strategy,
    )
    model = MicroModel(;
        hours, site, soil, surface=SurfaceModel(; max_surface_pool), atmosphere=AtmosphereModel(; heights),
        snow=snow_model, canopy=canopy_model, config,
    )
    inputs = MicroInputs(;
        environment_minmax, environment_daily, environment_hourly,
        initial_soil_temperature, initial_soil_moisture,
    )
    MicroProblem(model, inputs; days)
end

"""
    MicroState

Mutable simulation state that evolves during `solve!`. Lives on `MicroCache.state`.
Holds only the fields that need persistent storage across hour/day iterations
(arrays mutated in place). Immutable snow state lives as a local in `solve_soil!`.
"""
mutable struct MicroState{SM,SP}
    soil_moisture::SM
    ∑phase::SP
    # Snapshot of ∑phase at the start of each day, restored before every spinup
    # iteration of that day. Pre-allocated so the day-loop never calls `copy(∑phase)`.
    ∑phase_day_start::SP
end

"""
    MicroBuffers

Pre-allocated workspace. Every buffer the hot path touches lands here once
in `init`. Lives on `MicroCache.buffers`.
"""
struct MicroBuffers{SO,SOB,P,PB,SEB,SP,PT,SWB,SS,IB,CB}
    solar_out::SO                  # SolarRadiation output (NamedTuple of arrays)
    solar::SOB                     # SolarRadiation internal buffers (NamedTuple)
    condensation_profile::P        # atmospheric profile scratch used by ground condensation calculation
    air_profile::PB                # atmospheric profile scratch used by solve_air!
    soil_energy_balance::SEB
    soil_properties::SP
    phase_transition::PT
    soil_water_balance::SWB
    snow::SS                       # snow buffers (NamedTuple for SnowModel; nothing for NoSnow)
    interpolation::IB              # unused; diel `evaluate!` needs no scratch
    canopy::CB                     # canopy buffers (nested NamedTuple for MultilayerCanopy; nothing for NoCanopy)
end

"""
    MicroCache

Pre-allocated workspace for solving a `MicroProblem` without repeated allocation.
Created by `CommonSolve.init(problem::MicroProblem)` and solved in-place with
`CommonSolve.solve!(cache::MicroCache)`.

Use `reinit!(cache, new_inputs::MicroInputs)` to swap the input data while
keeping the model and all pre-allocated arrays, then call `solve!(cache)`
again.

# Example
```julia
cache = init(problem)
output = solve!(cache)

# Run the same model on a different set of inputs
reinit!(cache, new_inputs)
output = solve!(cache)
```
"""
mutable struct MicroCache{MP<:MicroProblem,O<:MicroResult,S<:MicroState,B<:MicroBuffers,F<:Forcing,I,Nsoil}
    problem::MP
    output::O
    state::S
    buffers::B
    forcing::F
    ode_integrator::I
    # Type-level soil-node count, used for `SVector(ntuple(f, Val(N_soil)))` in
    # the day loop so T0 reconstruction is type-stable.
    num_soil_nodes_val::Val{Nsoil}
end
