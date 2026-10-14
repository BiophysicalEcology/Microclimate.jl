"""
    SoilProfile(; bulk_density, mineral_density, mineral_conductivity, mineral_heat_capacity, hydraulics)

Per-depth soil column profile read by the soil properties model (thermal)
and the soil hydraulic model. Sized to match `MicroModel.depths`.

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
    MicroConfig(; ...)

Solver/iteration/data-delivery strategy. Lives on `MicroModel.config`.
Physical-process models live directly on `MicroModel`; this struct is
strictly the "how we iterate and how data is delivered" side of the model.

- `convergence`: `FixedIterationConvergence(3)` or `IterationToleranceConvergence(; tolerance, max_iterations_per_day)`
- `rainfall_schedule`: `DailyRainfall()` (default) or `HourlyRainfall()`
- `soil_moisture_strategy`: `PrescribedSoilMoisture()` or `DynamicSoilMoisture(; ...)`
- `max_surface_pool`: numerical clamp on the surface-pool state variable
  (not a physical limit — keeps the pool integration from running away)
- `max_pond_depth`: physical surface-ponding limit; standing water above this
  depth is tracked as `runoff` rather than left available to infiltrate.
  Defaults to `max_surface_pool`'s own value, so `runoff` stays `0` unless
  set to something physically meaningful for the site.
- `canopy_soil_convergence`: per-hour convergence between the canopy solve
  and the soil-heat ODE, iterating the two to a jointly self-consistent state
  within the hour by default (`IterationToleranceConvergence(; tolerance=0.05u"K",
  max_iterations_per_day=80)`); `FixedIterationConvergence(1)` for a single
  pass using the previous hour's soil temperature instead.
- `canopy_soil_relaxation`: under-relaxation on each `canopy_soil_convergence`
  pass's `ground_temperature` update (`x_new = relaxation*x_solved +
  (1-relaxation)*x_prev`), `1.0` (no relaxation) by default.
"""
@kwdef struct MicroConfig{CV,RFS,SMM,MSP,MPD,CSC,CSR}
    convergence::CV = FixedIterationConvergence(3)
    rainfall_schedule::RFS = DailyRainfall()
    soil_moisture_strategy::SMM = PrescribedSoilMoisture()
    max_surface_pool::MSP = 1.0e4u"kg/m^2"
    max_pond_depth::MPD = 1.0e4u"kg/m^2"
    canopy_soil_convergence::CSC = IterationToleranceConvergence(; tolerance=0.05u"K", max_iterations_per_day=80)
    canopy_soil_relaxation::CSR = 1.0
end

"""
    MicroModel(; hours, depths, heights,
                 soil_properties_model, soil_hydraulic_model,
                 radiation=RadiationModel(),
                 snow_model=NoSnow(),
                 vapour_pressure_equation=GoffGratch(),
                 boundary_layer_model=MoninObukhov(),
                 evaporation_model=BulkTransferEvaporation(),
                 soil_energy_model=SoilHeatTransport1D(),
                 config=MicroConfig())

Constant-across-runs scientific description of the simulation:

- solver geometry (`hours`, `depths`, `heights`) # TODO generalise to sub-hourly
- physical-process models:
    - `soil_properties_model::AbstractSoilProperties` (e.g. `CampbelldeVriesSoilProperties(...)`)
    - `soil_hydraulic_model::AbstractSoilHydraulicsModel` (e.g. `CampbellSoilHydraulics(...)`)
    - `radiation::RadiationModel` — bundle of `solar_radiation_model`,
      `longwave_model`, and `shortwave_model`Anchored 
    - `snow_model::AbstractSnowModel` — `NoSnow()` (default) or `SnowModel(...)`
    - `vapour_pressure_equation` — cross-cutting (`GoffGratch()` / `Teten()` / `Huang()`)
    - `boundary_layer_model` — cross-cutting (`MoninObukhov()`)
    - `evaporation_model` — surface latent flux
    - `soil_energy_model::SoilHeatTransportModel` — soil column energy ODE
      (carries the phase-transition `freezing_model` and ODE solver settings)
    - `canopy_model::AbstractCanopyModel` — `NoCanopy()` (default, today's
      scalar-`shade` vegetation handling, unchanged) or `MultilayerCanopy(...)`
      (layer-resolved two-stream radiative transfer through the canopy)
- iteration/data-delivery strategy in `config::MicroConfig`

Combine with a `MicroInputs` via `MicroProblem(model, inputs; days)` to run.
"""
@kwdef struct MicroModel{H,Dep,Ht,SPM,SHM,RAD,SNM,VPE,BLM,EVM,SEM,CAN,C}
    hours::H = DEFAULT_HOURS # hour of day for solar_radiation
    depths::Dep = DEFAULT_DEPTHS # soil nodes - keep spacing close near the surface
    heights::Ht = [0.01, 2]u"m" # air nodes for temperature, wind speed and humidity profile, last height is reference height for weather data
    soil_properties_model::SPM
    soil_hydraulic_model::SHM
    radiation::RAD = RadiationModel()
    snow_model::SNM = NoSnow()
    vapour_pressure_equation::VPE = GoffGratch()
    boundary_layer_model::BLM = MoninObukhov()
    evaporation_model::EVM = BulkTransferEvaporation()
    soil_energy_model::SEM = SoilHeatTransport1D()
    canopy_model::CAN = NoCanopy()
    config::C = MicroConfig()
end

"""
    MicroInputs(; site, soil_profile, environment_minmax, environment_daily,
                  environment_hourly=nothing,
                  initial_soil_temperature=nothing, initial_soil_moisture,
                  initial_snow_depth=0.0u"cm", initial_snow_temperature=u"K"(0.0u"°C"),
                  initial_snow_density=nothing)

Per-run input data: site, soil column profile, environment forcings, and
initial conditions.

- `site::Site` — properties of the place
- `soil_profile::SoilProfile` — per-depth `bulk_density`, `mineral_density`,
  `mineral_conductivity`, and `mineral_heat_capacity` read by both the soil
  properties and soil hydraulic models. Per-location structural soil-column
  data; lives here (not on `MicroModel`) because it varies between sites
  without changing the physics.
- `environment_minmax`, `environment_daily`, `environment_hourly` — input forcings
- `initial_*` — initial conditions (`initial_soil_temperature=nothing` falls back to
  the day-mean reference air temperature)
- `initial_snow_density=nothing` → use `model.snow_model.snow_density`

Pass a fresh `MicroInputs` to `reinit!(cache, inputs)` to re-solve the same
model on different data without rebuilding the cache.
"""
@kwdef struct MicroInputs{S,SP<:SoilProfile,EMM,EH,ED,IST,ISM,ISD,ISNT,ISND}
    site::S
    soil_profile::SP
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

function example_microclimate_problem(;
    days = DEFAULT_DAYS,
    hours = DEFAULT_HOURS,
    depths = DEFAULT_DEPTHS,
    heights = [0.01, 2]u"m",
    site = example_site(),
    soil_profile = example_soil_profile(depths),
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model = example_soil_hydraulic_model(),
    snow_model = NoSnow(),
    canopy_model = NoCanopy(),
    environment_minmax = example_monthly_weather(),
    environment_daily = example_daily_environment(days),
    environment_hourly = example_hourly_environment(days, hours; elevation=site.elevation),
    config = MicroConfig(),
    initial_soil_temperature = fill(u"K"(7.741667u"°C"), length(depths)),
    initial_soil_moisture = fill(0.42 * 0.25, length(depths)),
)
    model = MicroModel(;
        hours, depths, heights, soil_properties_model, soil_hydraulic_model, snow_model, canopy_model, config
    )
    inputs = MicroInputs(;
        site, soil_profile, environment_minmax, environment_daily, environment_hourly,
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
    soil_water_profile::P          # soil moisture profile scratch used by the moisture solver
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
