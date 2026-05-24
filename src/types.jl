"""
    SoilProfile(; bulk_density, mineral_density)

Per-depth soil column profile read by both the soil properties model
(thermal) and the soil hydraulic model. Sized to match `MicroModel.depths`.

- `bulk_density` — dry soil bulk density (kg/m³) at each depth node
- `mineral_density` — soil mineral density (kg/m³) at each depth node
"""
@kwdef struct SoilProfile{BD,MD}
    bulk_density::BD
    mineral_density::MD
end

function example_soil_profile(depths=DEFAULT_DEPTHS;
    bulk_density = 2.56u"Mg/m^3",
    mineral_density = 2.560u"Mg/m^3",
)
    SoilProfile(;
        bulk_density = bulk_density isa AbstractVector ? bulk_density : fill(bulk_density, length(depths)),
        mineral_density = mineral_density isa AbstractVector ? mineral_density : fill(mineral_density, length(depths)),
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

- `time_mode`: `NonConsecutiveDayMode()` or `ConsecutiveDayMode(; spinup_first_day=false)`
- `convergence`: `FixedSoilTemperatureIterations(3)` or `SoilTemperatureConvergenceTolerance(; tolerance, max_iterations_per_day)`
- `rainfall_schedule`: `DailyRainfall()` (default) or `HourlyRainfall()`
- `soil_moisture_strategy`: `PrescribedSoilMoisture()` or `DynamicSoilMoisture(; ...)`
- `max_surface_pool`: numerical clamp on the surface-pool state variable
  (not a physical limit — keeps the pool integration from running away)
"""
@kwdef struct MicroConfig{TM,CV,RFS,SMM,MSP}
    time_mode::TM = NonConsecutiveDayMode()
    convergence::CV = FixedSoilTemperatureIterations(3)
    rainfall_schedule::RFS = DailyRainfall()
    soil_moisture_strategy::SMM = PrescribedSoilMoisture()
    max_surface_pool::MSP = 1.0e4u"kg/m^2"
end

"""
    MicroModel(; days, hours, depths, heights,
                 soil_profile, soil_properties_model, soil_hydraulic_model,
                 radiation=RadiationModel(),
                 snow_model=NoSnow(),
                 vapour_pressure_equation=GoffGratch(),
                 boundary_layer_model=MoninObukhov(),
                 evaporation_model=BulkTransferEvaporation(),
                 soil_energy_model=SoilHeatTransport1D(),
                 config=MicroConfig())

Constant-across-runs scientific description of the simulation:

- solver geometry (`days`, `hours`, `depths`, `heights`)
- physical-process models:
    - `soil_profile::SoilProfile` — per-depth `bulk_density` and `mineral_density`
      profiles read by both the soil properties and soil hydraulic models
    - `soil_properties_model::AbstractSoilProperties` (e.g. `CampbelldeVriesSoilProperties(...)`)
    - `soil_hydraulic_model::AbstractSoilHydraulicsModel` (e.g. `CampbellSoilHydraulics(...)`)
    - `radiation::RadiationModel` — bundle of `solar_radiation_model`,
      `longwave_model`, and `shortwave_model`
    - `snow_model::AbstractSnowModel` — `NoSnow()` (default) or `SnowModel(...)`
    - `vapour_pressure_equation` — cross-cutting (`GoffGratch()` / `Teten()` / `Huang()`)
    - `boundary_layer_model` — cross-cutting (`MoninObukhov()`)
    - `evaporation_model` — surface latent flux
    - `soil_energy_model::SoilHeatTransportModel` — soil column energy ODE
      (carries the phase-transition `freezing_model` and ODE solver settings)
- iteration/data-delivery strategy in `config::MicroConfig`

Combine with a `MicroInputs` via `MicroProblem(model, inputs)` to run.
"""
@kwdef struct MicroModel{D,H,Dep,Ht,SPR,SPM,SHM,RAD,SNM,VPE,BLM,EVM,SEM,C}
    days::D = DEFAULT_DAYS # days of year to simulate - TODO leap years - why not use real dates?
    hours::H = DEFAULT_HOURS # hour of day for solar_radiation
    depths::Dep = DEFAULT_DEPTHS # soil nodes - keep spacing close near the surface
    heights::Ht = [0.01, 2]u"m" # air nodes for temperature, wind speed and humidity profile, last height is reference height for weather data
    soil_profile::SPR
    soil_properties_model::SPM
    soil_hydraulic_model::SHM
    radiation::RAD = RadiationModel()
    snow_model::SNM = NoSnow()
    vapour_pressure_equation::VPE = GoffGratch()
    boundary_layer_model::BLM = MoninObukhov()
    evaporation_model::EVM = BulkTransferEvaporation()
    soil_energy_model::SEM = SoilHeatTransport1D()
    config::C = MicroConfig()
end

"""
    MicroInputs(; site, environment_minmax, environment_daily, environment_hourly=nothing,
                  initial_soil_temperature=nothing, initial_soil_moisture,
                  initial_snow_depth=0.0u"cm", initial_snow_temperature=u"K"(0.0u"°C"),
                  initial_snow_density=nothing)

Per-run input data: site, environment forcings, and initial conditions.

- `site::Site` — properties of the place
- `environment_minmax`, `environment_daily`, `environment_hourly` — input forcings
- `initial_*` — initial conditions (`initial_soil_temperature=nothing` falls back to
  the day-mean reference air temperature)
- `initial_snow_density=nothing` → use `model.snow_model.snow_density`

Pass a fresh `MicroInputs` to `reinit!(cache, inputs)` to re-solve the same
model on different data without rebuilding the cache.
"""
@kwdef struct MicroInputs{S,EMM,EH,ED,IST,ISM,ISD,ISNT,ISND}
    site::S
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
    MicroProblem(model::MicroModel, inputs::MicroInputs)

Pairing of a model description with a set of input data. Pass to
`init(problem)` to allocate a `MicroCache` or `solve(problem)` to run.
"""
struct MicroProblem{M<:MicroModel,I<:MicroInputs}
    model::M
    inputs::I
end

function example_microclimate_problem(;
    days = DEFAULT_DAYS,
    hours = DEFAULT_HOURS,
    depths = DEFAULT_DEPTHS,
    heights = [0.01, 2]u"m",
    site = example_site(),
    soil_profile = example_soil_profile(depths),
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model = example_soil_hydraulic_model(depths),
    snow_model = NoSnow(),
    environment_minmax = example_monthly_weather(),
    environment_daily = example_daily_environment(days),
    environment_hourly = example_hourly_environment(days, hours; elevation=site.elevation),
    config = MicroConfig(),
    initial_soil_temperature = fill(u"K"(7.741667u"°C"), length(depths)),
    initial_soil_moisture = fill(0.42 * 0.25, length(depths)),
)
    model = MicroModel(; days, hours, depths, heights,
        soil_profile, soil_properties_model, soil_hydraulic_model, snow_model, config)
    inputs = MicroInputs(;
        site, environment_minmax, environment_daily, environment_hourly,
        initial_soil_temperature, initial_soil_moisture,
    )
    MicroProblem(model, inputs)
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
struct MicroBuffers{SO,SOB,P,PB,SEB,SP,PT,SWB,SS,IB}
    solar_out::SO                  # SolarRadiation output (NamedTuple of arrays)
    solar::SOB                     # SolarRadiation internal buffers (NamedTuple)
    soil_water_profile::P          # atmospheric profile scratch used by the moisture solver
    air_profile::PB                # atmospheric profile scratch used by solve_air!
    soil_energy_balance::SEB
    soil_properties::SP
    phase_transition::PT
    soil_water_balance::SWB
    snow::SS                       # snow buffers (NamedTuple for SnowModel; nothing for NoSnow)
    interpolation::IB              # 24-element scratch buffers for hourly_from_min_max!
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
