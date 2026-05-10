# Microclimate.jl structural refactor plan

## Goals

1. **Modular & swappable formulations.** Every physics choice (soil thermal, hydraulics, snow, radiation, boundary layer, …) is an abstract type with concrete implementations the user picks. Adding a new formulation = adding one new file.
2. **All allocation up-front, in one object.** `MicroCache` owns every buffer the simulation will ever need. The hot path allocates nothing.
3. **One parameters object.** All physical parameters live in a single `MicroParameters` struct (which composes per-physics-concept structs).
4. **One config object.** All formulation choices and solver tuning live in a single `MicroConfig` struct.
5. **One site object.** All site/place properties (latitude, slope, roughness, …) live in a single `Site` struct, replacing the duplicated `MicroTerrain` + `SolarTerrain`.
6. **Clean naming.** Drop the `Micro` prefix from anything that isn't a top-level container. Drop Fortran-isms. Drop `2`-suffixes. No abstract-typed or untyped struct fields.

---

## Conventions

### Struct field types — always parameterize

Every Julia struct in this package must give each user-data field a type parameter. Never leave fields as `Any` (implicit) or as an abstract type — that boxes the value and tanks performance in hot paths. The whole `ALLOCATION_AUDIT.md` effort exists to remove exactly this dynamic dispatch.

```julia
# WRONG — fields box, callsite is type-unstable
struct Site
    latitude
    elevation
end

# RIGHT — parametric, concrete at each call site
@kwdef struct Site{LAT,EL}
    latitude::LAT
    elevation::EL
end
```

Concrete fixed types (`Int`, `Bool`, `Float64`) can be written directly with no parameter. The rule is: never `Any`, never abstract.

### File-per-formulation

Each `<concept>/` folder contains one `abstract.jl` declaring the abstract type and the interface methods, plus one file per concrete formulation. Each formulation file holds: the parameter struct, its constructor, any `allocate_*` for buffers it owns, and every dispatched method on the abstract type. Adding a new formulation = adding one file + one `include` + one `export`.

### Naming

- Drop the `Micro` prefix on non-container types (`MicroTerrain` → gone, `MicroForcing` → `Forcing`, `MicroProfile` → `AtmosphericProfile`).
- Keep `Micro` on top-level containers: `MicroProblem`, `MicroCache`, `MicroResult`, `MicroConfig`, `MicroParameters`.
- No `2` suffixes (`soil_bulk_density2` → `bulk_density`, deduped).
- No Fortran-isms in user-facing names (`ndmax` → `iterations_per_day`, `qfreze` → `Q_freeze`, `moist_step` → `moisture_timestep`).
- Distinguish *model* (parameters) from *strategy* (solver behaviour). `AbstractSoilMoistureMode` → `AbstractSoilMoistureStrategy`.

---

## Target file layout

```
src/
  Microclimate.jl
  constants.jl
  problem.jl                       # MicroProblem
  config.jl                        # MicroConfig
  parameters.jl                    # MicroParameters (composes per-concept structs)
  cache.jl                         # MicroCache, init, reinit!
  solve.jl                         # solve!, day loop
  outputs.jl                       # MicroResult, AtmosphericProfile
  forcing.jl                       # Forcing + update_forcing_day!
  energy_balance.jl                # the ODE RHS
  phase_transition.jl
  boundary_layer.jl                # atmospheric profile, Obukhov

  inputs/                          # forcings (one per shape)
    monthly_minmax.jl
    daily_timeseries.jl
    hourly_timeseries.jl

  site/
    abstract.jl
    site.jl                        # the unified Site struct

  soil_thermal/
    abstract.jl                    # AbstractSoilThermalModel + interface
    campbell_devries.jl            # struct + soil_properties + soil_properties!

  soil_hydraulics/
    abstract.jl
    campbell.jl                    # struct + soil_water_balance + buffers

  soil_moisture/                   # the strategy (was "mode")
    abstract.jl
    prescribed.jl
    dynamic.jl

  snow/
    abstract.jl                    # AbstractSnowModel
    no_snow.jl                     # NoSnow + every no-op method
    snow_model.jl                  # SnowModel + algorithms (split further if needed)

  time_mode/
    abstract.jl
    consecutive.jl
    non_consecutive.jl

  convergence/
    abstract.jl
    fixed.jl
    tolerance.jl

  boundary_layer_model/
    abstract.jl
    monin_obukhov.jl               # MoninObukhov + κ, γ constants

  atmospheric_radiation/
    abstract.jl
    swinbank.jl
    campbell_norman.jl

  cloud_adjust/
    abstract.jl
    angstrom.jl                    # currently hard-coded constants

  diffuse/
    abstract.jl
    erbs.jl

  rainfall/
    abstract.jl
    daily_rainfall.jl
    hourly_rainfall.jl
```

A folder is overkill when there's truly one implementation and the concept is small (`phase_transition.jl`, `boundary_layer.jl`, `forcing.jl` stay flat). Promote to a folder the moment a second implementation appears.

---

## The five top-level user-facing structs

```julia
@kwdef struct MicroProblem{G,S,P,I,IC,C}
    grid::G                        # MicroGrid: days, hours, depths, heights
    site::S                        # Site (issue 6 below)
    parameters::P                  # MicroParameters
    inputs::I                      # MicroInputs: env_minmax, env_daily, env_hourly
    initial::IC                    # MicroInitial: T_soil, θ_soil, snow_*
    config::C                      # MicroConfig
end
```

This shrinks `MicroProblem` from its current 21 type parameters to 6, each grouping a clear category.

---

## Issue 1 — `MicroTerrain` and `SolarTerrain` collapse into `Site`

Current `MicroTerrain` jumbles three different things: site geometry (`elevation`, `viewfactor`), surface property (`roughness_height`), and boundary-layer model coefficients (`karman_constant`, `dyer_constant`). The last two are properties of the *Monin–Obukhov formulation*, not of the site. `elevation` is also duplicated on `SolarTerrain`. `albedo` lives on `SolarTerrain` but is overridden by snow.

Replace both with one `Site`:

```julia
@kwdef struct Site{LAT,LON,EL,SL,AS,HA,SVF,AB,RH}
    latitude::LAT
    longitude::LON
    elevation::EL
    slope::SL
    aspect::AS
    horizon_angles::HA
    sky_view_fraction::SVF         # was viewfactor
    albedo::AB
    roughness_height::RH           # z₀ — a property of the patch of ground
end
```

Lift the boundary-layer constants out into a model that lives in `MicroConfig`:

```julia
abstract type AbstractBoundaryLayerModel end

@kwdef struct MoninObukhov{KC,DC} <: AbstractBoundaryLayerModel
    karman_constant::KC = 0.4
    dyer_constant::DC = 16.0
end
```

Snow's albedo override stops being weird: the override is temporary state in `cache.state.snow`, the baseline lives on `Site`.

---

## Issue 2 — `MicroConfig` extraction

All the current `MicroProblem` fields that are formulation choices or solver tuning move into one struct:

```julia
@kwdef struct MicroConfig{VPE,SOS,SOK,TM,CV,DFM,BLM,ARM,CAM,RFS,SNM,MST}
    vapour_pressure_equation::VPE
    soil_ode_solver::SOS
    soil_ode_kwargs::SOK
    time_mode::TM
    convergence::CV
    diffuse_fraction_model::DFM
    boundary_layer_model::BLM
    atmospheric_radiation_model::ARM        # currently hard-coded inside precompute_longwave_sky
    cloud_adjust_model::CAM                 # currently hard-coded a, b, gamma constants
    rainfall_schedule::RFS                  # replaces hourly_rainfall::Bool with a dispatched type
    snow_model::SNM                         # may move to parameters depending on framing
    maximum_surface_temperature::MST
end
```

Solver tuning currently inside `CampbellSoilHydraulics` (`moist_error`, `moist_count`, `moist_step`, `maxpool`) also moves here — it's not Campbell-specific.

---

## Issue 3 — `MicroParameters` dedup

Currently `bulk_density` and `mineral_density` are stored in *both* `CampbelldeVriesSoilThermal` *and* `CampbellSoilHydraulics` (as `soil_bulk_density2` / `soil_mineral_density2`). The user has to set the same numbers twice.

```julia
@kwdef struct MicroParameters{ST,SH,SN}
    soil_thermal::ST                # AbstractSoilThermalModel
    soil_hydraulics::SH             # AbstractSoilHydraulicsModel (renamed from "moisture model")
    snow::SN                        # AbstractSnowModel (alternatively in MicroConfig)
end
```

`bulk_density` and `mineral_density` live in *one* concept's struct (probably `soil_hydraulics`, which has the per-depth profile) and the thermal model reads them via composition or via fields shared on `MicroParameters`. Drop the `2` suffixes. Drop the `mode` field from `CampbellSoilHydraulics` (it belongs in config, not in parameters).

---

## Issue 4 — modularity gaps to close

- **Snow has no abstract type.** Currently `Union{NoSnow, SnowModel}`. Introduce `AbstractSnowModel` with `NoSnow{} <: AbstractSnowModel` and `SnowModel{N,…} <: AbstractSnowModel`.
- **Atmospheric radiation model is hard-coded.** `precompute_longwave_sky` defaults to `CampbellNormanAtmosphericRadiation()`. `Swinbank` exists but is unreachable. Thread through `MicroConfig.atmospheric_radiation_model`.
- **Cloud-adjustment is hard-coded.** `cloud_adjust_radiation`'s `a, b, gamma` Ångström constants are buried as kwargs and never exposed. Lift to `AbstractCloudAdjustModel` with `Angstrom{A,B,G}` as the first concrete implementation.
- **Rainfall is a Bool.** `hourly_rainfall::Bool` should be a dispatched `AbstractRainfallSchedule` with `DailyRainfall` / `HourlyRainfall` so it follows the same trait pattern as everything else.

---

## Issue 5 — `MicroCache` reorganisation

Today `MicroCache` mixes buffers, working state, and one-of-each helpers as flat fields:

```
problem, output, solar_radiation_out, solar_buffers, buffers, ode_integrator,
soil_moisture, nodes, nodes_day, ∑phase, profile_buffers, snow_state,
snow_scratch, forcing, num_soil_nodes_val
```

Reorganise into clear halves:

```julia
mutable struct MicroCache{P,O,S,B,F,I}
    problem::P
    output::O
    state::S                       # MicroState: snow_state, soil_moisture, nodes, ∑phase, surface_pool, …
    buffers::B                     # MicroBuffers: solar, profile, soil_props, soil_energy, soil_water, snow_scratch, phase
    forcing::F                     # Forcing
    ode_integrator::I
end
```

`MicroBuffers(problem)` is a thin assembler that calls each concept's `allocate_*`. "All allocation up-front, in one object" then becomes literally true and provable — it's exactly that one constructor call.

---

## Issue 6 — output split

`MicroResult` currently mixes:
- echoed inputs (`pressure`, `reference_temperature`, `cloud_cover`)
- solar radiation NamedTuple (a SolarRadiation.jl scratch+output)
- modelled outputs (`soil_temperature`, profile, snow_*)

It also subtypes `AbstractEnvironment` alongside `DailyTimeseries` / `HourlyTimeseries`, conflating "what you give the model" with "what you get back".

Split: `MicroResult` holds modelled outputs only. Input echo accessed via `cache.problem.inputs`. Give input environments their own abstract supertype distinct from outputs.

---

## Issue 7 — naming sweep

| Current | New |
|---|---|
| `soil_bulk_density2` | `bulk_density` (deduped) |
| `soil_mineral_density2` | `mineral_density` (deduped) |
| `moist_error` | `moisture_tolerance` |
| `moist_count` | `moisture_max_iterations` |
| `moist_step` | `moisture_timestep` |
| `maxpool` | `max_surface_pool` |
| `ndmax` (on `NonConsecutiveDayMode`) | `iterations_per_day` |
| `qfreze` | `Q_freeze` |
| `hourly_rainfall::Bool` | `rainfall_schedule::AbstractRainfallSchedule` |
| `viewfactor` | `sky_view_fraction` |
| `MicroTerrain` | absorbed into `Site` |
| `MicroForcing` | `Forcing` |
| `MicroProfile` | `AtmosphericProfile` |
| `AbstractSoilMoistureMode` | `AbstractSoilMoistureStrategy` |
| `interpolate_solar`, `interpolate_zenith`, … (fields of `Forcing`) | `solar`, `zenith`, … |

---

## Issue 8 — broken `example_*` constructors

`example_monthly_weather()` calls `MonthlyMinMaxWeather(...)` but only `MonthlyMinMaxEnvironment` exists. `example_daily_environmental(...)` calls `EnvironmentTimeseries(; albedo, ...)` but only `DailyTimeseries` exists, and `albedo` isn't an argument. These have rotted since a rename and aren't covered by tests.

Fix the names, then add a smoke test that calls each `example_*` constructor so they can't rot again.

---

## Issue 9 — `allocate_*` API uniformity

Today: some `allocate_*` return NamedTuples, one returns a struct (`Forcing`), one returns a 2-tuple (`allocate_solar`), and `allocate_snow_scratch(::NoSnow, …)` returns `nothing`.

Standardise: every `allocate_*(model_or_problem)` returns a struct (parametric per field per the convention above), dispatched on the model so `NoSnow` returns an empty zero-cost struct rather than `nothing`. This removes one of the patterns the allocation audit is still working through.

---

## Migration order — smallest blast radius first

Each step is a self-contained refactor that leaves the existing solver intact. Run the test suite (`micro_testrun_monthly.jl`, `micro_testrun_daily.jl`, snow tests) between each.

1. **Fix broken `example_*` constructors.** Add a smoke test. (Zero risk; unblocks downstream renames.)
2. **Split `landscape.jl` into the conceptual file layout.** Pure mechanical move, no behaviour change. `landscape.jl` disappears.
3. **Introduce `AbstractSnowModel`.** Migrate the `Union{NoSnow, SnowModel}` dispatch sites. Aligns snow with the rest of the codebase pattern.
4. **Collapse `MicroTerrain` + `SolarTerrain` → `Site`.** Lift `karman_constant` / `dyer_constant` into `MoninObukhov` (in config). Dedupe `elevation`, `albedo`.
5. **Dedupe `bulk_density` / `mineral_density`** between thermal and hydraulics structs. Drop `2` suffixes.
6. **Extract `MicroConfig`.** Move all formulation/solver fields off `MicroProblem`. Move `mode`, `moist_*`, `maxpool` out of `CampbellSoilHydraulics`.
7. **Extract `MicroParameters`.** Group the per-concept parameter structs.
8. **Reorganise `MicroCache` into `state` / `buffers` halves.** Standardise `allocate_*` to return structs.
9. **Lift atmospheric radiation, cloud-adjust constants, and rainfall scheduling** into the trait pattern so they're swappable from the top.
10. **Split `MicroResult`.** Stop subtyping `AbstractEnvironment`; give inputs their own abstract supertype.
11. **Naming sweep.** Apply the renames in the table above. Drop the `Micro` prefix from non-container types.

Steps 1–5 are mostly mechanical and can land quickly. Steps 6–10 reshape the public API and want one combined release / migration note.
