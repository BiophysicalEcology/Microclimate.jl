"""
    MicroConfig(; ...)

Single home for the simulation's formulation choices and solver tuning.
Lives on `MicroProblem.config`. Kept separate from `MicroProblem`'s grid,
site, parameters, and inputs so users can swap any one category
independently.

# Formulation choices (all type-dispatched)
- `vapour_pressure_equation`: `GoffGratch()`, `Teten()`, or `Huang()`
- `boundary_layer_model`: `MoninObukhov()` (κ, γ constants)
- `time_mode`: `NonConsecutiveDayMode()` or `ConsecutiveDayMode(; spinup_first_day=false)`
- `convergence`: `FixedSoilTemperatureIterations(3)` or `SoilTemperatureConvergenceTolerance(; tolerance, max_iterations_per_day)`
- `diffuse_fraction_model`: `ErbsDiffuseFraction()`
- `atmospheric_radiation_model`: `CampbellNormanAtmosphericRadiation()` or `SwinbankAtmosphericRadiation()`
- `cloud_adjust_model`: `Angstrom(; a, b, gamma)` — Ångström–Prescott scaling
- `rainfall_schedule`: `DailyRainfall()` (default) or `HourlyRainfall()`
- `soil_moisture_strategy`: `PrescribedSoilMoisture()` or `DynamicSoilMoisture()`

# Solver options
- `soil_ode_solver`: any SciML algorithm (`Tsit5()`, `RK4()`, …)
- `soil_ode_kwargs`: NamedTuple of kwargs forwarded to the integrator

# Soil-moisture solver tuning (used only when `soil_moisture_strategy = DynamicSoilMoisture()`)
- `moisture_tolerance`: maximum mass-balance residual
- `moisture_max_iterations`: maximum iterations of the moisture solver
- `moisture_timestep`: timestep of the moisture solver (≤ 1 hour)
- `max_surface_pool`: maximum surface-pool depth

# Other
- `maximum_surface_temperature`: surface-temperature safety clamp (Fortran microinput[74])
"""
@kwdef struct MicroConfig{VPE,BLM,TM,CV,DFM,ARM,CAM,RFS,SMM,SOS,SOK,MSF,MT,MTS,MSP}
    vapour_pressure_equation::VPE = GoffGratch()
    boundary_layer_model::BLM = MoninObukhov()
    time_mode::TM = NonConsecutiveDayMode()
    convergence::CV = FixedSoilTemperatureIterations(3)
    diffuse_fraction_model::DFM = ErbsDiffuseFraction()
    atmospheric_radiation_model::ARM = CampbellNormanAtmosphericRadiation()
    cloud_adjust_model::CAM = Angstrom()
    rainfall_schedule::RFS = DailyRainfall()
    soil_moisture_strategy::SMM = PrescribedSoilMoisture()
    soil_ode_solver::SOS = Tsit5()
    soil_ode_kwargs::SOK = (; reltol=1e-6u"K", abstol=1e-3u"K")
    maximum_surface_temperature::MSF = 85.0u"°C"
    moisture_tolerance::MT = 1e-6u"kg/m^2/s"
    moisture_max_iterations::Int = 500
    moisture_timestep::MTS = 360.0u"s"
    max_surface_pool::MSP = 1.0e4u"kg/m^2"
end
