module Microclimate

using CommonSolve: CommonSolve
using ConstructionBase
using Interpolations, Statistics, Dates
using SciMLBase, OrdinaryDiffEqTsit5, OrdinaryDiffEqAdamsBashforthMoulton
using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqAdamsBashforthMoulton: AB3, AB4, AB5, ABM32, ABM43, ABM54
using Unitful, UnitfulMoles
using ModelParameters, DelimitedFiles
using SpecialFunctions, StaticArrays

using FluidProperties: atmospheric_pressure, wet_air_properties, dry_air_properties, vapour_pressure
using FluidProperties: enthalpy_of_vaporisation, molar_enthalpy_of_vaporisation, water_properties
using FluidProperties: g_n, σ, atm, R
using FluidProperties: GoffGratch, Teten, Huang
using Unitful: °, rad, °C
using Interpolations: AbstractInterpolation
using SolarRadiation


export MicroProblem, MicroCache, MicroConfig, MicroParameters

export GoffGratch, Teten, Huang
export Tsit5
export AB3, AB4, AB5, ABM32, ABM43, ABM54

# Snow model
export AbstractSnowModel, NoSnow, SnowModel, SnowState

# Soil thermal model
export CampbelldeVriesSoilThermal

# Soil hydraulics
export AbstractSoilHydraulicsModel, CampbellSoilHydraulics

# Soil moisture strategy (lives on MicroConfig.soil_moisture_strategy)
export AbstractSoilMoistureStrategy, PrescribedSoilMoisture, DynamicSoilMoisture

# Convergence strategies
export AbstractSoilTemperatureConvergence, FixedSoilTemperatureIterations, SoilTemperatureConvergenceTolerance

# Time modes
export AbstractTimeMode, NonConsecutiveDayMode, ConsecutiveDayMode

# Diffuse fraction models
export AbstractDiffuseFractionModel, ErbsDiffuseFraction

# Atmospheric radiation models
export AbstractAtmosphericRadiationModel, SwinbankAtmosphericRadiation, CampbellNormanAtmosphericRadiation, atmospheric_radiation

# Cloud adjustment models
export AbstractCloudAdjustModel, Angstrom, sunshine_fraction

# Rainfall schedule
export AbstractRainfallSchedule, DailyRainfall, HourlyRainfall, is_hourly

export AbstractEnvironment, MonthlyMinMaxEnvironment, DailyMinMaxEnvironment, DailyTimeseries, HourlyTimeseries
export Site, AbstractSite
export AbstractBoundaryLayerModel, MoninObukhov
export Forcing, AtmosphericProfile


export cloud_adjust_radiation, longwave_radiation, precompute_longwave_sky

export atmospheric_surface_profile, calc_convection

export soil_properties, soil_properties!, allocate_soil_properties

export soil_energy_balance, evaporation, soil_water_balance!, phase_transition

export example_site, example_monthly_weather,
    example_daily_environment, example_hourly_environment,
    example_soil_hydraulics, example_soil_thermal_parameters, example_microclimate_problem

import CommonSolve: solve, solve!, init
export solve, solve!, init, reinit!


include("constants.jl")

# Soil thermal parameters
include("soil_thermal/abstract.jl")
include("soil_thermal/campbell_devries.jl")

# Soil moisture strategy (PrescribedSoilMoisture is the default for CampbellSoilHydraulics.mode,
# so it must be loaded before soil_hydraulics).
include("soil_moisture/abstract.jl")
include("soil_moisture/prescribed.jl")
include("soil_moisture/dynamic.jl")

# Soil hydraulics parameters (Mode field defaults to PrescribedSoilMoisture)
include("soil_hydraulics/abstract.jl")
include("soil_hydraulics/campbell.jl")

# Site (place properties)
include("site/abstract.jl")
include("site/site.jl")

# Boundary-layer model
include("boundary_layer_model/abstract.jl")
include("boundary_layer_model/monin_obukhov.jl")

# Environment inputs
include("inputs/abstract_environment.jl")
include("inputs/monthly_minmax.jl")
include("inputs/daily_minmax.jl")
include("inputs/daily_timeseries.jl")
include("inputs/hourly_timeseries.jl")

# Convergence (loaded before time_mode — iterations_for_day calls max_iterations)
include("convergence/abstract.jl")
include("convergence/fixed.jl")
include("convergence/tolerance.jl")

# Time mode
include("time_mode/abstract.jl")
include("time_mode/non_consecutive.jl")
include("time_mode/consecutive.jl")

# Diffuse fraction
include("diffuse/abstract.jl")
include("diffuse/erbs.jl")

# Atmospheric radiation
include("atmospheric_radiation/abstract.jl")
include("atmospheric_radiation/swinbank.jl")
include("atmospheric_radiation/campbell_norman.jl")

# Cloud-cover adjustment of clear-sky radiation
include("cloud_adjust/abstract.jl")
include("cloud_adjust/angstrom.jl")

# Rainfall schedule
include("rainfall/abstract.jl")
include("rainfall/daily_rainfall.jl")
include("rainfall/hourly_rainfall.jl")

# Outputs and forcing (depend on AbstractEnvironment / AbstractInterpolation)
include("outputs.jl")
include("forcing.jl")

# Algorithms
include("interpolation.jl")
include("soil_properties.jl")
include("radiation.jl")
include("boundary_layer.jl")
include("soil_balance.jl")
include("snow.jl")

# Top-level config, parameters, and problem
include("config.jl")
include("parameters.jl")
include("simulation.jl")

end
