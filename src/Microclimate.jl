module Microclimate

using CommonSolve: CommonSolve
using ConstructionBase
using Interpolations, Statistics, Dates
using SciMLBase, OrdinaryDiffEqTsit5
using OrdinaryDiffEqTsit5: Tsit5
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


export MicroProblem, MicroCache

export GoffGratch, Teten, Huang
export Tsit5

# Snow model
export NoSnow, SnowModel, SnowState

# Soil thermal model
export CampbelldeVriesSoilThermal

# Soil hydraulics and moisture mode
export CampbellSoilHydraulics
export PrescribedSoilMoisture, DynamicSoilMoisture

# Convergence strategies
export AbstractSoilTemperatureConvergence, FixedSoilTemperatureIterations, SoilTemperatureConvergenceTolerance

# Time modes
export AbstractTimeMode, NonConsecutiveDayMode, ConsecutiveDayMode

# Diffuse fraction models
export AbstractDiffuseFractionModel, ErbsDiffuseFraction

export MonthlyMinMaxEnvironment, DailyMinMaxEnvironment, DailyTimeseries, HourlyTimeseries, MicroTerrain

export daily_cycle_sine_exponential, daily_cycle_linear, hourly_from_min_max

export cloud_adjust_radiation, longwave_radiation, precompute_longwave_sky

export atmospheric_surface_profile, calc_convection

export soil_properties, soil_properties!, allocate_soil_properties

export soil_energy_balance, evaporation, soil_water_balance!, phase_transition

export example_micro_terrain, example_monthly_weather, example_daily_environmental, example_soil_hydraulics, example_soil_thermal_parameters, example_microclimate_problem

import CommonSolve: solve, solve!, init
export solve, solve!, init, reinit!


include("constants.jl")
include("landscape.jl")
include("interpolation.jl")
include("soil_properties.jl")
include("radiation.jl")
include("boundary_layer.jl")
include("soil_balance.jl")
include("snow.jl")
include("simulation.jl")

end
