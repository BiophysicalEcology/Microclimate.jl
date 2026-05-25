module Microclimate

using CommonSolve: CommonSolve
using ConstructionBase
using Interpolations, Statistics, Dates
using SciMLBase, OrdinaryDiffEqTsit5, OrdinaryDiffEqAdamsBashforthMoulton
using OrdinaryDiffEqCore: ismultistep
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


export MicroProblem, MicroModel, MicroInputs, MicroCache, MicroConfig
export SoilProfile, RadiationModel
export example_soil_profile

export GoffGratch, Teten, Huang
export Tsit5
export AB3, AB4, AB5, ABM32, ABM43, ABM54

# Snow model
export AbstractSnowModel, NoSnow, SnowModel, SnowState

# Apparent heat capacity (latent-heat-of-fusion treatment in snow)
export AbstractApparentHeatCapacity, BonacinaStep, TanhSmoothed, Gaussian, WestermannSigmoid, apparent_heat_capacity

# Soil thermal model
export AbstractSoilProperties, CampbelldeVriesSoilProperties

# Soil hydraulics
export AbstractSoilHydraulicsModel, CampbellSoilHydraulics

# Soil moisture strategy
export AbstractSoilMoistureStrategy, PrescribedSoilMoisture, DynamicSoilMoisture

# Convergence strategies
export AbstractSoilTemperatureConvergence, FixedSoilTemperatureIterations, SoilTemperatureConvergenceTolerance

# Time modes
export AbstractTimeMode, NonConsecutiveDayMode, ConsecutiveDayMode

# Diffuse fraction models
export AbstractDiffuseFractionModel, ErbsDiffuseFraction

# Atmospheric radiation models (clear-sky longwave building block)
export AbstractAtmosphericRadiationModel, SwinbankAtmosphericRadiation, CampbellNormanAtmosphericRadiation, atmospheric_radiation

# Sunshine fraction models (Ångström building block)
export AbstractSunshineFractionModel, Angstrom, sunshine_fraction

# Longwave budget algorithms
export AbstractLongwaveModel, ViewFactorLongwave

# Shortwave budget algorithms
export AbstractShortwaveModel, AngstromMaxwellShortwave

# Surface evaporation models
export AbstractEvaporationModel, BulkTransferEvaporation, surface_convection_evaporation

# Soil energy balance models
export SoilHeatTransportModel, SoilHeatTransport1D

# Soil freezing models
export SoilPhaseTransitionModel, PhaseTransitionLatentHeat, allocate_phase_transition

# Rainfall schedule
export AbstractRainfallSchedule, DailyRainfall, HourlyRainfall, is_hourly

export AbstractEnvironment, MonthlyMinMaxEnvironment, DailyMinMaxEnvironment, DailyTimeseries, HourlyTimeseries
export Site, AbstractSite
export AbstractBoundaryLayerModel, MoninObukhov
export Forcing, AtmosphericProfile


export shortwave_radiation!, longwave_radiation, precompute_longwave_sky

export calc_convection

export soil_properties, soil_properties!, allocate_soil_properties

export soil_energy_balance, evaporation, soil_water_balance!, phase_transition

export example_site, example_monthly_weather,
    example_daily_environment, example_hourly_environment,
    example_soil_hydraulic_model, example_soil_properties_model, example_microclimate_problem

import CommonSolve: solve, solve!, init
export solve, solve!, init, reinit!


include("constants.jl")

# Soil properties
include("soil_properties/abstract.jl")
include("soil_properties/campbell_devries.jl")

# Soil moisture strategy
include("soil_moisture/abstract.jl")
include("soil_moisture/prescribed.jl")
include("soil_moisture/dynamic.jl")

# Surface evaporation
include("evaporation/abstract.jl")
include("evaporation/bulk_transfer.jl")

# Soil hydraulics
include("soil_hydraulics/abstract.jl")
include("soil_hydraulics/campbell.jl")

# Site
include("site/abstract.jl")
include("site/site.jl")

# Boundary layer
include("boundary_layer/abstract.jl")
include("boundary_layer/monin_obukhov.jl")

# Environment inputs
include("inputs/abstract_environment.jl")
include("inputs/monthly_minmax.jl")
include("inputs/daily_minmax.jl")
include("inputs/daily_timeseries.jl")
include("inputs/hourly_timeseries.jl")

# Soil temperature convergence
include("soil_temperature_convergence/abstract.jl")
include("soil_temperature_convergence/fixed.jl")
include("soil_temperature_convergence/tolerance.jl")

# Time mode
include("time_mode/abstract.jl")
include("time_mode/non_consecutive.jl")
include("time_mode/consecutive.jl")

# Radiation building blocks and surface budgets
include("radiation/diffuse/abstract.jl")
include("radiation/diffuse/erbs.jl")

include("radiation/atmospheric/abstract.jl")
include("radiation/atmospheric/swinbank.jl")
include("radiation/atmospheric/campbell_norman.jl")

include("radiation/sunshine_fraction/abstract.jl")
include("radiation/sunshine_fraction/angstrom.jl")

include("radiation/longwave/abstract.jl")
include("radiation/longwave/view_factor.jl")

include("radiation/shortwave/abstract.jl")
include("radiation/shortwave/angstrom_maxwell.jl")

# Rainfall schedule
include("rainfall/abstract.jl")
include("rainfall/daily_rainfall.jl")
include("rainfall/hourly_rainfall.jl")

# Outputs and forcing (depend on AbstractEnvironment / AbstractInterpolation)
include("outputs.jl")
include("forcing.jl")

# Algorithms
include("interpolation.jl")

# Soil balance: freezing correction is referenced as a default by the energy
# ODE, so freezing/ must come before energy/.
include("soil_balance/freezing/abstract.jl")
include("soil_balance/freezing/latent_heat.jl")
include("soil_balance/energy/abstract.jl")
include("soil_balance/energy/heat_transport_1d.jl")

# Apparent heat capacity (latent-heat treatment near phase change; used by snow)
include("apparent_heat_capacity/abstract.jl")
include("apparent_heat_capacity/bonacina.jl")
include("apparent_heat_capacity/tanh.jl")
include("apparent_heat_capacity/gaussian.jl")
include("apparent_heat_capacity/westermann.jl")

# Snow models (one file per variant)
include("snow/abstract.jl")
include("snow/no_snow.jl")
include("snow/snow_model.jl")

# Top-level config, parameters, problem, state, buffers, and cache types
include("types.jl")
include("simulation.jl")

end
