module Microclimate

using ConstructionBase
using Setfield: @set, @set!
using Interpolations, Statistics, Dates
using SciMLBase, OrdinaryDiffEqTsit5
using Unitful, UnitfulMoles
using ModelParameters, DelimitedFiles
using SpecialFunctions, StaticArrays

using FluidProperties: atmospheric_pressure, wet_air_properties, dry_air_properties, vapour_pressure 
using FluidProperties: enthalpy_of_vaporisation, molar_enthalpy_of_vaporisation, water_properties
using FluidProperties: g_n, σ, atm, R
using Unitful: °, rad, °C
using Interpolations: AbstractInterpolation
using SolarRadiation


export MicroProblem

export CampbelldeVriesSoilThermal, SoilMoistureModel

export MonthlyMinMaxEnvironment, DailyTimeseries, HourlyTimeseries, MicroTerrain

export daily_cycle_sine_exponential, daily_cycle_linear, hourly_from_min_max

export cloud_adjust_radiation, longwave_radiation

export atmospheric_surface_profile, calc_convection

export soil_properties, soil_properties!, allocate_soil_properties

export soil_energy_balance, evaporation, soil_water_balance!, phase_transition

export example_micro_terrain, example_monthly_weather, example_daily_environmental, example_soil_moisture_model, example_soil_thermal_parameters, example_microclimate_problem

# TODO replace this with CommonSolve.jl
export solve

export cold_air_pooling, ColdAirPoolingMethod, ColdAirFlow
export surface_water_flow, surface_water_event, SurfaceWaterMethod, SurfaceWaterFlow

# Snow models
export SnowModel, NoSnow, DegreeDaySnow, Snow17, UtahEnergyBalance, KearneySnow
export SnowAccumulation, NoAccumulation, ThresholdAccumulation, LinearTransitionAccumulation
export SnowState, SnowForcing
export snow_accumulation, snow_melt, update_snow_state
export snow_thermal_conductivity, snow_specific_heat, snow_density, snow_albedo, snow_depth
export seasonal_melt_factor, snow_energy_balance, snow_temperature, snow_refreeze
# Formula structs
export ThermalConductivityFormula, Djachkova
export SnowDensityFormula, SturmSnowDensity, CompactionSnowDensity, SimpleMixingDensity, AndersonSnowDensity, AndersonDensityEvolution
export AlbedoFormula, NoAlbedo, AndersonAlbedo, ExponentialSnowAlbedo
export SnowTemperatureFormula, NoTemperatureEvolution, RelaxationTemperature, ColdContentTemperature, EnergyContentTemperature

# Spatial snow
export SnowRedistributionMethod, WindDrivenSnowRedistribution, GravitationalSnowRedistribution
export snow_redistribution, update_snow_from_precipitation!, update_snow_melt!

# Spatial simulation
export SpatialMicroState, SpatialMicroTerrain, SpatialMicroProblem
export surface_temperature, surface_moisture

include("constants.jl")
include("snow.jl")
include("spatial/cold_air.jl")
include("spatial/surface_water.jl")
include("spatial/state.jl")
include("spatial/terrain.jl")
include("spatial/snow.jl")
include("spatial/coupling.jl")
include("spatial/solve.jl")
include("landscape.jl")   
include("interpolation.jl")
include("soil_properties.jl")
include("radiation.jl")
include("boundary_layer.jl")
include("soil_balance.jl")
include("simulation.jl")

end
