module Microclimate

using ConstructionBase
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

export sine_exponential!, vsine, hourly_vars

export cloud_adjust_radiation, longwave_radiation

export atmospheric_surface_profile, calc_convection

export soil_properties, soil_properties!, allocate_soil_properties

export soil_energy_balance, evaporation, soil_water_balance!, phase_transition

export example_micro_terrain, example_monthly_weather, example_daily_environmental, example_soil_moisture_model, example_soil_thermal_parameters, example_microclimate_problem

# TODO replace this with CommonSolve.jl
export solve


include("constants.jl")
include("landscape.jl")   
include("interpolation.jl")
include("soil_properties.jl")
include("radiation.jl")
include("boundary_layer.jl")
include("soil_balance.jl")
include("simulation.jl")

end
