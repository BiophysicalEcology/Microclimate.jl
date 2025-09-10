module Microclimate

using FluidProperties: wet_air, dry_air, vapour_pressure, get_λ_evap, waterprop, get_pressure

using OrdinaryDiffEq, Interpolations, Statistics, Dates

using Unitful, UnitfulMoles, ModelParameters, DelimitedFiles

using Unitful: °, rad, °C#, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R

using SpecialFunctions, StaticArrays

export gads

export sinec!, vsine, hourly_vars

export hour_angle, solar_geometry, elev_corr, dchxy, solrad, cloud_adjust_radiation, get_longwave, init_dchxy_buffers

export get_profile

export soil_properties

export init_soillayers, init_moistlayers, MicroInputs

export soil_energy_balance!, evap, soil_water_balance, phase_transition

export MicroParams, MicroForcing

export runmicro

export load_CRU_CL_v2

include("constants.jl")
include("gads.jl")
include("landscape.jl")   
include("interpolation.jl")
include("soil_properties.jl")
include("radiation.jl")
include("boundary_layer.jl")
include("soil_balance.jl")
include("simulation.jl")
include("weather_and_climate.jl")

function __init__()
    Unitful.register(Microclimate)
end

end
