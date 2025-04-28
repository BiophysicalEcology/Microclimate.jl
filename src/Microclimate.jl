module Microclimate

function __init__()\
    Unitful.register(Microclimate)
end

using DifferentialEquations, Interpolations, Statistics, Dates

using Unitful, UnitfulMoles, ModelParameters

using Unitful: °, rad, °C#, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R

export air_pressure, vapour_pressure, wet_air, dry_air

export MicroParams

export sinec!, vsine, hourly_vars

export hour_angle, solar_geometry, check_skylight, elev_corr, dchxy, dexpi, solrad

export get_profile

export soil_properties

export soil_energy_balance!, evap

export MicroParams, MicroForcing, MicroInput

include("air_properties.jl")
include("landscape.jl")      
include("interpolation.jl")
include("solar_radiation.jl")
include("boundary_layer.jl")
include("soil_properties.jl")
include("soil_balance.jl")

end
