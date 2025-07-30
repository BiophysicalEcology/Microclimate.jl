module Microclimate

function __init__()\
    Unitful.register(Microclimate)
end

using OrdinaryDiffEq, Interpolations, Statistics, Dates

using Unitful, UnitfulMoles, ModelParameters

using Unitful: °, rad, °C#, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R

export air_pressure, vapour_pressure, wet_air, dry_air, λ_evap, phase_transition!

export sinec!, vsine, hourly_vars

export hour_angle, solar_geometry, check_skylight, elev_corr, gamma, dchxy, dexpi, solrad

export get_profile

export soil_properties

export init_soillayers, MicroInputs

export soil_energy_balance!, evap, soil_water_balance

export MicroParams, MicroForcing, MicroInput

include("landscape.jl")   
include("interpolation.jl")
include("fluid_properties.jl")
include("soil_properties.jl")
include("radiation.jl")
include("boundary_layer.jl")
include("soil_balance.jl")
include("simulation.jl")

end
