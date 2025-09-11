module Microclimate

using OrdinaryDiffEq, Interpolations, Statistics, Dates
using Unitful, UnitfulMoles, ModelParameters, DelimitedFiles
using SpecialFunctions, StaticArrays

using FluidProperties: wet_air, dry_air, vapour_pressure, get_λ_evap, waterprop, get_pressure
using Unitful: °, rad, °C#, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R
using Interpolations: AbstractInterpolation


export sinec!, vsine, hourly_vars

export hour_angle, solar_geometry, elev_corr, dchxy, solrad, cloud_adjust_radiation, get_longwave, init_dchxy_buffers

export get_profile

export soil_properties

export init_soillayers, init_moistlayers, MicroInputs

export soil_energy_balance!, evap, soil_water_balance, phase_transition

export MicroParams, MicroForcing

export runmicro


include("constants.jl")
include("landscape.jl")   
include("interpolation.jl")
include("soil_properties.jl")
include("radiation.jl")
include("boundary_layer.jl")
include("soil_balance.jl")
include("simulation.jl")


function __init__()
    Unitful.register(Microclimate)
end

end
