"""
    PrescribedSoilMoisture(; precomputed_soil_moisture=nothing)

Use prescribed soil moisture values rather than running a dynamic soil moisture solver.
Soil wetness is taken from the daily environment timeseries.

When `precomputed_soil_moisture` is provided as a matrix (ndepths × ndays), those values
override `initial_soil_moisture` at the start of each day (monthly-representative mode only).
"""
struct PrescribedSoilMoisture{PSM} <: AbstractSoilMoistureStrategy
    precomputed_soil_moisture::PSM
end
PrescribedSoilMoisture(; precomputed_soil_moisture=nothing) =
    PrescribedSoilMoisture(precomputed_soil_moisture)

get_soil_wetness(::PrescribedSoilMoisture, environment_instant) = environment_instant.soil_wetness
init_soil_wetness!(::PrescribedSoilMoisture) = nothing
initialise_soil_humidity!(::PrescribedSoilMoisture, output, soil_water_potential, T0) = nothing

function reset_day_soil_moisture!(mode::PrescribedSoilMoisture, soil_moisture, initial_soil_moisture, day_index)
    if !isnothing(mode.precomputed_soil_moisture)
        soil_moisture .= mode.precomputed_soil_moisture[:, day_index]
    else
        soil_moisture .= initial_soil_moisture
    end
end
