"""
    RootFindLeafTemperature()

Exact root-find (`HeatExchange.zbrent`) of `leaf_heat_balance` each call.
More expensive than [`LinearizedLeafTemperature`](@ref), but doesn't depend
on the quality of a supplied guess.
"""
struct RootFindLeafTemperature <: AbstractLeafTemperatureSolver end

function leaf_temperature(::RootFindLeafTemperature, absorbed_radiation, air_temperature, relative_humidity,
    wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body;
    leaf_area=1.0u"m^2", leaf_temperature_guess=air_temperature,
    bracket=(air_temperature - 30.0u"K", air_temperature + 40.0u"K"),
)
    residual(T) = ustrip(u"W", leaf_heat_balance(T * u"K", absorbed_radiation, air_temperature, relative_humidity,
        wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body, leaf_area).net)
    lo, hi = ustrip(u"K", bracket[1]), ustrip(u"K", bracket[2])
    solved = try
        zbrent(residual, lo, hi, 1e-3)
    catch
        ustrip(u"K", leaf_temperature_guess)
    end
    return solved * u"K"
end
