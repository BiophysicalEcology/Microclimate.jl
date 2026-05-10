"""
    DynamicSoilMoisture()

Run the Campbell soil water balance solver to dynamically compute soil moisture
at each hourly timestep. Tracks the evolving surface soil wetness fraction
internally.
"""
mutable struct DynamicSoilMoisture <: AbstractSoilMoistureMode
    soil_wetness::Float64
end
DynamicSoilMoisture() = DynamicSoilMoisture(0.0)

get_soil_wetness(mode::DynamicSoilMoisture, _) = mode.soil_wetness
init_soil_wetness!(mode::DynamicSoilMoisture) = (mode.soil_wetness = 0.0; nothing)

function initialise_soil_humidity!(::DynamicSoilMoisture, output, soil_water_potential, T0)
    water_molar_mass = 0.01801528u"kg/mol"
    output.soil_humidity[1, :] = clamp.(exp.(water_molar_mass .* soil_water_potential ./ (R .* T0)), 0, 1)
end
function reset_day_soil_moisture!(::DynamicSoilMoisture, soil_moisture, initial_soil_moisture, day_index)
    soil_moisture .= initial_soil_moisture
end
