"""
    DynamicSoilMoisture(; moisture_tolerance=1e-6u"kg/m^2/s",
                          moisture_max_iterations=500,
                          moisture_timestep=360.0u"s")

Run the Campbell soil water balance solver to dynamically compute soil moisture
at each hourly timestep. Tracks the evolving surface soil wetness fraction
internally.

Solver tuning:
- `moisture_tolerance` — maximum mass-balance residual
- `moisture_max_iterations` — maximum iterations of the moisture solver
- `moisture_timestep` — timestep of the moisture solver (≤ 1 hour)
"""
mutable struct DynamicSoilMoisture{MT,MTS} <: AbstractSoilMoistureStrategy
    soil_wetness::Float64
    moisture_tolerance::MT
    moisture_max_iterations::Int
    moisture_timestep::MTS
end
function DynamicSoilMoisture(;
    moisture_tolerance = 1e-6u"kg/m^2/s",
    moisture_max_iterations = 500,
    moisture_timestep = 360.0u"s",
)
    DynamicSoilMoisture(0.0, moisture_tolerance, moisture_max_iterations, moisture_timestep)
end

get_soil_wetness(mode::DynamicSoilMoisture, _) = mode.soil_wetness
init_soil_wetness!(mode::DynamicSoilMoisture) = (mode.soil_wetness = 0.0; nothing)

function initialise_soil_humidity!(::DynamicSoilMoisture, output, soil_water_potential, T0)
    water_molar_mass = 0.01801528u"kg/mol"
    output.soil_humidity[1, :] = clamp.(exp.(water_molar_mass .* soil_water_potential ./ (R .* T0)), 0, 1)
end
function reset_day_soil_moisture!(::DynamicSoilMoisture, soil_moisture, initial_soil_moisture, day_index)
    soil_moisture .= initial_soil_moisture
end
