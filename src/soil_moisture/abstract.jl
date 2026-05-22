abstract type AbstractSoilMoistureStrategy end

"""
    get_soil_wetness(strategy, environment_instant)

Surface wetness (0–1) at the current step.
"""
function get_soil_wetness end

"""
    init_soil_wetness!(strategy)

Reset any mutable wetness state held by the strategy at the start of a run.
"""
function init_soil_wetness! end

"""
    initialise_soil_humidity!(strategy, output, soil_water_potential, T0)

Populate the initial soil-humidity row of `output` from the first-hour soil
water potential and temperature.
"""
function initialise_soil_humidity! end

"""
    reset_day_soil_moisture!(strategy, soil_moisture, initial_soil_moisture, day_index)

Restore per-day soil moisture for time modes that treat days independently.
"""
function reset_day_soil_moisture! end
