abstract type AbstractRainfallSchedule end

"""
    is_hourly(::AbstractRainfallSchedule)::Bool

Whether rainfall arrives one value per hour rather than as a daily total.
"""
function is_hourly end

"""
    current_rainfall(schedule; environment_instant, environment_hourly, step, i, midnight_i)

Rainfall (kg/m²) entering the surface at the given step. `step` is the
absolute hour index, `i` is the within-day hour, `midnight_i` is the hour
at which the daily total is released for `DailyRainfall`.
"""
function current_rainfall end
