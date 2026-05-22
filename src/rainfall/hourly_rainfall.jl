"""
    HourlyRainfall()

Rainfall is read hourly from `environment_hourly.rainfall`. The daily-total
field on `environment_daily` is ignored.
"""
struct HourlyRainfall <: AbstractRainfallSchedule end

is_hourly(::HourlyRainfall) = true

current_rainfall(::HourlyRainfall; environment_instant, environment_hourly, step, i, midnight_i) =
    environment_hourly.rainfall[step]
