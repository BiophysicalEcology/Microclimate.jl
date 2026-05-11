"""
    DailyRainfall()

The day's total rainfall is read from `environment_daily.rainfall` and
applied at the solar-midnight hour. The hourly forcing's `rainfall` array
is ignored.
"""
struct DailyRainfall <: AbstractRainfallSchedule end

is_hourly(::DailyRainfall) = false
