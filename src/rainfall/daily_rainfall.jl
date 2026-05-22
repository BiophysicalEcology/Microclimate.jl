"""
    DailyRainfall()

The day's total rainfall is read from `environment_daily.rainfall` and
applied at the solar-midnight hour. The hourly forcing's `rainfall` array
is ignored.
"""
struct DailyRainfall <: AbstractRainfallSchedule end

is_hourly(::DailyRainfall) = false

current_rainfall(::DailyRainfall; environment_instant, environment_hourly, step, i, midnight_i) =
    i == midnight_i ? environment_instant.rainfall : 0.0u"kg/m^2"
