"""
    DailyRainfall(; spread_hours=1)

The day's total rainfall (`environment_daily.rainfall`) is split evenly
across `spread_hours` hours starting at solar midnight (24-hour day);
`spread_hours=1` is the original single-hour dump. Ignores hourly `rainfall`.
"""
@kwdef struct DailyRainfall{SH} <: AbstractRainfallSchedule
    spread_hours::SH = 1
end

is_hourly(::DailyRainfall) = false

current_rainfall(schedule::DailyRainfall; environment_instant, environment_hourly, step, i, midnight_i) =
    mod(i - midnight_i, 24) < schedule.spread_hours ? environment_instant.rainfall / schedule.spread_hours : 0.0u"kg/m^2"
