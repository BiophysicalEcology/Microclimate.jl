"""
    AbstractEnvironment

Supertype for *input* environment forcings (`MonthlyMinMaxEnvironment`,
`DailyMinMaxEnvironment`, `DailyTimeseries`, `HourlyTimeseries`). Output
data structures (`MicroResult`) deliberately do NOT subtype this — inputs
and outputs are different concepts.
"""
abstract type AbstractEnvironment end

"""
    consecutive_days(env::AbstractEnvironment) -> Bool

Whether successive day-indices are adjacent calendar days to trigger
appropriate `DielCurve` shape algorithm.
"""
consecutive_days(::AbstractEnvironment) = false
