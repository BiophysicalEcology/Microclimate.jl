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

Whether successive day-indices are genuinely adjacent calendar days (so a
midnight-spanning `DielCurve` shape can use index `i+1`'s values for its far
endpoint) rather than non-adjacent representative days (e.g.
`MonthlyMinMaxEnvironment`), where index `i+1` isn't tomorrow.
"""
consecutive_days(::AbstractEnvironment) = false
