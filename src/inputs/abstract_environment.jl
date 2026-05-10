"""
    AbstractEnvironment

Supertype for *input* environment forcings (`MonthlyMinMaxEnvironment`,
`DailyMinMaxEnvironment`, `DailyTimeseries`, `HourlyTimeseries`). Output
data structures (`MicroResult`) deliberately do NOT subtype this — inputs
and outputs are different concepts.
"""
abstract type AbstractEnvironment end
