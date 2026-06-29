"""
    DailyMinMaxEnvironment(; forcings)

Per-day analogue of [`MonthlyMinMaxEnvironment`](@ref) for consecutive-day
simulations (ERA5, station data, etc.). Each `Forcing`'s value series carries
one entry per actual calendar day. Passing this to `simulate_microclimate`
automatically sets `daily=true` so consecutive days inherit soil state and
iterate once.

`forcings` is a `NamedTuple` keyed by output variable, each a [`Forcing`](@ref);
`minmax_forcings` builds the two-anchor min/max case.
"""
struct DailyMinMaxEnvironment{NT<:NamedTuple} <: AbstractEnvironment
    forcings::NT
end
DailyMinMaxEnvironment(; forcings) = DailyMinMaxEnvironment(forcings)
