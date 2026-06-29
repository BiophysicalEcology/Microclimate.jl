# A monthly min/max (or any timed-anchor) environment: one representative day
# per month, each variable a `Forcing` (declarative `DielCurve` + per-day data).
"""
    MonthlyMinMaxEnvironment(; forcings)

Forcing for non-consecutive representative days (Fortran "monthly mode").
`forcings` is a `NamedTuple` keyed by output variable
(`reference_temperature`, `reference_wind_speed`, `reference_humidity`,
`cloud_cover`, …); each is a [`Forcing`](@ref). `minmax_forcings` builds the
classic min/max case — see `example_monthly_weather`.
"""
struct MonthlyMinMaxEnvironment{NT<:NamedTuple} <: AbstractEnvironment
    forcings::NT
end
MonthlyMinMaxEnvironment(; forcings) = MonthlyMinMaxEnvironment(forcings)

"""
    example_monthly_weather

Example monthly weather input for the middle day of each month for Madison
Wisconsin, USA, originally extracted from the CRU CL v. 2.0 dataset.
"""
function example_monthly_weather(;
    reference_temperature_min = [-14.3, -12.1, -5.1, 1.2, 6.9, 12.3, 15.2, 13.6, 8.9, 3, -3.2, -10.6]u"°C",
    reference_temperature_max = [-3.2, 0.1, 6.8, 14.6, 21.3, 26.4, 29, 27.7, 23.3, 16.6, 7.8, -0.4]u"°C",
    reference_wind_min = [4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8] * 0.1u"m/s",
    reference_wind_max = [4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8]u"m/s",
    reference_humidity_min = [50.2, 48.4, 48.7, 40.8, 40, 42.1, 45.5, 47.3, 47.6, 45, 51.3, 52.8] ./ 100.0,
    reference_humidity_max = fill(1.0, 12),
    cloud_min = [50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1] ./ 100.0,
    cloud_max = [50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1] ./ 100.0,
)
    MonthlyMinMaxEnvironment(; forcings = minmax_forcings(;
        reference_temperature_min, reference_temperature_max,
        reference_wind_min, reference_wind_max,
        reference_humidity_min, reference_humidity_max,
        cloud_min, cloud_max,
    ))
end
