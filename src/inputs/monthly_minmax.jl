# TODO: this should be more generic.
# We could possible make a field type that is either interpolated or indexed
# so we just mix min-max fields with e.g. daily fields in a single environment object
@kwdef struct MonthlyMinMaxEnvironment{AT,W,H,C,M} <: AbstractEnvironment
    reference_temperature_min::AT
    reference_temperature_max::AT
    reference_wind_min::W
    reference_wind_max::W
    reference_humidity_min::H
    reference_humidity_max::H
    cloud_min::C
    cloud_max::C
    minima_times::M
    maxima_times::M
end

"""
    example_monthly_weather

Example monthly weather input for the middle day of each month for Madison Wisconsin, USA, originally extracted from the CRU CL v. 2.0 dataset.
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
    # Per-variable hour offsets indexed [air_temp, wind, humidity, cloud].
    # Air temp & wind: minima offset from sunrise, maxima offset from solar noon.
    # Humidity & cloud: minima offset from solar noon, maxima offset from sunrise.
    minima_times = [0, 0, 1, 1],
    maxima_times = [1, 1, 0, 0],
)
    MonthlyMinMaxEnvironment(;
        reference_temperature_min, reference_temperature_max,
        reference_wind_min, reference_wind_max,
        reference_humidity_min, reference_humidity_max,
        cloud_min, cloud_max,
        minima_times, maxima_times,
    )
end
