"""
    HourlyTimeseries(; pressure, reference_temperature, reference_humidity,
                       reference_wind_speed, global_radiation, longwave_radiation,
                       cloud_cover, rainfall, zenith_angle)

Directly-supplied hourly forcing, one value per `(day, hour)` combination
(length `length(days) * length(hours)`), for driving the model from 
weather-station or reanalysis data instead of interpolated daily min/max.

- `pressure` — atmospheric pressure (Pa)
- `reference_temperature`, `reference_humidity`, `reference_wind_speed` — air
  temperature (°C), relative humidity (0–1) and wind speed (m/s) at the
  site's reference height
- `global_radiation` — global horizontal solar irradiance (W/m²)
- `longwave_radiation` — incoming sky longwave irradiance (W/m²), if
  measured; when supplied it is used in place of the cloud-cover-derived
  estimate (see [`precompute_longwave_sky`](@ref))
- `cloud_cover` — fractional cloud cover (0–1)
- `rainfall` — hourly rainfall total (kg/m²)
- `zenith_angle` — solar zenith angle (°), if supplied directly rather than
  computed from the site and time
"""
@kwdef struct HourlyTimeseries{P,RT,RH,RWS,GR,LW,CC,R,ZA} <: AbstractEnvironment
    pressure::P
    reference_temperature::RT
    reference_humidity::RH
    reference_wind_speed::RWS
    global_radiation::GR
    longwave_radiation::LW
    cloud_cover::CC
    rainfall::R
    zenith_angle::ZA
end

"""
    example_hourly_environment(days=DEFAULT_DAYS, hours=DEFAULT_HOURS; kwargs...)

Example [`HourlyTimeseries`](@ref) with only `pressure` populated (from
`elevation`); every other field defaults to `nothing`. Note that
`reference_temperature`/`reference_humidity`/`reference_wind_speed`/`cloud_cover` are
only read from here when `MicroInputs.environment_minmax` is `nothing` — otherwise the
monthly/daily min-max forcing takes priority for those fields regardless of what's
supplied here (`global_radiation` is the exception: it's used whenever non-`nothing`).
 `kwargs` override any field.
"""
function example_hourly_environment(days=DEFAULT_DAYS, hours=DEFAULT_HOURS;
    elevation = 226.0u"m",
    pressure = fill(atmospheric_pressure(elevation), length(days) * length(hours)),
    reference_temperature = nothing,
    reference_humidity = nothing,
    reference_wind_speed = nothing,
    global_radiation = nothing,
    cloud_cover = nothing,
    rainfall = nothing,
    zenith_angle = nothing,
    longwave_radiation = nothing,
)
    HourlyTimeseries(;
        pressure, reference_temperature, reference_humidity, reference_wind_speed,
        global_radiation, cloud_cover, rainfall, zenith_angle, longwave_radiation,
    )
end
