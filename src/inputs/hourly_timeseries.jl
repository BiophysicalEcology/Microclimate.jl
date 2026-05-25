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
