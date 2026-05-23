"""
    Forcing

Per-day in-place interpolated forcings driving the soil-energy ODE. Each
field is a `ScaledInterpolation` over a 25-element coefficient buffer that
gets refilled per day-iter via `update_forcing_day!`.
"""
@kwdef struct Forcing{
    S<:AbstractInterpolation,ZE<:AbstractInterpolation,ZS<:AbstractInterpolation,T<:AbstractInterpolation,
    V<:AbstractInterpolation,RH<:AbstractInterpolation,CL<:AbstractInterpolation,P<:AbstractInterpolation,
}
    solar::S
    zenith::ZE
    slope_zenith::ZS
    temperature::T
    wind::V
    humidity::RH
    cloud::CL
    pressure::P
end

function interpolate_forcings(f::Forcing, t)
    t_m = ustrip(u"minute", t)
    return (;
        atmospheric_pressure = f.pressure(t_m),
        air_temperature = f.temperature(t_m),
        wind_speed = f.wind(t_m),
        zenith_angle = min(90.0u"°", u"°"(round(f.zenith(t_m), digits=3))),
        solar_radiation = max(0.0u"W/m^2", f.solar(t_m)),
        cloud_cover = clamp(f.cloud(t_m), 0.0, 1.0),
        relative_humidity = clamp(f.humidity(t_m), 0.0, 1.0),
        slope_zenith_angle = min(90.0u"°", f.slope_zenith(t_m)),
    )
end
