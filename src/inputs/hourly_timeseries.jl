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
