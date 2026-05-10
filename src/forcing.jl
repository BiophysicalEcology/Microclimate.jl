@kwdef struct MicroForcing{
    S<:AbstractInterpolation,ZE<:AbstractInterpolation,ZS<:AbstractInterpolation,T<:AbstractInterpolation,
    V<:AbstractInterpolation,RH<:AbstractInterpolation,CL<:AbstractInterpolation,P<:AbstractInterpolation,
}
    interpolate_solar::S
    interpolate_zenith::ZE
    interpolate_slope_zenith::ZS
    interpolate_temperature::T
    interpolate_wind::V
    interpolate_humidity::RH
    interpolate_cloud::CL
    interpolate_pressure::P
end
