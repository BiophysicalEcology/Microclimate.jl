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
