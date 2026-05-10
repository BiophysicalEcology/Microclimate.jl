"""
    DailyMinMaxEnvironment

Per-day analogue of `MonthlyMinMaxEnvironment` for consecutive-day simulations
(ERA5, station data, etc.).  Each entry corresponds to one actual calendar day.
Passing this to `simulate_microclimate` automatically sets `daily=true` so that
consecutive days inherit soil state and iterate once.
"""
@kwdef struct DailyMinMaxEnvironment{AT,W,H,C,M}
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
