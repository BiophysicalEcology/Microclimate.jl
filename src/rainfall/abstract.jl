abstract type AbstractRainfallSchedule end

"""
    is_hourly(::AbstractRainfallSchedule)::Bool

Whether rainfall arrives one value per hour rather than as a daily total.
Implementations dispatch the day-loop's rainfall handling and the snow
model's snowfall accounting.
"""
function is_hourly end
