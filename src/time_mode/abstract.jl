abstract type AbstractTimeMode end

"""
    independent_days(time_mode) -> Bool

Whether each day is solved independently (no carry-over of state).
"""
function independent_days end

"""
    is_reset_day(time_mode, day_index) -> Bool

Whether this day starts from reset initial state rather than continuing.
"""
function is_reset_day end

"""
    reset_phase_per_iter(time_mode) -> Bool

Whether the soil phase-change accumulator should reset each iteration.
"""
function reset_phase_per_iter end

"""
    reset_moisture_per_day(time_mode) -> Bool

Whether soil moisture should reset to its initial value at the start of each day.
"""
function reset_moisture_per_day end

"""
    reset_snow_per_day(time_mode) -> Bool

Whether snow scratch should reset at the start of each day.
"""
function reset_snow_per_day end

"""
    iter_resets_T(time_mode) -> Bool

Whether soil temperature is reset between within-day iterations.
"""
function iter_resets_T end

"""
    iterations_for_day(time_mode, convergence, day_index) -> Int

Number of iterations to run for this day under the given convergence strategy.
"""
function iterations_for_day end
