"""
    ConsecutiveDayMode(; spinup_first_day=false)

Consecutive-day simulation where the end state of each day becomes the
initial condition for the next. Only one iteration per day is performed.
When `spinup_first_day=true`, the first day is iterated multiple times to
establish a quasi-equilibrium initial state.
"""
struct ConsecutiveDayMode <: AbstractTimeMode
    spinup_first_day::Bool
end
ConsecutiveDayMode(; spinup_first_day=false) = ConsecutiveDayMode(spinup_first_day)

independent_days(::ConsecutiveDayMode) = false
is_reset_day(::ConsecutiveDayMode, j::Int) = false
reset_phase_per_iter(::ConsecutiveDayMode) = false
reset_moisture_per_day(::ConsecutiveDayMode) = false
reset_snow_per_day(::ConsecutiveDayMode) = false
iter_resets_T(::ConsecutiveDayMode) = true

function iterations_for_day(mode::ConsecutiveDayMode, convergence::AbstractSoilTemperatureConvergence, day_index)
    return (mode.spinup_first_day && day_index == 1) ? max_iterations(convergence) : 1
end
