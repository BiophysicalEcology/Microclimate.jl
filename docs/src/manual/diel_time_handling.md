# Diel curves and time handling

Driving weather data is often only available as daily minimum and maximum values —
air temperature, humidity, wind speed, cloud cover. However, the model needs hourly values, so
these daily extremes have to be turned into a within-day cycle: for example air temperature is
often modelled as a sine wave from the daily minimum (near sunrise) up through the maximum, and an
exponential decay overnight back down; humidity, wind and cloud might be modelled piecewise-linear.
The timing of the minima and maxima can be anchored tosolar reference points (relative to dawn or midday).

## A declarative curve model

Rather than hard-coding those shapes, Microclimate.jl expresses them with a small
declarative vocabulary:

- [`TimeOfDay`](@ref) — a within-day instant: [`Sunrise`](@ref), [`Sunset`](@ref) and
  [`Midday`](@ref) are solar events, resolved from that day's solar geometry;
  [`Midnight`](@ref) and [`ClockTime`](@ref) are fixed clock times, independent of the
  sun. All four take an optional hour offset.
- [`Shape`](@ref) — how a curve travels between two `TimeOfDay`s: [`Sine`](@ref)
  (space-filling — pinned by its trough and peak, fills whatever span the other shapes
  leave), [`Decay`](@ref) and [`Linear`](@ref) (bounded — cover exactly their own
  `[from, to)` span).
- [`DielCurve`](@ref) — a tuple of shapes covering the full day, checked at
  construction for gaps and overlaps.
- [`DielForcing`](@ref) — a `DielCurve` paired with the actual per-day min/max (or
  other) value series.

[`minmax_forcings`](@ref) builds this standard curve set — temperature a sine to an
afternoon peak with an overnight decay; wind, humidity and cloud piecewise-linear
between their daily extremes and the daily mean at true midnight:

```@example diel
using Microclimate, Unitful

forcings = minmax_forcings(;
    reference_temperature_min = [-14.3, -12.1]u"°C",
    reference_temperature_max = [-3.2, 0.1]u"°C",
    reference_wind_speed_min = [0.49, 0.48]u"m/s",
    reference_wind_speed_max = [4.9, 4.8]u"m/s",
    reference_humidity_min = [0.502, 0.484],
    reference_humidity_max = [1.0, 1.0],
    cloud_cover_min = [0.503, 0.47],
    cloud_cover_max = [0.503, 0.47],
)
keys(forcings)
```

This declarative structure — write down *when* the curve changes shape and *how*, not
a fixed algorithm — is why timezone, leap-year and daylight-saving handling can be
automatic rather than something the user manages by hand: every `TimeOfDay` resolves
against the actual solar geometry and calendar of the day being simulated, so the same
`DielCurve` definition is correct regardless of location or year. `Midnight` and
`ClockTime` are the two fixed anchors that never move with solar geometry, by design.

## Non-standard curves and derived quantities

[`ForcingSpec`](@ref) is the structure-only half of a `DielForcing` (a `DielCurve`
plus the names of the per-day quantities that feed it), letting a whole model be
declared once and bound to data later. [`Derived`](@ref) computes a quantity from
other named quantities rather than interpolating a curve directly — for example
[`RelativeHumidityFromVapourPressureAndTemperature`](@ref) and its inverse
[`VapourPressureFromRelativeHumidityAndTemperature`](@ref), for models that carry
vapour pressure rather than relative humidity as the forcing variable.

## Time modes

Whether a day's initial soil-temperature profile carries forward from the previous
day, or starts fresh from the daily mean, is controlled separately, by
[`AbstractTimeMode`](@ref) on [`MicroProblem`](@ref):

- [`NonConsecutiveDayMode`](@ref) — independent representative days (e.g. one day per
  month), each iterated to a steady-periodic solution from a uniform initial profile.
- [`ConsecutiveDayMode`](@ref) — a continuous run where state carries from one day to
  the next, with an optional `spinup_first_day`.

## References

See the [References](references.md) page for the full bibliography.
