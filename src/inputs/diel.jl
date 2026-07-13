# Declarative within-day curve model.
#
# A forcing variable's daily curve is defined by *shapes* over *times* — purely
# structural, no data. Per-slice values are supplied later, against this fixed
# definition. Three axes, kept separate:
#
#   Time  — when in the day. Either movable, resolved from the day's solar
#           geometry (Sunrise, Sunset, Midday), or fixed to a clock hour
#           regardless of solar geometry (ClockTime, Midnight).
#   Shape — how the curve travels; each shape's two Time args mean whatever
#           defines that shape (Sine: trough/peak; Decay/Linear: from/to).
#   Curve — a coverage-checked set of shapes plus the times fed from outside.
#
# Nothing reaches for solar geometry except through a Time's own `resolve_time`.

# ---------------------------------------------------------------------------
# Time
# ---------------------------------------------------------------------------

"""
    TimeOfDay

A within-day instant, resolved to a clock time (decimal hours) for a given day
by `resolve_time(::TimeOfDay, day)`. Solar events carry an optional offset in hours.
(Named `TimeOfDay`, not `Time`, to avoid colliding with `Dates.Time`.)
"""
abstract type TimeOfDay end

"""
    Sunrise(offset = 0)

Local sunrise, plus `offset` hours.
"""
struct Sunrise <: TimeOfDay
    offset::Float64
end
Sunrise() = Sunrise(0.0)

"""
    Sunset(offset = 0)

Local sunset, plus `offset` hours.
"""
struct Sunset <: TimeOfDay
    offset::Float64
end
Sunset() = Sunset(0.0)

"""
    Midday(offset = 0)

Solar noon (sun highest), plus `offset` hours.
"""
struct Midday <: TimeOfDay
    offset::Float64
end
Midday() = Midday(0.0)

"""
    Midnight(offset = 0)

Fixed clock midnight (00:00), plus `offset` hours — not tied to solar geometry
(unlike [`Sunrise`](@ref)/[`Sunset`](@ref)/[`Midday`](@ref)).
"""
struct Midnight <: TimeOfDay
    offset::Float64
end
Midnight() = Midnight(0.0)

"""
    ClockTime(hour)

A fixed local clock time, `hour` in decimal hours (`9.0`, `15.5`).
"""
struct ClockTime <: TimeOfDay
    hour::Float64
end

# Resolve to a clock time (decimal hours) for one day. `day` carries the day's
# solar geometry as decimal hours: `sunrise`, `sunset`, `solar_noon`.
@inline resolve_time(t::ClockTime, day) = t.hour
@inline resolve_time(t::Sunrise, day) = day.sunrise + t.offset
@inline resolve_time(t::Sunset, day) = day.sunset + t.offset
@inline resolve_time(t::Midday, day) = day.solar_noon + t.offset
@inline resolve_time(t::Midnight, day) = t.offset

# Nominal cyclic position (hours, 0–24) — for the *structural* coverage check
# only, which depends on the fixed event order, not the day's exact times.
@inline nominal_hour(t::Midnight) = mod(0.0 + t.offset, 24.0)
@inline nominal_hour(t::Sunrise) = mod(6.0 + t.offset, 24.0)
@inline nominal_hour(t::Midday) = mod(12.0 + t.offset, 24.0)
@inline nominal_hour(t::Sunset) = mod(18.0 + t.offset, 24.0)
@inline nominal_hour(t::ClockTime) = mod(t.hour, 24.0)

# ---------------------------------------------------------------------------
# Shape
# ---------------------------------------------------------------------------

"""
    Shape

How a curve travels between its defining times. Each shape's two `Time` fields
mean whatever defines it. A shape is either *bounded* — it covers exactly the
span between its two times — or *space-filling* — it covers whatever the bounded
shapes leave.
"""
abstract type Shape end

"""
    Sine(trough, peak)

A sinusoid pinned by its minimum at `trough` and maximum at `peak`. Space-filling
— the sinusoid is fully fixed by the two extrema and is evaluated across whatever
span the bounded shapes don't claim.
"""
struct Sine{A<:TimeOfDay, B<:TimeOfDay} <: Shape
    trough::A
    peak::B
end

"""
    Decay(from, to)

Monotonic relaxation from `from` to `to`. Bounded — covers exactly `[from, to)`.
"""
struct Decay{A<:TimeOfDay, B<:TimeOfDay} <: Shape
    from::A
    to::B
end

"""
    Linear(from, to)

Straight line from `from` to `to`. Bounded — covers exactly `[from, to)`.
"""
struct Linear{A<:TimeOfDay, B<:TimeOfDay} <: Shape
    from::A
    to::B
end

# Whether a shape fills leftover space (Sine) or claims a fixed span (Linear, Decay).
@inline is_space_filling(::Shape) = false
@inline is_space_filling(::Sine) = true

# The claimed span of a bounded shape, as a (from, to) Time pair.
@inline span(s::Linear) = (s.from, s.to)
@inline span(s::Decay) = (s.from, s.to)

# Both defining times of a shape (used to validate that inputs name real endpoints).
@inline endpoints(s::Sine) = (s.trough, s.peak)
@inline endpoints(s::Decay) = (s.from, s.to)
@inline endpoints(s::Linear) = (s.from, s.to)

# ---------------------------------------------------------------------------
# DielCurve
# ---------------------------------------------------------------------------

"""
    DielCurve(shapes; inputs)

A variable's within-day curve: a tuple of `shapes` (cyclic over the day) and
`inputs`, the times fed from outside. Holds no data — only structure.

The constructor checks, on the fixed solar-event order:

- the bounded shapes don't overlap, and with no space-filling shape they tile
  the whole day (a gap is an error);
- at most one space-filling shape (two leave coverage ambiguous);
- every `inputs` time is an endpoint of some shape.

Endpoints not listed in `inputs` are seams — their height is handed off from the
abutting shape at evaluation, not supplied.
"""
struct DielCurve{S<:Tuple, I<:Tuple}
    shapes::S
    inputs::I
    function DielCurve(shapes::Tuple, inputs::Tuple)
        _validate_coverage(shapes)
        _validate_inputs(shapes, inputs)
        return new{typeof(shapes), typeof(inputs)}(shapes, inputs)
    end
end
DielCurve(shapes::Tuple; inputs::Tuple) = DielCurve(shapes, inputs)
DielCurve(shapes...; inputs::Tuple) = DielCurve(shapes, inputs)

# Coverage check on nominal event positions: bounded spans must not overlap, and
# absent a space-filler must tile the 24 h cycle.
function _validate_coverage(shapes::Tuple)
    nfill = count(is_space_filling, shapes)
    nfill <= 1 || error("DielCurve: $nfill space-filling shapes; coverage is ambiguous (allow at most one)")

    bounded = filter(!is_space_filling, collect(shapes))
    # (start, length) per bounded span, lengths in (0, 24).
    arcs = map(bounded) do s
        f, t = span(s)
        a, b = nominal_hour(f), nominal_hour(t)
        len = mod(b - a, 24.0)
        len == 0.0 && error("DielCurve: shape $(s) spans zero length ($(f) to $(t))")
        (a, len)
    end
    isempty(arcs) && return nothing                # lone space-filler covers all
    sort!(arcs; by = first)
    n = length(arcs)
    for i in 1:n
        a, len = arcs[i]
        cover_end = a + len
        next_start = i < n ? arcs[i + 1][1] : arcs[1][1] + 24.0
        if cover_end > next_start + 1e-9
            error("DielCurve: shapes overlap near hour $(mod(next_start, 24.0))")
        elseif nfill == 0 && cover_end < next_start - 1e-9
            error("DielCurve: day not covered — gap from hour $(mod(cover_end, 24.0)) to $(mod(next_start, 24.0)); add a shape or a space-filling one")
        end
    end
    return nothing
end

# Every input must name an endpoint of some shape.
function _validate_inputs(shapes::Tuple, inputs::Tuple)
    eps = mapreduce(endpoints, (a, b) -> (a..., b...), shapes; init = ())
    for inp in inputs
        any(==(inp), eps) || error("DielCurve: input $(inp) is not an endpoint of any shape")
    end
    return nothing
end

# Endpoints not fed from outside — handed off from the abutting shape.
seams(c::DielCurve) =
    Tuple(e for s in c.shapes for e in endpoints(s) if !any(==(e), c.inputs))

# ---------------------------------------------------------------------------
# DielForcing: a curve plus its per-day input data
# ---------------------------------------------------------------------------

"""
    DielForcing(curve, values)

Pairs a declarative [`DielCurve`](@ref) with the per-day value series for its
input times. `values[i]` is the series for `curve.inputs[i]`. Data lives here,
never in the curve.
"""
struct DielForcing{C<:DielCurve, V<:Tuple}
    curve::C
    values::V
end

@inline n_days(f::DielForcing) = length(first(f.values))

"""
    ForcingSpec(curve, inputs::Tuple{Symbol, …})

The declarative half of a [`DielForcing`](@ref): a [`DielCurve`](@ref) plus the
*names* of the per-day quantities that feed the curve's input times, in
`curve.inputs` order. Holds no data. `bind_forcings` pairs it with value series to
produce a `DielForcing`.

A forcing *model* is a `NamedTuple` keyed by output quantity, each node a
`ForcingSpec` or a [`Derived`](@ref). It is the structure of a within-day model,
free of any per-slice data — see [`MINMAX_FORCING_MODEL`](@ref).
"""
struct ForcingSpec{C<:DielCurve, Inputs}
    curve::C
end
ForcingSpec(curve::DielCurve, inputs::Tuple) =
    ForcingSpec{typeof(curve), inputs}(curve)

@inline forcing_inputs(::ForcingSpec{C, Inputs}) where {C, Inputs} = Inputs

"""
    Derived(functor, inputs::Tuple{Symbol, …})

A quantity computed from other named quantities rather than interpolated from a
curve: `out[k] = functor(inputs[1][k], inputs[2][k], …)`, applied elementwise
over the hourly buffers. `inputs` are the source quantity names, in the functor's
argument order; they may be forcings or other derived quantities.

Relative humidity from vapour pressure and temperature is one instance —
`Derived(RelativeHumidityFromVapourPressureAndTemperature(method), (:actual_vapour_pressure, :reference_temperature))`
— but nothing here is humidity-specific.
"""
struct Derived{Names, F}
    functor::F
    Derived(functor::F, inputs::Tuple) where {F} = new{inputs, F}(functor)
end

@inline derived_inputs(::Derived{Names}) where {Names} = Names

"""
    RelativeHumidityFromVapourPressureAndTemperature(saturation_method)

Functor `(e, T) -> e / e_sat(T)`, clamped to `[0, 1]`. The vapour-pressure→humidity
relation, as a value to pass to [`Derived`](@ref).
"""
struct RelativeHumidityFromVapourPressureAndTemperature{M}
    saturation_method::M
end
@inline (rh::RelativeHumidityFromVapourPressureAndTemperature)(e, T) =
    clamp(NoUnits(e / vapour_pressure(rh.saturation_method, T)), 0.0, 1.0)

"""
    VapourPressureFromRelativeHumidityAndTemperature(saturation_method)

Functor `(rh, T) -> rh · e_sat(T)`, the inverse of
[`RelativeHumidityFromVapourPressureAndTemperature`](@ref). Turns a relative
humidity and temperature into an actual vapour pressure, as a value to pass to
[`Derived`](@ref).
"""
struct VapourPressureFromRelativeHumidityAndTemperature{M}
    saturation_method::M
end
@inline (f::VapourPressureFromRelativeHumidityAndTemperature)(rh, T) =
    rh * vapour_pressure(f.saturation_method, T)

# ---------------------------------------------------------------------------
# Evaluation — curve + day's solar geometry + input values → hourly profile
# ---------------------------------------------------------------------------

# The day's solar geometry in decimal hours, from the solar-radiation output.
# `hour_angle_sunrise` is an hour angle (hours from noon to sunrise), not a clock time.
@inline function solar_day(solar, iday)
    solar_noon = solar.hour_solar_noon[iday]
    sunrise = solar_noon - solar.hour_angle_sunrise[iday]
    sunset = 2 * solar_noon - sunrise
    return (; sunrise, sunset, solar_noon)
end

"""
    evaluate!(out, forcing, solar) -> out

Fill `out` (length `24 * n_days`) with the variable's hourly curve, evaluating
each shape over its covered span and handing seam heights off from the abutting
shape.
"""
function evaluate!(out, f::DielForcing, solar)
    # Sample in the output buffer's unit so the shape arithmetic stays in one
    # (non-affine) unit — mixing °C inputs with K differences is an affine error.
    for iday in 1:n_days(f)
        day = solar_day(solar, iday)
        vals = map(s -> convert(eltype(out), s[iday]), f.values)
        base = (iday - 1) * 24
        for h in 0:23
            out[base + h + 1] = sample(f.curve, Float64(h), day, vals)
        end
    end
    return out
end

# Curve value at clock hour `h`: the covering shape, evaluated.
sample(curve::DielCurve, h, day, vals) =
    sample_shape(covering_shape(curve, h, day), h, day, curve, vals)

# The shape covering hour `h`: a bounded shape whose span holds it, else the filler.
function covering_shape(curve::DielCurve, h, day)
    filler = nothing
    for s in curve.shapes
        if is_space_filling(s)
            filler = s
        else
            f, t = span(s)
            in_span(h, resolve_time(f, day), resolve_time(t, day)) && return s
        end
    end
    return filler
end

# Cyclic membership of `h` in [t1, t2).
@inline function in_span(h, t1, t2)
    t2c = t2 <= t1 ? t2 + 24 : t2
    hc = h < t1 ? h + 24 : h
    return t1 <= hc < t2c
end

# Evaluate one shape at `h`, resolving each endpoint's height (input or seam).
function sample_shape(shape, h, day, curve, vals)
    a, b = endpoints(shape)
    ta = resolve_time(a, day)
    tb = resolve_time(b, day)
    va = endpoint_height(a, ta, day, curve, vals)
    vb = endpoint_height(b, tb, day, curve, vals)
    return shape_value(shape, h, ta, va, tb, vb)
end

# Height at an endpoint: its input value, or — for a seam — handed off from the
# shape covering the instant just before it.
function endpoint_height(t, th, day, curve, vals)
    i = input_index(curve.inputs, t)
    i === nothing || return vals[i]
    prev = covering_shape(curve, mod(th - 1e-6, 24.0), day)
    return sample_shape(prev, th, day, curve, vals)
end

@inline function input_index(inputs::Tuple, t)
    for i in eachindex(inputs)
        inputs[i] == t && return i
    end
    return nothing
end

# Shape geometry between (t1, v1) and (t2, v2), evaluated at h.
@inline shape_value(::Linear, h, t1, v1, t2, v2) = v1 + (v2 - v1) * frac(h, t1, t2)
@inline shape_value(::Decay, h, t1, v1, t2, v2) = v2 + (v1 - v2) * exp(-3 * frac(h, t1, t2))
@inline shape_value(::Sine, h, t1, v1, t2, v2) = v1 + (v2 - v1) * (1 - cospi(frac(h, t1, t2))) / 2

# Fractional progress of `h` from t1 to t2, cyclic.
@inline function frac(h, t1, t2)
    t2c = t2 <= t1 ? t2 + 24 : t2
    hc = h < t1 ? h + 24 : h
    return (hc - t1) / (t2c - t1)
end

"""
    populate_quantities!(buffer_for, forcings, solar)

Fill every quantity's hourly buffer from `forcings` (a `NamedTuple` of
`DielForcing`/`Derived`). `buffer_for(::Val{name})` returns the buffer for a
quantity. Interpolated forcings are filled first, then derived quantities — so a
`Derived` reads its inputs after they exist. (A derived quantity that depends on
another derived quantity must be declared after it.)
"""
function populate_quantities!(buffer_for, forcings::NamedTuple{K}, solar) where {K}
    keyvals = map(Val, K)
    map(keyvals) do vk
        q = getproperty(forcings, _name(vk))
        q isa DielForcing && evaluate!(buffer_for(vk), q, solar)
        nothing
    end
    map(keyvals) do vk
        q = getproperty(forcings, _name(vk))
        q isa Derived && fill_derived!(buffer_for(vk), q, buffer_for)
        nothing
    end
    return nothing
end

@inline _name(::Val{k}) where {k} = k

# out[k] = functor(in₁[k], in₂[k], …) over each input quantity's buffer.
function fill_derived!(out, d::Derived, buffer_for)
    ins = map(buffer_for, map(Val, derived_inputs(d)))
    @inbounds for k in eachindex(out)
        out[k] = d.functor(map(b -> b[k], ins)...)
    end
    return out
end

# ---------------------------------------------------------------------------
# Binding a forcing model to value series
# ---------------------------------------------------------------------------

"""
    bind_forcings(model, valuesource) -> NamedTuple

Turn a declarative forcing `model` (a `NamedTuple` of [`ForcingSpec`](@ref) and
[`Derived`](@ref) nodes) into one of [`DielForcing`](@ref) / `Derived`, ready for
an environment. `valuesource(name)` returns the per-day value series for an input
quantity `name`. `ForcingSpec`s are bound to their input series; `Derived` nodes
pass through (they resolve against the populated buffers at evaluation).
"""
@inline bind_forcings(model::NamedTuple, valuesource) =
    map(node -> _bind_node(node, valuesource), model)

@inline _bind_node(spec::ForcingSpec, valuesource) =
    DielForcing(spec.curve, map(valuesource, forcing_inputs(spec)))
@inline _bind_node(d::Derived, valuesource) = d

# ---------------------------------------------------------------------------
# Standard min/max forcing model
# ---------------------------------------------------------------------------

"""
    minmax_forcing_model(; temperature_sunrise_offset=0, temperature_midday_offset=1,
                            wind_sunrise_offset=0, wind_midday_offset=1,
                            humidity_sunrise_offset=1, humidity_midday_offset=0,
                            cloud_sunrise_offset=1, cloud_midday_offset=0) -> NamedTuple

The classic daily min/max within-day structure, as a data-free forcing model
(`NamedTuple` of [`ForcingSpec`](@ref) keyed by output variable). Temperature is a
sine to an afternoon peak with an overnight decay; wind, humidity and cloud are
piecewise-linear between their daily extremes and the daily mean at true midnight
(humidity and cloud high near sunrise, low at midday; wind the reverse).

The `_offset` arguments (hours) shift each variable's dawn-side and midday-side
anchor away from `Sunrise()`/`Midday()`; the defaults match NicheMapR's standard
`TIMINS`/`TIMAXS`. Types don't depend on the offset values, so passing non-default
offsets carries no runtime cost or type instability — the curve is built once, here.

Bind value series to it with `bind_forcings`, or use the [`minmax_forcings`](@ref)
convenience, which accepts the same offset keywords.
"""
function minmax_forcing_model(;
    temperature_sunrise_offset = 0, temperature_midday_offset = 1,
    wind_sunrise_offset = 0, wind_midday_offset = 1,
    humidity_sunrise_offset = 1, humidity_midday_offset = 0,
    cloud_sunrise_offset = 1, cloud_midday_offset = 0,
)
    return (;
        reference_temperature = ForcingSpec(
            DielCurve((Sine(Sunrise(temperature_sunrise_offset), Midday(temperature_midday_offset)),
                       Decay(Sunset(), Sunrise(temperature_sunrise_offset))),
                inputs = (Sunrise(temperature_sunrise_offset), Midday(temperature_midday_offset))),
            (:reference_temperature_min, :reference_temperature_max)),
        reference_wind_speed = ForcingSpec(
            DielCurve((Linear(Midnight(), Sunrise(wind_sunrise_offset)),
                       Linear(Sunrise(wind_sunrise_offset), Midday(wind_midday_offset)),
                       Linear(Midday(wind_midday_offset), Midnight())),
                inputs = (Midnight(), Sunrise(wind_sunrise_offset), Midday(wind_midday_offset))),
            (:reference_wind_speed_mean, :reference_wind_speed_min, :reference_wind_speed_max)),
        reference_humidity = ForcingSpec(
            DielCurve((Linear(Midnight(), Sunrise(humidity_sunrise_offset)),
                       Linear(Sunrise(humidity_sunrise_offset), Midday(humidity_midday_offset)),
                       Linear(Midday(humidity_midday_offset), Midnight())),
                inputs = (Midnight(), Sunrise(humidity_sunrise_offset), Midday(humidity_midday_offset))),
            (:reference_humidity_mean, :reference_humidity_max, :reference_humidity_min)),
        cloud_cover = ForcingSpec(
            DielCurve((Linear(Midnight(), Sunrise(cloud_sunrise_offset)),
                       Linear(Sunrise(cloud_sunrise_offset), Midday(cloud_midday_offset)),
                       Linear(Midday(cloud_midday_offset), Midnight())),
                inputs = (Midnight(), Sunrise(cloud_sunrise_offset), Midday(cloud_midday_offset))),
            (:cloud_cover_mean, :cloud_cover_max, :cloud_cover_min)),
    )
end

"""
    MINMAX_FORCING_MODEL

[`minmax_forcing_model`](@ref) built with its default offsets.
"""
const MINMAX_FORCING_MODEL = minmax_forcing_model()

"""
    minmax_forcings(; reference_temperature_min, reference_temperature_max,
                      reference_wind_speed_min, reference_wind_speed_max,
                      reference_humidity_min, reference_humidity_max,
                      cloud_cover_min, cloud_cover_max, kwargs...) -> NamedTuple

Bind the standard min/max forcing model to per-day min/max value series,
returning a `NamedTuple` of [`DielForcing`](@ref) keyed by output variable.
`kwargs...` are the anchor-offset keywords of [`minmax_forcing_model`](@ref).
"""
function minmax_forcings(;
    reference_temperature_min, reference_temperature_max,
    reference_wind_speed_min, reference_wind_speed_max,
    reference_humidity_min, reference_humidity_max,
    cloud_cover_min, cloud_cover_max,
    kwargs...,
)
    # Lazy: callers bind forcings before buffers are populated, then mutate
    # min/max in place afterward. `.+`/`./` would freeze a mean of zeros now
    # instead of tracking the later mutation.
    _mean(a, b) = Base.Broadcast.broadcasted(/, Base.Broadcast.broadcasted(+, a, b), 2)
    reference_wind_speed_mean = _mean(reference_wind_speed_min, reference_wind_speed_max)
    reference_humidity_mean = _mean(reference_humidity_min, reference_humidity_max)
    cloud_cover_mean = _mean(cloud_cover_min, cloud_cover_max)
    series = (; reference_temperature_min, reference_temperature_max,
        reference_wind_speed_min, reference_wind_speed_max, reference_wind_speed_mean,
        reference_humidity_min, reference_humidity_max, reference_humidity_mean,
        cloud_cover_min, cloud_cover_max, cloud_cover_mean)
    model = isempty(kwargs) ? MINMAX_FORCING_MODEL : minmax_forcing_model(; kwargs...)
    return bind_forcings(model, name -> getproperty(series, name))
end
