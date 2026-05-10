"""
    NonConsecutiveDayMode(; ndmax=3)

Each output day j is integrated `ndmax` times back-to-back with its own forcings,
starting fresh from a uniform TINS reset. T continues across iterations *within*
a day; the first iteration always resets to `TINS(:,j) = mean(daily air temperature)`.
Output is written only on the last iter. Days are independent — no T continuity
across output days.

Replicates Fortran NicheMapR's `microdaily=0` flow. The two `IFINAL` counters in
the Fortran source are *separate local variables* (no COMMON block), not one shared
counter:
- `MICROCLIMATE.f:209` `DATA IFINAL/1/`, incremented at `:984`, clamped at `:988-989`
  (`IF IFINAL > ND THEN IFINAL = 1`). Drives the `T = TINS` reset at `:915` when
  `IFINAL == 1`. Cycles `1 → 2 → … → ND → 1`, so the reset fires every ND passes
  starting from pass 1.
- `OSUB.f:158` `DATA IFINAL/0/`, incremented once per SFODE call at `:316` (only at
  `TIME=0`). Output gated by `IFINAL ≥ ND` at `:353`; reset to 0 at `:372` (only
  reached on output passes — skipped passes return at `:353` before `:372`). So
  output fires every ND passes, on pass `ND, 2*ND, 3*ND, …`.

The two counters are aligned so each NDMAX-pass block has T reset on its first pass
and output written on its last. Effect for `ndmax=3`: every output day gets 3
integrations of itself, starting from TINS, with the third writing output. Day 1
has *no* special case — there's nothing in the source code that distinguishes it,
since the two IFINAL counters reach their fire-conditions on the same passes.
"""
struct NonConsecutiveDayMode <: AbstractTimeMode
    ndmax::Int
end
NonConsecutiveDayMode(; ndmax::Int=3) = NonConsecutiveDayMode(ndmax)

# `independent_days` is the historical predicate — true means *every* day starts
# fresh, false means days run continuously. NonConsecutiveDayMode is "true" — every
# output day resets T to its own TINS profile (matching Fortran's `MICROCLIMATE.f:915`
# firing on the first pass of every NDMAX-pass block).
independent_days(::NonConsecutiveDayMode) = true
is_reset_day(::NonConsecutiveDayMode, j::Int) = true
reset_phase_per_iter(::NonConsecutiveDayMode) = true
reset_moisture_per_day(::NonConsecutiveDayMode) = true
reset_snow_per_day(::NonConsecutiveDayMode) = true
iter_resets_T(::NonConsecutiveDayMode) = false

# Each output day gets `ndmax` SFODE passes — pass 1 starts from TINS, passes
# 2..ndmax continue with the same forcings. Output is written only on the last
# (`ndmax`th) pass. The convergence kwarg is ignored — NDMAX drives the count.
iterations_for_day(mode::NonConsecutiveDayMode, ::AbstractSoilTemperatureConvergence, day_index) =
    mode.ndmax
