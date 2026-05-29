# snow_obs_comparison.jl
#
# Compares the Julia microclimate model against NicheMapR predictions and SNOTEL 329
# field observations (snow depth, soil temperatures, soil moisture) for a 4-year run
# at SNOTEL site 329, Utah (39.14°N, 111.56°W, 2435 m), 2010-01-01 → 2013-12-31.
#
# Required data in test/data/scan329_long/:
#   - forcing.csv      (1461 days of daily Tmin/Tmax/RH/wind/rain/cloud)
#   - initial_sm.csv   (initial soil moisture, 10 coarse nodes)
#   - nmr_hourly.csv   (NicheMapR hourly output, unified format)
#   - obs.csv          (SNOTEL ground-truth observations)
#
# Run manually (not in CI):
#   julia test/snow_obs_comparison.jl

using Microclimate
using SolarRadiation
using FluidProperties
using Unitful
using CSV, DataFrames, Dates

testdir = @__DIR__
datadir = joinpath(testdir, "data", "scan329_long")

# ── Load NicheMapR hourly comparison data ────────────────────────────────────
nmr = DataFrame(CSV.File(joinpath(datadir, "nmr_hourly.csv")))

# ── Hard-wired site constants (SNOTEL 329, Utah) ─────────────────────────────
const ELEVATION   = 2435.35u"m"
const LATITUDE    = (39.0 + 8.2098 / 60.0)u"°"       # 39.1368°N
const LONGITUDE   = -(111.0 + 33.4878 / 60.0)u"°"    # 111.5580°W
const RUF         = 0.004u"m"
const USRHYT      = 0.01u"m"
const REFHYT      = 2.0u"m"
const ALBEDO      = 0.15
const EMISSIVITY  = 0.95
const START_DATE  = DateTime(2010, 1, 1)

# 10 coarse → 19 fine soil nodes (depths in cm: 0, 1.25, 2.5, 3.75, 5, 7.5, 10, 12.5, 15, 17.5,
# 20, 25, 30, 40, 50, 75, 100, 150, 200)
const DEPTHS = ([0.0, 1.25, 2.5, 3.75, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5,
                 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0] ./ 100.0) .* u"m"

# Soil properties (19 fine nodes)
const BULK_DENSITY        = 1.3u"Mg/m^3"
const SATURATION_MOISTURE = 0.4922u"m^3/m^3"
const MINERAL_DENSITY     = 2.56u"Mg/m^3"
const MINERAL_CONDUCTIVITY  = [0.2, 0.2, 0.2, 1.35, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5,
                                2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]u"W/m/K"
const MINERAL_HEAT_CAPACITY = [1920.0, 1920.0, 1920.0, 1395.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0,
                                870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0]u"J/kg/K"
const DE_VRIES_SHAPE_FACTOR = 0.1
const AIR_ENTRY_POTENTIAL   = 1.1u"J/kg"
const SAT_HYDRAULIC_COND    = 0.0037u"kg*s/m^3"
const CAMPBELL_B            = 4.5
const ROOT_DENSITY = [0.0, 0.0, 82000.0, 80000.0, 78000.0, 74000.0, 71000.0, 64000.0, 58000.0, 48000.0,
                      40000.0, 18000.0, 9000.0, 6000.0, 8000.0, 4000.0, 4000.0, 0.0, 0.0]u"m/m^3"

# Initial soil moisture: 10 coarse layers from initial_sm.csv → expand to 19 fine nodes
const _INITIAL_SM_COARSE = let df = DataFrame(CSV.File(joinpath(datadir, "initial_sm.csv")))
    # The file is keyed by `layer` (1..10) with `moisture` values
    [df[df.layer .== i, :moisture][1] for i in 1:10]
end
const INITIAL_SM = let coarse = _INITIAL_SM_COARSE, n = length(coarse)
    result = Vector{Float64}(undef, 2n - 1)
    for i in 1:n;     result[2i-1] = coarse[i]; end
    for i in 1:n-1;   result[2i]   = (coarse[i] + coarse[i+1]) / 2; end
    result
end
INITIAL_ST = fill(274.15u"K", length(DEPTHS))
INITIAL_ST[19] = (273.15 + 6.4)u"K"
INITIAL_ST[18] = (273.15 + (6.4 + 3.0) / 2.0)u"K"
INITIAL_ST[17] = (273.15 + 3)u"K"
INITIAL_ST[16] = (273.15 + 3)u"K"
INITIAL_ST[15] = (273.15 + 3)u"K"

# ── Load daily forcing ───────────────────────────────────────────────────────
forcing = DataFrame(CSV.File(joinpath(datadir, "forcing.csv")))
const NDAYS = nrow(forcing)
@info "Loaded forcing" NDAYS
days2do = 1:NDAYS

# ── Site ─────────────────────────────────────────────────────────────────────
site = Site(;
    latitude       = LATITUDE,
    longitude      = LONGITUDE,
    elevation      = ELEVATION,
    slope          = 0.0u"°",
    aspect         = 0.0u"°",
    horizon_angles = zeros(24)u"°",
    sky_view_fraction = 1.0,
    albedo         = ALBEDO,
    roughness_height = RUF,
    atmospheric_pressure = atmospheric_pressure(ELEVATION),
)
boundary_layer_model = MoninObukhov(; karman_constant=0.4, dyer_constant=16.0)

# ── Soil models ──────────────────────────────────────────────────────────────
soil_properties_model = CampbelldeVriesSoilProperties(;
    de_vries_shape_factor = DE_VRIES_SHAPE_FACTOR,
    mineral_conductivity  = MINERAL_CONDUCTIVITY,
    mineral_heat_capacity = MINERAL_HEAT_CAPACITY,
    recirculation_power   = 4.0,
    return_flow_threshold = 0.162,
)

soil_profile = SoilProfile(;
    bulk_density    = fill(BULK_DENSITY, length(DEPTHS)),
    mineral_density = fill(MINERAL_DENSITY, length(DEPTHS)),
)

soil_hydraulic_model = example_soil_hydraulic_model(DEPTHS;
    air_entry_water_potential        = fill(AIR_ENTRY_POTENTIAL, length(DEPTHS)),
    saturated_hydraulic_conductivity = fill(SAT_HYDRAULIC_COND, length(DEPTHS)),
    campbell_b_parameter             = fill(CAMPBELL_B, length(DEPTHS)),
    root_density                     = ROOT_DENSITY,
)

# ── Daily min/max environment ────────────────────────────────────────────────
environment_minmax = DailyMinMaxEnvironment(;
    reference_temperature_min = forcing.TMINN[days2do] .* u"°C",
    reference_temperature_max = forcing.TMAXX[days2do] .* u"°C",
    reference_wind_min        = forcing.WNMINN[days2do] .* u"m/s",
    reference_wind_max        = forcing.WNMAXX[days2do] .* u"m/s",
    reference_humidity_min    = forcing.RHMINN[days2do] ./ 100.0,
    reference_humidity_max    = forcing.RHMAXX[days2do] ./ 100.0,
    cloud_min                 = zeros(NDAYS),
    cloud_max                 = forcing.CCMAXX[days2do] ./ 100.0,
    minima_times              = [0.0, 0.0, 1.0, 1.0],
    maxima_times              = [1.0, 1.0, 0.0, 0.0],
)

# ── Daily timeseries — deep soil temp from NMR (D200cm at midnight each day) ─
environment_daily = DailyTimeseries(;
    shade                 = zeros(NDAYS),
    soil_wetness          = zeros(NDAYS),
    surface_emissivity    = fill(EMISSIVITY, NDAYS),
    cloud_emissivity      = fill(EMISSIVITY, NDAYS),
    rainfall              = forcing.RAINFALL[days2do] .* u"kg/m^2",
    deep_soil_temperature = nmr.D50cm[1:24:size(nmr, 1)][days2do] .* u"°C",
    leaf_area_index       = fill(0.1, NDAYS),
)

# ── Hourly placeholder ───────────────────────────────────────────────────────
environment_hourly = HourlyTimeseries(;
    pressure              = fill(atmospheric_pressure(ELEVATION), NDAYS * 24),
    reference_temperature = nothing,
    reference_humidity    = nothing,
    reference_wind_speed  = nothing,
    global_radiation      = nothing,
    cloud_cover           = nothing,
    rainfall              = nothing,
    zenith_angle          = nothing,
    longwave_radiation    = nothing,
)

radiation = RadiationModel(;
    solar_radiation_model = SolarProblem(; diffuse_model = SolarRadiation.NoScattering()),
)

snow_model = SnowModel(;
    snow_temperature_threshold = 1.5u"°C",
    snow_density = 0.375u"g/cm^3",
    snow_melt_factor = 1.0,
    undercatch = 1.0,
    rain_multiplier = 1.0,
    rain_melt_factor = 0.0125,
    density_function = (0.5979, 0.2178, 0.001, 0.0038),
    snow_conductivity = 0.0u"W/m/K",
    canopy_interception = 0.0,
)

# ── Build and solve ──────────────────────────────────────────────────────────
config = MicroConfig(;
    time_mode = ConsecutiveDayMode(; spinup_first_day=false),
    convergence = FixedSoilTemperatureIterations(3),
    rainfall_schedule = DailyRainfall(),
    soil_moisture_strategy = DynamicSoilMoisture(),
)

model = MicroModel(;
    hours   = collect(0.0:1:23.0),
    depths  = DEPTHS,
    heights = [USRHYT, REFHYT],
    radiation,
    soil_profile,
    soil_properties_model,
    soil_hydraulic_model,
    snow_model,
    boundary_layer_model,
    config,
)
inputs = MicroInputs(;
    site,
    environment_minmax,
    environment_daily,
    environment_hourly,
    initial_soil_temperature = INITIAL_ST,
    initial_soil_moisture = Vector{Float64}(INITIAL_SM),
)
problem = MicroProblem(model, inputs; days = forcing.DOY[days2do])

println("Running Julia microclimate model for $NDAYS days (2010–2013)...")
@time micro_out = Microclimate.solve(problem)
println("Simulation complete.")

# ── Load SNOTEL 329 observations ─────────────────────────────────────────────
# obs columns: DateTime, SNWD_mm (mm), WTEQ_in10 (inches × 10),
#              STO_5cm / STO_20cm / STO_50cm (°C), SMS_5cm / SMS_20cm / SMS_50cm (%)
const OBS_PATH = joinpath(datadir, "obs.csv")
obs = DataFrame(CSV.File(OBS_PATH; missingstring = ["NA", ""]))
# Rename WTEQ_in10 → WTEQ_in for forward compatibility (it's actually inches × 10)
if "WTEQ_in10" in String.(propertynames(obs)) && !("WTEQ_in" in String.(propertynames(obs)))
    rename!(obs, :WTEQ_in10 => :WTEQ_in)
end

# ── Extract Julia and NMR series at standard depths ──────────────────────────
# 19-node DEPTHS indices: 5=5cm, 7=10cm, 11=20cm, 15=50cm, 17=100cm
julia_D5cm  = ustrip.(u"°C", micro_out.soil_temperature[:, 5])
julia_D20cm = ustrip.(u"°C", micro_out.soil_temperature[:, 11])
julia_D50cm = ustrip.(u"°C", micro_out.soil_temperature[:, 15])
julia_WC5cm  = micro_out.soil_moisture[:, 5]
julia_WC20cm = micro_out.soil_moisture[:, 11]
julia_WC50cm = micro_out.soil_moisture[:, 15]
snow_depth_julia   = ustrip.(u"cm",     micro_out.snow_depth)
snow_density_julia = ustrip.(u"g/cm^3", micro_out.snow_density)

nhours_out = NDAYS * 24
nmr_sub = nmr[1:nhours_out, :]

# ── Hourly DateTime axis for mapping obs onto model hours ────────────────────
t_full = collect(range(START_DATE, step=Hour(1), length=nhours_out))

# Build a same-length series of obs values aligned to t_full by exact-hour match.
# `obs.DateTime` is at hourly resolution starting at `START_DATE`, so this is just
# a 1:1 mapping (no interpolation needed).
function obs_aligned(col::Symbol, t::Vector{DateTime}, obs::DataFrame)
    # Build a Dict (DateTime → value) and look up. Missing observations remain `missing`.
    d = Dict{DateTime, Union{Float64, Missing}}()
    for r in eachrow(obs)
        d[r.DateTime] = r[col]
    end
    result = Vector{Union{Float64, Missing}}(missing, length(t))
    for i in eachindex(t)
        if haskey(d, t[i]); result[i] = d[t[i]]; end
    end
    return result
end

# SNOTEL column-name lore: despite the "_mm" suffix, SNWD values are in INCHES
# (×2.54 → cm). Same for WTEQ_in10 — inches, not inches × 10. The misleading
# names are preserved for compatibility with the source data files.
obs_SNWD_cm   = obs_aligned(:SNWD_mm, t_full, obs)
obs_SNWD_cm  .= obs_SNWD_cm .* 2.54
obs_SWE_cm    = obs_aligned(:WTEQ_in, t_full, obs)
obs_SWE_cm   .= obs_SWE_cm .* 2.54
obs_STO_5cm   = obs_aligned(:STO_5cm, t_full, obs)
obs_STO_20cm  = obs_aligned(:STO_20cm, t_full, obs)
obs_STO_50cm  = obs_aligned(:STO_50cm, t_full, obs)

# Helper: skip-missing diff stats. Returns (mean, max_abs, n_valid) over rows where
# both series are present.
function diff_stats(a::AbstractVector, b::AbstractVector)
    n = min(length(a), length(b))
    s = 0.0; smax = 0.0; nv = 0
    for i in 1:n
        ai = a[i]; bi = b[i]
        ai === missing && continue
        bi === missing && continue
        d = ai - bi
        s += d
        if abs(d) > smax; smax = abs(d); end
        nv += 1
    end
    return nv == 0 ? (NaN, NaN, 0) : (s/nv, smax, nv)
end

# ── Headline comparison: Julia-vs-obs and NMR-vs-obs ─────────────────────────
println("\n══════════ Julia vs OBSERVATIONS vs NicheMapR (3-way diff stats) ══════════")
println(rpad("variable", 18), rpad("mean Δ(Jul−Obs)", 17), rpad("max|J−O|", 11),
        rpad("mean Δ(NMR−Obs)", 17), rpad("max|N−O|", 11), "n_obs")
for (name, jl, nmr_col, ob) in [
    ("snow depth (cm)",  snow_depth_julia, nmr_sub.SNOWDEP, obs_SNWD_cm),
    ("soil T 5cm (°C)",  julia_D5cm,       nmr_sub.D5cm,    obs_STO_5cm),
    ("soil T 20cm (°C)", julia_D20cm,      nmr_sub.D20cm,   obs_STO_20cm),
    ("soil T 50cm (°C)", julia_D50cm,      nmr_sub.D50cm,   obs_STO_50cm),
]
    mJ, xJ, nv = diff_stats(jl, ob)
    mN, xN, _  = diff_stats(nmr_col, ob)
    println(rpad(name, 18),
            rpad(round(mJ; digits=3), 17), rpad(round(xJ; digits=3), 11),
            rpad(round(mN; digits=3), 17), rpad(round(xN; digits=3), 11), nv)
end

# Snow water equivalent (SWE) — Julia and NMR compute SWE = depth × density
julia_SWE = snow_depth_julia .* snow_density_julia
nmr_SWE   = nmr_sub.SNOWDEP .* nmr_sub.SNOWDENS
mJ, xJ, nv = diff_stats(julia_SWE, obs_SWE_cm)
mN, xN, _  = diff_stats(nmr_SWE,   obs_SWE_cm)
println(rpad("SWE (cm w.e.)", 18),
        rpad(round(mJ; digits=3), 17), rpad(round(xJ; digits=3), 11),
        rpad(round(mN; digits=3), 17), rpad(round(xN; digits=3), 11), nv)

# ── Per-year breakdown (which model wins each year?) ─────────────────────────
println("\n── Per-year max-abs deviation vs observations ──")
println(rpad("year", 6), rpad("var", 18),
        rpad("|J−O|max", 11), rpad("|N−O|max", 11), "winner (closer to obs)")
for yr in 2010:2013
    yr_idx = findall(i -> year(t_full[i]) == yr, eachindex(t_full))
    for (name, jl, nmr_col, ob) in [
        ("snow depth",  snow_depth_julia, nmr_sub.SNOWDEP, obs_SNWD_cm),
        ("soil T 5cm",  julia_D5cm,       nmr_sub.D5cm,    obs_STO_5cm),
        ("soil T 20cm", julia_D20cm,      nmr_sub.D20cm,   obs_STO_20cm),
        ("soil T 50cm", julia_D50cm,      nmr_sub.D50cm,   obs_STO_50cm),
    ]
        _, xJ, _ = diff_stats(jl[yr_idx],      ob[yr_idx])
        _, xN, _ = diff_stats(nmr_col[yr_idx], ob[yr_idx])
        winner = isnan(xJ) || isnan(xN) ? "—" : (xJ < xN ? "Julia" : "NicheMapR")
        println(rpad(yr, 6), rpad(name, 18),
                rpad(round(xJ; digits=3), 11), rpad(round(xN; digits=3), 11), winner)
    end
end

# ── Spring melt analysis (the divergence regime) ─────────────────────────────
# For each year, find first hour where obs snow depth dropped to < 1 cm after peak
# and compare what Julia / NMR show on that date.
println("\n── Spring melt-out date (first day obs snow < 1 cm after peak) ──")
for yr in 2010:2013
    yr_idx = findall(i -> year(t_full[i]) == yr, eachindex(t_full))
    obs_yr = obs_SNWD_cm[yr_idx]
    has_data = collect(skipmissing(obs_yr))
    if isempty(has_data)
        println("  $yr: no obs")
        continue
    end
    peak_val, peak_off = findmax(replace(obs_yr, missing => -Inf))
    melt_out = nothing
    for off in (peak_off + 1):length(obs_yr)
        v = obs_yr[off]
        v === missing && continue
        if v < 1.0
            melt_out = off
            break
        end
    end
    if isnothing(melt_out)
        println("  $yr: peak=$(round(peak_val; digits=1)) cm, never melts below 1 cm in record")
        continue
    end
    melt_step = yr_idx[melt_out]
    println("  $yr: obs peak=$(round(peak_val; digits=1)) cm on $(t_full[yr_idx[peak_off]]), ",
            "melt-out on $(t_full[melt_step]) — Julia=$(round(snow_depth_julia[melt_step]; digits=1)) cm, ",
            "NMR=$(round(nmr_sub.SNOWDEP[melt_step]; digits=1)) cm")
end

# ── Visual comparisons (run with Plots available) ────────────────────────────
plot_year  = nothing    # set e.g. 2013 to zoom
plot_month = nothing    # set e.g. 3 to zoom

const PLOTS_AVAILABLE = try
    @eval using Plots
    true
catch
    @info "Plots not available — skipping visualisations"
    false
end

PLOTS_AVAILABLE && let
    # Hourly index mask based on plot_year/plot_month
    if !isnothing(plot_year) && !isnothing(plot_month)
        title_suffix = " — $(Dates.monthname(plot_month)) $plot_year"
        hmask = findall(i -> year(t_full[i]) == plot_year && month(t_full[i]) == plot_month,
            eachindex(t_full))
    elseif !isnothing(plot_year)
        title_suffix = " — $plot_year"
        hmask = findall(i -> year(t_full[i]) == plot_year, eachindex(t_full))
    elseif !isnothing(plot_month)
        title_suffix = " — $(Dates.monthname(plot_month)) (all years)"
        hmask = findall(i -> month(t_full[i]) == plot_month, eachindex(t_full))
    else
        title_suffix = " — 2010–2013"
        hmask = eachindex(t_full)
    end
    t = t_full[hmask]

    # ── Snow panel ───────────────────────────────────────────────────────────
    pa1 = plot(t, snow_depth_julia[hmask], label="Julia", color=:black,
               title="Snow depth (cm)", ylabel="cm")
    plot!(pa1, t, nmr_sub.SNOWDEP[hmask], label="NicheMapR", color=:blue)
    plot!(pa1, t, obs_SNWD_cm[hmask], label="SNOTEL obs", color=:red, lw=1)

    pa2 = plot(t, julia_SWE[hmask], label="Julia", color=:black,
               title="Snow water equivalent (cm w.e.)", ylabel="cm")
    plot!(pa2, t, nmr_SWE[hmask], label="NicheMapR", color=:blue)
    plot!(pa2, t, obs_SWE_cm[hmask], label="SNOTEL obs", color=:red, lw=1)

    pa3 = plot(t, snow_density_julia[hmask], label="Julia", color=:black,
               title="Snow density (g/cm³)", ylabel="g/cm³")
    plot!(pa3, t, nmr_sub.SNOWDENS[hmask], label="NicheMapR", color=:blue)

    panel_snow = plot(pa1, pa2, pa3; layout=(3, 1), size=(1000, 800),
                     plot_title="Snow$title_suffix", link=:x)
    display(panel_snow)

    # ── Soil temperature panel ──────────────────────────────────────────────
    pb1 = plot(t, julia_D5cm[hmask], label="Julia", color=:black,
               title="Soil temp 5 cm (°C)", ylabel="°C")
    plot!(pb1, t, nmr_sub.D5cm[hmask], label="NicheMapR", color=:blue)
    plot!(pb1, t, obs_STO_5cm[hmask],  label="SNOTEL obs", color=:red, lw=1)

    pb3 = plot(t, julia_D20cm[hmask], label="Julia", color=:black,
               title="Soil temp 20 cm (°C)", ylabel="°C")
    plot!(pb3, t, nmr_sub.D20cm[hmask], label="NicheMapR", color=:blue)
    plot!(pb3, t, obs_STO_20cm[hmask],  label="SNOTEL obs", color=:red, lw=1)

    pb4 = plot(t, julia_D50cm[hmask], label="Julia", color=:black,
               title="Soil temp 50 cm (°C)", ylabel="°C")
    plot!(pb4, t, nmr_sub.D50cm[hmask], label="NicheMapR", color=:blue)
    plot!(pb4, t, obs_STO_50cm[hmask],  label="SNOTEL obs", color=:red, lw=1)

    panel_temp = plot(pb1, pb3, pb4; layout=(3, 1), size=(1000, 800),
                     plot_title="Soil temperature$title_suffix", link=:x)
    display(panel_temp)
end
