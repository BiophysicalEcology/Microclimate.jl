# snow_obs_comparison.jl
#
# Compares the Julia microclimate model against NicheMapR predictions and SNOTEL 329
# field observations (snow depth, soil temperatures, soil moisture) for a 4-year run
# at SNOTEL site 329, Utah (39.14°N, 111.56°W, 2435 m), 2010-01-01 to 2013-12-31.
#
# Run manually (not in CI — requires test/data/scan329/ data files):
#   julia test/snow_obs_comparison.jl
#
# Generate test/data/scan329/ first by running:
#   julia test/prepare_scan329_data.jl

using Microclimate
using SolarRadiation
using FluidProperties
using Unitful
using CSV, DataFrames, Dates
using Test

testdir = @__DIR__
datadir = joinpath(testdir, "data", "scan329")

# ── Hard-wired site constants (SNOTEL 329, Utah) ─────────────────────────────
# Derived from micro csv input/microinput.csv for the 2010–2013 run
const ELEVATION   = 2435.35u"m"
const LATITUDE    = (39.0 + 8.2098 / 60.0)u"°"       # 39.1368°N
const LONGITUDE   = -(111.0 + 33.4878 / 60.0)u"°"    # 111.5580°W
const RUF         = 0.004u"m"
const USRHYT      = 0.01u"m"
const REFHYT      = 2.0u"m"
const ALBEDO      = 0.25                               # REFLS (constant all days)
const EMISSIVITY  = 0.95                               # SLES (constant)
const TANNUL      = 6.43                               # deep soil temperature (°C, constant)
const NDAYS       = 1461

# Soil depths: [0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200] cm → metres
const DEPTHS = ([0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0] ./ 100.0) .* u"m"

# Soil thermal properties (from micro csv input/soilprop.csv, first layer used uniformly)
const BULK_DENSITY        = 1.3u"Mg/m^3"
const SATURATION_MOISTURE = 0.4922u"m^3/m^3"
const MINERAL_CONDUCTIVITY = 0.2u"W/m/K"
const MINERAL_HEAT_CAPACITY = 1920.0u"J/kg/K"
const MINERAL_DENSITY     = 2.56u"Mg/m^3"

# Initial soil moisture (m³/m³) extracted from micro csv input/moists.csv column 2
# (update by running prepare_scan329_data.jl if needed)
const INITIAL_SM = DataFrame(CSV.File(joinpath(datadir, "initial_sm.csv"))).moisture

# ── Load daily forcing ────────────────────────────────────────────────────────
forcing = DataFrame(CSV.File(joinpath(datadir, "forcing.csv")))
days2do = 1:NDAYS

# ── Terrain ───────────────────────────────────────────────────────────────────
micro_terrain = MicroTerrain(;
    elevation      = ELEVATION,
    roughness_height = RUF,
    karman_constant  = 0.4,
    dyer_constant    = 16.0,
    viewfactor       = 1.0,
)

solar_terrain = SolarTerrain(;
    slope              = 0.0u"°",
    aspect             = 0.0u"°",
    elevation          = ELEVATION,
    horizon_angles     = zeros(24)u"°",
    albedo             = ALBEDO,
    atmospheric_pressure = atmospheric_pressure(ELEVATION),
    latitude           = LATITUDE,
    longitude          = LONGITUDE,
)

# ── Soil models ───────────────────────────────────────────────────────────────
soil_thermal_model = CampbelldeVriesSoilThermal(;
    bulk_density         = BULK_DENSITY,
    mineral_density      = MINERAL_DENSITY,
    de_vries_shape_factor = 0.1,
    mineral_conductivity = MINERAL_CONDUCTIVITY,
    mineral_heat_capacity = MINERAL_HEAT_CAPACITY,
    saturation_moisture  = SATURATION_MOISTURE,
    recirculation_power  = 4.0,
    return_flow_threshold = 0.162,
)

soil_moisture_model = example_soil_hydraulics(DEPTHS;
    bulk_density    = BULK_DENSITY,
    mineral_density = MINERAL_DENSITY,
    root_density    = fill(0.0, length(DEPTHS))u"m/m^3",
    mode            = DynamicSoilMoisture(),
)

# ── Daily min/max environment (1461 days) ─────────────────────────────────────
# RHMINN/RHMAXX and CCMAXX are in %; divide by 100 for fractional
environment_minmax = DailyMinMaxEnvironment(;
    reference_temperature_min = forcing.TMINN[days2do] .* u"°C",
    reference_temperature_max = forcing.TMAXX[days2do] .* u"°C",
    reference_wind_min        = forcing.WNMINN[days2do] .* u"m/s",
    reference_wind_max        = forcing.WNMAXX[days2do] .* u"m/s",
    reference_humidity_min    = forcing.RHMINN[days2do] ./ 100.0,
    reference_humidity_max    = forcing.RHMAXX[days2do] ./ 100.0,
    cloud_min                 = zeros(NDAYS),
    cloud_max                 = forcing.CCMAXX[days2do] ./ 100.0,
    minima_times              = [1.0, 1.0, 1.0, 0.0],  # TIMINS1-4
    maxima_times              = [1.0, 1.0, 0.0, 0.0],  # TIMAXS1-4
)

# ── Daily timeseries ──────────────────────────────────────────────────────────
environment_daily = DailyTimeseries(;
    shade                 = zeros(NDAYS),
    soil_wetness          = zeros(NDAYS),
    surface_emissivity    = fill(EMISSIVITY, NDAYS),
    cloud_emissivity      = fill(EMISSIVITY, NDAYS),
    rainfall              = forcing.RAINFALL[days2do] .* u"kg/m^2",
    deep_soil_temperature = fill(TANNUL, NDAYS) .* u"°C",
    leaf_area_index       = fill(0.1, NDAYS),
)

# ── Hourly placeholder (pressure only) ───────────────────────────────────────
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

solar_model = SolarProblem(; scattered_uv = false)

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
    grass_shade = false,
)

# ── Build and solve ───────────────────────────────────────────────────────────
problem = MicroProblem(;
    latitude              = LATITUDE,
    days                  = forcing.DOY[days2do],
    hours                 = collect(0.0:1:23.0),
    depths                = DEPTHS,
    heights               = [USRHYT, REFHYT],
    solar_model,
    solar_terrain,
    micro_terrain,
    soil_moisture_model,
    soil_thermal_model,
    environment_minmax,
    environment_daily,
    environment_hourly,
    time_mode             = ConsecutiveDayMode(; spinup_first_day=true),
    convergence           = FixedSoilTemperatureIterations(3),
    hourly_rainfall       = false,
    snow_model,
    initial_soil_temperature = nothing,
    initial_soil_moisture = Vector{Float64}(INITIAL_SM),
)

println("Running Julia microclimate model for $NDAYS days (2010–2013)...")
@time micro_out = Microclimate.solve(problem)
println("Simulation complete.")

# ── Load NicheMapR comparison data ────────────────────────────────────────────
nmr = DataFrame(CSV.File(joinpath(datadir, "nmr_hourly.csv")))

# ── Load SNOTEL 329 observations ──────────────────────────────────────────────
obs = DataFrame(CSV.File(joinpath(datadir, "obs.csv"),
    missingstring = ["NA", ""],
    types = Dict(:SNWD_mm => Float64)
))
# DateTime column is auto-parsed by CSV
# Unit notes: SNWD_mm = snow depth in mm, WTEQ_in10 = SWE in tenths-of-inches
#   Convert SWE to cm: WTEQ_in10 × 2.54 / 10
#   Soil temps STO_*cm in °C, soil moisture SMS_*cm in %

# ── Comparisons ───────────────────────────────────────────────────────────────
# Extract Julia soil temperatures at key depths
# DEPTHS index: 3=5cm, 4=10cm, 6=20cm, 8=50cm  (1-indexed)
julia_D5cm  = ustrip.(u"°C", micro_out.soil_temperature[:, 3])
julia_D20cm = ustrip.(u"°C", micro_out.soil_temperature[:, 6])
julia_D50cm = ustrip.(u"°C", micro_out.soil_temperature[:, 8])

julia_WC5cm  = micro_out.soil_moisture[:, 3]
julia_WC20cm = micro_out.soil_moisture[:, 6]
julia_WC50cm = micro_out.soil_moisture[:, 8]

nhours_out = NDAYS * 24
nmr_sub = nmr[1:nhours_out, :]

# ── Detailed node diagnostics at a snow-covered hour ─────────────────────────
# Pick hour 12 of day 60 (deep snow period, week 9)
diag_step = (60 - 1) * 24 + 12
println("\n── Node temperatures at step $diag_step (day 60, hour 12) ──")
println("Snow depth: ", ustrip(u"cm", micro_out.snow_depth[diag_step]), " cm")
println("Snow density: ", ustrip(u"g/cm^3", micro_out.snow_density[diag_step]), " g/cm³")
println("Soil node temperatures (all depths):")
for (k, d) in enumerate(DEPTHS)
    st = ustrip(u"°C", micro_out.soil_temperature[diag_step, k])
    tc = ustrip(u"W/m/K", micro_out.soil_thermal_conductivity[diag_step, k])
    hc = ustrip(u"J/kg/K", micro_out.soil_heat_capacity[diag_step, k])
    bd = ustrip(u"kg/m^3", micro_out.soil_bulk_density[diag_step, k])
    sm = micro_out.soil_moisture[diag_step, k]
    println("  depth=$(ustrip(u"cm",d))cm: T=$(round(st,digits=3))°C  k=$(round(tc,digits=4))  cp=$(round(hc,digits=1))  rho=$(round(bd,digits=1))  SM=$(round(sm,digits=4))")
end
println("NicheMapR at same step: D5cm=$(nmr.D5cm[diag_step])  D20cm=$(nmr.D20cm[diag_step])  D50cm=$(nmr.D50cm[diag_step])  snow=$(nmr.SNOWDEP[diag_step])cm")

# ── Diagnostic summary ──────────────────────────────────────────────────────
snow_depth_julia   = ustrip.(u"cm",     micro_out.snow_depth)
snow_density_julia = ustrip.(u"g/cm^3", micro_out.snow_density)

println("\n── Diagnostic summary (first $NDAYS days) ──")
println("Soil temp 5cm  — Julia range: ", extrema(julia_D5cm), "  NicheMapR range: ", extrema(nmr_sub.D5cm))
println("Soil temp 20cm — Julia range: ", extrema(julia_D20cm), "  NicheMapR range: ", extrema(nmr_sub.D20cm))
println("Soil temp 50cm — Julia range: ", extrema(julia_D50cm), "  NicheMapR range: ", extrema(nmr_sub.D50cm))
println("Snow depth     — Julia range: ", extrema(snow_depth_julia), "  NicheMapR range: ", extrema(nmr_sub.SNOWDEP))
println("Snow density   — Julia range: ", extrema(snow_density_julia), "  NicheMapR range: ", extrema(nmr_sub.SNOWDENS))

# Check where biggest divergences occur
diff_5cm = julia_D5cm .- nmr_sub.D5cm
diff_20cm = julia_D20cm .- nmr_sub.D20cm
diff_50cm = julia_D50cm .- nmr_sub.D50cm
diff_snow = snow_depth_julia .- nmr_sub.SNOWDEP

println("\nSoil temp 5cm  diff — mean: ", round(sum(diff_5cm)/length(diff_5cm), digits=3), "  max abs: ", round(maximum(abs.(diff_5cm)), digits=3))
println("Soil temp 20cm diff — mean: ", round(sum(diff_20cm)/length(diff_20cm), digits=3), "  max abs: ", round(maximum(abs.(diff_20cm)), digits=3))
println("Soil temp 50cm diff — mean: ", round(sum(diff_50cm)/length(diff_50cm), digits=3), "  max abs: ", round(maximum(abs.(diff_50cm)), digits=3))
println("Snow depth diff    — mean: ", round(sum(diff_snow)/length(diff_snow), digits=3), "  max abs: ", round(maximum(abs.(diff_snow)), digits=3))

# Show a few time slices where divergence is large
worst_idx_5cm = argmax(abs.(diff_5cm))
worst_day_5cm = div(worst_idx_5cm - 1, 24) + 1
worst_hr_5cm = mod(worst_idx_5cm - 1, 24)
println("\nWorst 5cm divergence at step $worst_idx_5cm (day $worst_day_5cm, hour $worst_hr_5cm):")
println("  Julia: $(julia_D5cm[worst_idx_5cm])°C  NicheMapR: $(nmr_sub.D5cm[worst_idx_5cm])°C  Snow: $(snow_depth_julia[worst_idx_5cm])cm vs $(nmr_sub.SNOWDEP[worst_idx_5cm])cm")

worst_idx_snow = argmax(abs.(diff_snow))
worst_day_snow = div(worst_idx_snow - 1, 24) + 1
worst_hr_snow = mod(worst_idx_snow - 1, 24)
println("\nWorst snow depth divergence at step $worst_idx_snow (day $worst_day_snow, hour $worst_hr_snow):")
println("  Julia: $(snow_depth_julia[worst_idx_snow])cm  NicheMapR: $(nmr_sub.SNOWDEP[worst_idx_snow])cm")

# Show weekly averages to see trends
println("\n── Weekly soil temp 5cm comparison ──")
for w in 1:min(26, div(NDAYS, 7))
    r = ((w-1)*7*24+1):min(w*7*24, nhours_out)
    j_avg = round(sum(julia_D5cm[r])/length(r), digits=2)
    n_avg = round(sum(nmr_sub.D5cm[r])/length(r), digits=2)
    s_j = round(sum(snow_depth_julia[r])/length(r), digits=2)
    s_n = round(sum(nmr_sub.SNOWDEP[r])/length(r), digits=2)
    println("  Week $w: Julia=$(j_avg)C  NMR=$(n_avg)C  (snow: J=$s_j  N=$s_n cm)")
end

# ── Daily melt summary during divergence period ──
println("\n── Daily snow depth and melt summary (days 78-112) ──")
println("  day | J_snow  N_snow | J_daily_Δ  N_daily_Δ | J_T0surf  N_T5cm")
for day in 78:112
    step_start = (day - 1) * 24 + 1
    step_end = day * 24
    (step_start > nhours_out || step_end > nhours_out) && continue
    j_snow_end = snow_depth_julia[step_end]
    n_snow_end = nmr_sub.SNOWDEP[step_end]
    j_snow_prev = day > 1 ? snow_depth_julia[(day-2)*24 + 24] : 0.0
    n_snow_prev = day > 1 ? nmr_sub.SNOWDEP[(day-2)*24 + 24] : 0.0
    j_daily = j_snow_end - j_snow_prev
    n_daily = n_snow_end - n_snow_prev
    # Soil surface temp at noon
    noon_step = (day - 1) * 24 + 13
    j_tsurf = ustrip(u"°C", micro_out.soil_temperature[noon_step, 1])
    n_t5 = nmr_sub.D5cm[noon_step]
    println("  $(lpad(day,3)) | $(lpad(round(j_snow_end,digits=1),6))  $(lpad(round(n_snow_end,digits=1),6)) | $(lpad(round(j_daily,digits=1),9))  $(lpad(round(n_daily,digits=1),9)) | $(lpad(round(j_tsurf,digits=2),7))  $(lpad(round(n_t5,digits=1),5))")
end

# ── Soil moisture and thermal property comparison ──
println("\n── Soil moisture & conductivity comparison ──")
println("  day | SM5_J SM5_N | SM20_J SM20_N | SM50_J SM50_N | k5_J   k20_J   k50_J")
for day in [90, 180, 270, 360, 450, 540, 730, 1000]
    step = (day-1)*24 + 13
    (step > nhours_out || step > size(nmr_sub, 1)) && continue
    j5 = round(micro_out.soil_moisture[step, 3], digits=3)
    j20 = round(micro_out.soil_moisture[step, 6], digits=3)
    j50 = round(micro_out.soil_moisture[step, 8], digits=3)
    n5 = round(nmr_sub.WC5cm[step], digits=3)
    n20 = round(nmr_sub.WC20cm[step], digits=3)
    n50 = round(nmr_sub.WC50cm[step], digits=3)
    k5 = round(ustrip(u"W/m/K", micro_out.soil_thermal_conductivity[step, 3]), digits=3)
    k20 = round(ustrip(u"W/m/K", micro_out.soil_thermal_conductivity[step, 6]), digits=3)
    k50 = round(ustrip(u"W/m/K", micro_out.soil_thermal_conductivity[step, 8]), digits=3)
    println("  $(lpad(day,4)) | $j5  $n5 | $j20  $n20 | $j50  $n50 | $k5  $k20  $k50")
end

@testset "SNOTEL 329 — Julia vs NicheMapR ($(NDAYS)d)" begin
    @test size(micro_out.soil_temperature, 1) == nhours_out
    @test size(micro_out.soil_temperature, 2) == length(DEPTHS)

    @test julia_D5cm  ≈ nmr_sub.D5cm  rtol=0.5
    @test julia_D20cm ≈ nmr_sub.D20cm rtol=0.2
    @test julia_D50cm ≈ nmr_sub.D50cm rtol=0.1

    @test snow_depth_julia   ≈ nmr_sub.SNOWDEP  rtol=0.2
    @test snow_density_julia ≈ nmr_sub.SNOWDENS rtol=0.2
end

#── Visual comparisons — run manually (not in CI) ─────────────────────────────
# using Plots, Dates
# let
#     t = range(DateTime(2010,1,1), step=Hour(1), length=NDAYS*24)

#     # Snow depth: Julia vs NicheMapR vs observations
#     p1 = plot(t, nmr.SNOWDEP, label="NicheMapR", title="Snow depth (cm)", ylabel="cm")
#     plot!(p1, t, ustrip.(u"cm", micro_out.snow_depth), label="Julia", color=:red)
#     plot!(p1, obs.DateTime, obs.SNWD_mm .* 2.54, label="SNOTEL 329 obs", ms=1, color=:red)
#     display(p1)

#     # Soil temperature at 5 cm
#     p2 = plot(t, nmr.D5cm, label="NicheMapR", title="Soil temp 5 cm (°C)", ylabel="°C")
#     plot!(p2, t, julia_D5cm, label="Julia", color=:red)
#     plot!(p2, obs.DateTime, obs.STO_5cm, label="SNOTEL 329 obs", ms=1, color=:red)
#     display(p2)

#     # Soil temperature at 20 cm
#     p3 = plot(t, nmr.D20cm, label="NicheMapR", title="Soil temp 20 cm (°C)", ylabel="°C")
#     plot!(p3, t, julia_D20cm, label="Julia", color=:red)
#     plot!(p3, obs.DateTime, obs.STO_20cm, label="SNOTEL 329 obs", ms=1, color=:red)
#     display(p3)

#     # Soil moisture at 5 cm
#     p4 = plot(t, nmr.WC5cm, label="NicheMapR", title="Soil moisture 5 cm (m³/m³)", ylabel="m³/m³")
#     plot!(p4, t, julia_WC5cm, label="Julia", color=:red)
#     plot!(p4, obs.DateTime, obs.SMS_5cm ./ 100.0, label="SNOTEL 329 obs", ms=1, color=:red)
#     display(p4)
# end
