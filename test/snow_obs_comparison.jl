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

soil_moisture_model = example_soil_moisture_model(DEPTHS;
    bulk_density    = BULK_DENSITY,
    mineral_density = MINERAL_DENSITY,
    root_density    = fill(0.0, length(DEPTHS))u"m/m^3",
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
    rainfall              = (forcing.RAINFALL[days2do] ./ 1000.0) .* u"kg/m^2",
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
    iterate_day           = 3,
    daily                 = true,
    runmoist              = true,
    hourly_rainfall       = false,
    spinup                = true,
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
obs.DateTime = DateTime.(obs.DateTime, dateformat"yyyy-mm-dd HH:MM:SS")
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

@assert nrow(nmr) == NDAYS * 24 "Expected $(NDAYS*24) NicheMapR rows, got $(nrow(nmr))"

@testset "SNOTEL 329 — Julia vs NicheMapR" begin
    # Sanity checks: output dimensions correct
    @test size(micro_out.soil_temperature, 1) == NDAYS * 24
    @test size(micro_out.soil_temperature, 2) == length(DEPTHS)

    # TODO: uncomment soil temperature tests once snow model is fully implemented.
    # Snow insulates soil — without snow model, large discrepancies expected in winter.
    # @test julia_D5cm  ≈ nmr.D5cm  rtol=0.5
    # @test julia_D20cm ≈ nmr.D20cm rtol=0.2
    # @test julia_D50cm ≈ nmr.D50cm rtol=0.1

    # TODO: uncomment once snow model is implemented
    # snow_depth_julia   = ustrip.(u"cm",     micro_out.snow_depth)
    # snow_density_julia = ustrip.(u"g/cm^3", micro_out.snow_density)
    # @test snow_depth_julia   ≈ nmr.SNOWDEP  rtol=0.2
    # @test snow_density_julia ≈ nmr.SNOWDENS rtol=0.2
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
