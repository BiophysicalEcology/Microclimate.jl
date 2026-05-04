# prepare_scan329_data.jl
#
# Generates the NicheMapR comparison data for the SNOTEL 329 snow test.
# Run manually from the repo root:
#   julia test/prepare_scan329_data.jl
#
# Reads forcing, initial soil moisture, and observations from
# test/data/scan329/ (all self-contained in the package — no Dropbox needed),
# writes NicheMapR input CSVs, calls test/R/run_nmr.R via Rscript, then saves
# the NicheMapR outputs as nmr_hourly.csv for use by snow_obs_comparison.jl.

using CSV, DataFrames, Dates

const NDAYS = 365   # 2013 only

# Soil properties (19 fine nodes) — must match snow_obs_comparison.jl
const _MINERAL_CONDUCTIVITY  = [0.2, 0.2, 0.2, 1.35, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5,
                                 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
const _MINERAL_HEAT_CAPACITY = [1920.0, 1920.0, 1920.0, 1395.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0,
                                 870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0, 870.0]
const _ROOT_DENSITY = [0.0, 0.0, 82000.0, 80000.0, 78000.0, 74000.0, 71000.0, 64000.0, 58000.0, 48000.0,
                       40000.0, 18000.0, 9000.0, 6000.0, 8000.0, 4000.0, 4000.0, 0.0, 0.0]
# Indices of the 10 NMR coarse nodes in the 19-node fine array
const _NMR_IDX = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
# Initial soil temperature at each of the 10 NMR nodes (°C) — matches INITIAL_ST in snow_obs_comparison.jl
const _INIT_ST = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0, 3.0, 6.4]

testdir = realpath(joinpath(@__DIR__))
outdir  = joinpath(testdir, "data", "scan329")
println("Working in: $outdir")

# ── 1. Load daily forcing (already trimmed to 365 days) ───────────────────────
println("\nReading forcing.csv...")
forcing = DataFrame(CSV.File(joinpath(outdir, "forcing.csv")))
@assert nrow(forcing) == NDAYS "forcing.csv must have $NDAYS rows, got $(nrow(forcing))"
println("  $(nrow(forcing)) days loaded")

# ── 2. Load initial soil moisture ─────────────────────────────────────────────
println("\nReading initial_sm.csv...")
initial_sm = DataFrame(CSV.File(joinpath(outdir, "initial_sm.csv"))).moisture
println("  Initial soil moisture (m³/m³): $initial_sm")

# ── 3. Write NicheMapR input CSVs ─────────────────────────────────────────────
println("\nWriting NicheMapR input files...")

# nmr_forcing.csv — same forcing as Julia, in the column names run_nmr.R expects
CSV.write(joinpath(outdir, "nmr_forcing.csv"), DataFrame(
    doy       = forcing.DOY,
    Tmin_C    = forcing.TMINN,
    Tmax_C    = forcing.TMAXX,
    RHmin_pct = forcing.RHMINN,
    RHmax_pct = forcing.RHMAXX,
    Wind_ms   = forcing.WNMAXX,   # run_nmr.R uses this for both WNMINN and WNMAXX
    Rain_mm   = forcing.RAINFALL,
    CCmax_pct = forcing.CCMAXX,
    CCmin_pct = zeros(NDAYS),
    tannul_C  = fill(6.43, NDAYS),
))
println("  Wrote nmr_forcing.csv")

# nmr_params.csv — scalar site and model parameters
CSV.write(joinpath(outdir, "nmr_params.csv"), DataFrame(
    latitude    = [39.1368],
    longitude   = [-111.5580],
    elevation_m = [2435.35],
    dstart      = ["01/01/2013"],
    dfinish     = ["31/12/2013"],
    nyears      = [1],
    snowtemp    = [1.5],
    snowdens    = [0.375],
    snowmelt    = [1.0],
    undercatch  = [1.0],
    rainmelt    = [0.0125],
    densfun1    = [0.5979],
    densfun2    = [0.2178],
    densfun3    = [0.001],
    densfun4    = [0.0038],
    snowcond    = [0.0],
    intercept   = [0.0],
    REFL        = [0.15],
    SLE         = [0.95],
    RUF         = [0.004],
))
println("  Wrote nmr_params.csv")

# nmr_soil.csv — 10 NMR coarse thermal nodes
CSV.write(joinpath(outdir, "nmr_soil.csv"), DataFrame(
    DEP            = [0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200],
    Thcond         = _MINERAL_CONDUCTIVITY[_NMR_IDX],
    SpecHeat       = _MINERAL_HEAT_CAPACITY[_NMR_IDX],
    BulkDensity    = fill(1.3, 10),
    MineralDensity = fill(2.56, 10),
))
println("  Wrote nmr_soil.csv")

# nmr_soil19.csv — 19-node hydraulic properties
CSV.write(joinpath(outdir, "nmr_soil19.csv"), DataFrame(
    node = 1:19,
    PE   = fill(1.1,    19),
    KS   = fill(0.0037, 19),
    BB   = fill(4.5,    19),
    BD   = fill(1.3,    19),
    DD   = fill(2.56,   19),
    L    = _ROOT_DENSITY,
))
println("  Wrote nmr_soil19.csv")

# nmr_initial.csv — initial T (°C) and SM (m³/m³) at 10 NMR nodes
init_row = DataFrame(
    reshape(vcat(_INIT_ST, Float64.(initial_sm)), 1, 20),
    vcat(["ST$i" for i in 1:10], ["SM$i" for i in 1:10]),
)
CSV.write(joinpath(outdir, "nmr_initial.csv"), init_row)
println("  Wrote nmr_initial.csv")

# ── 4. Run NicheMapR via R ────────────────────────────────────────────────────
println("\nCalling run_nmr.R...")
rscript = joinpath(testdir, "R", "NicheMapR_micro_testrun_scan329.R")
run(`Rscript $rscript $outdir`)
println("  R run complete.")

# ── 5. Read NicheMapR outputs and save nmr_hourly.csv ────────────────────────
println("\nReading NicheMapR output files...")
metout = DataFrame(CSV.File(joinpath(outdir, "metout.csv"),     normalizenames=true))
soil   = DataFrame(CSV.File(joinpath(outdir, "soil.csv"),       normalizenames=true))
smoist = DataFrame(CSV.File(joinpath(outdir, "soilmoist.csv"),  normalizenames=true))
println("  metout rows: $(nrow(metout)),  soil rows: $(nrow(soil)),  soilmoist rows: $(nrow(smoist))")

nhours = NDAYS * 24
nmr = DataFrame(
    DOY      = Int.(metout.DOY[1:nhours]),
    HOUR     = Int.(round.(metout.TIME[1:nhours] ./ 60.0)),
    SNOWDEP  = round.(Float64.(metout.SNOWDEP[1:nhours]),  digits=1),
    SNOWDENS = round.(Float64.(metout.SNOWDENS[1:nhours]), digits=3),
    D5cm     = round.(Float64.(soil.D5cm[1:nhours]),       digits=1),
    D10cm    = round.(Float64.(soil.D10cm[1:nhours]),      digits=1),
    D20cm    = round.(Float64.(soil.D20cm[1:nhours]),      digits=1),
    D50cm    = round.(Float64.(soil.D50cm[1:nhours]),      digits=1),
    D200cm   = round.(Float64.(soil.D200cm[1:nhours]),     digits=1),
    WC5cm    = round.(Float64.(smoist.WC5cm[1:nhours]),    digits=3),
    WC10cm   = round.(Float64.(smoist.WC10cm[1:nhours]),   digits=3),
    WC20cm   = round.(Float64.(smoist.WC20cm[1:nhours]),   digits=3),
    WC50cm   = round.(Float64.(smoist.WC50cm[1:nhours]),   digits=3),
)
CSV.write(joinpath(outdir, "nmr_hourly.csv"), nmr)
println("  Wrote nmr_hourly.csv ($(nrow(nmr)) rows × $(ncol(nmr)) cols)")

println("\nDone.")
