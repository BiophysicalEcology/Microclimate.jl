# prepare_scan329_data.jl
#
# One-time data preparation script for the SNOTEL 329 snow comparison test.
# Run manually from the repo root:
#   julia test/prepare_scan329_data.jl
#
# Reads forcing inputs and NicheMapR outputs from the external SCAN_SNOTEL_TEST
# directory and saves compact CSV files to test/data/scan329/.
# Also reads observations from the SNOTEL 329 data file and subsets to the
# simulation period (2013-01-01 to 2016-12-31).

using CSV, DataFrames, Dates

const DATADIR  = "C:/Users/mrke/Dropbox/Completed Research Projects/SCAN_SNOTEL_TEST"
const INPUTDIR = joinpath(DATADIR, "micro csv input")

testdir = realpath(joinpath(@__DIR__))
outdir  = joinpath(testdir, "data", "scan329")
mkpath(outdir)
println("Writing output to: $outdir")

# ── 1. Daily forcing ──────────────────────────────────────────────────────────
println("\nReading daily forcing inputs...")

doy    = DataFrame(CSV.File(joinpath(INPUTDIR, "doy.csv")))[:, 2]
tminn  = DataFrame(CSV.File(joinpath(INPUTDIR, "TMINN.csv")))[:, 2]
tmaxx  = DataFrame(CSV.File(joinpath(INPUTDIR, "TMAXX.csv")))[:, 2]
wnminn = DataFrame(CSV.File(joinpath(INPUTDIR, "WNMINN.csv")))[:, 2]
wnmaxx = DataFrame(CSV.File(joinpath(INPUTDIR, "WNMAXX.csv")))[:, 2]
rhminn = DataFrame(CSV.File(joinpath(INPUTDIR, "RHMINN.csv")))[:, 2]
rhmaxx = DataFrame(CSV.File(joinpath(INPUTDIR, "RHMAXX.csv")))[:, 2]
ccmaxx = DataFrame(CSV.File(joinpath(INPUTDIR, "CCMAXX.csv")))[:, 2]
rain   = DataFrame(CSV.File(joinpath(INPUTDIR, "rain.csv")))[:, 2]

ndays = length(doy)
println("  ndays = $ndays")

forcing = DataFrame(
    DOY      = Int.(round.(doy)),
    TMINN    = round.(tminn,   digits=2),
    TMAXX    = round.(tmaxx,   digits=2),
    WNMINN   = round.(wnminn,  digits=3),
    WNMAXX   = round.(wnmaxx,  digits=3),
    RHMINN   = round.(rhminn,  digits=2),
    RHMAXX   = round.(rhmaxx,  digits=2),
    CCMAXX   = round.(ccmaxx,  digits=2),
    RAINFALL = round.(rain,    digits=2),
)
CSV.write(joinpath(outdir, "forcing.csv"), forcing)
println("  Wrote forcing.csv ($(nrow(forcing)) rows × $(ncol(forcing)) cols)")

# ── 2. Initial soil moisture from moists.csv ──────────────────────────────────
# moists.csv: 10 rows (soil layers) × ndays columns; col 2 = day-1 initial values
println("\nReading moists.csv (large file, may be slow)...")
moists_df = DataFrame(CSV.File(joinpath(INPUTDIR, "moists.csv")))
initial_sm = round.(Float64.(moists_df[1:10, 2]), digits=4)
println("  Initial soil moisture (m³/m³):")
println("  $initial_sm")
CSV.write(joinpath(outdir, "initial_sm.csv"),
    DataFrame(; layer=1:10, moisture=initial_sm))
println("  Wrote initial_sm.csv")

# ── 3. NicheMapR hourly outputs ───────────────────────────────────────────────
# Columns needed: SNOWDEP, SNOWDENS from metout; DEP3/4/6/8 from soil; WC3/4/6/8 from soilmoist
# Depth mapping (DEP.csv = [0,2.5,5,10,15,20,30,50,100,200] cm):
#   DEP3/WC3 = 5 cm, DEP4/WC4 = 10 cm, DEP6/WC6 = 20 cm, DEP8/WC8 = 50 cm
println("\nReading NicheMapR output files...")
# R write.csv format: row index + dates + DOY + TIME(minutes) + named depth columns
# normalizenames=true handles dots in names like "WC2.5cm" → "WC2_5cm"
metout = DataFrame(CSV.File(joinpath(DATADIR, "metout.csv"),     normalizenames=true))
soil   = DataFrame(CSV.File(joinpath(DATADIR, "soil.csv"),       normalizenames=true))
smoist = DataFrame(CSV.File(joinpath(DATADIR, "soilmoist.csv"),  normalizenames=true))
println("  metout rows: $(nrow(metout)),  soil rows: $(nrow(soil)),  soilmoist rows: $(nrow(smoist))")

# ── 3b. Snow node temperatures from sunsnow.csv ───────────────────────────────
# sunsnow matrix (from OSUB.f): 11 cols = JULDAY, TIME, TT(1)..TT(8) [snow], TT(9) [soil surface]
# Nodes 1-8 are snow nodes (top→bottom), node 9 is the soil surface.
# With maxsnode=4 and thresholds [2,5,10,20]cm: inactive nodes 1-4 = surface temp,
# active nodes 5-8 = 2, 5, 10, 20 cm depth from snow surface.
sunsnow_path = joinpath(DATADIR, "sunsnow.csv")
if isfile(sunsnow_path)
    println("\nReading sunsnow.csv (snow node temperatures)...")
    sunsnow_raw = DataFrame(CSV.File(sunsnow_path, normalizenames=true))
    println("  sunsnow rows: $(nrow(sunsnow_raw)),  cols: $(ncol(sunsnow_raw))")
    println("  Column names: $(propertynames(sunsnow_raw))")

    # sunsnow.csv written by NicheMapR write.csv: col1=row_index, col2=dates (string),
    # col3=DOY, col4=TIME (minutes 0-1380), col5-13=SN1-SN9 (8 snow nodes + soil surface).
    # Column names after normalizenames=true: Column1, dates, DOY, TIME, SN1-SN9.
    println("  Reading DOY/TIME from named columns :DOY and :TIME")

    snowtemp_df = DataFrame(
        DOY  = Int.(sunsnow_raw[!, :DOY]),
        HOUR = Int.(round.(Float64.(sunsnow_raw[!, :TIME]) ./ 60.0)),  # TIME (min) → hour
    )
    # SN1-SN8 are snow nodes (top→bottom), SN9 is the soil surface temperature
    for k in 1:8
        col_sym = Symbol("SN$(k)")
        snowtemp_df[!, col_sym] = round.(Float64.(sunsnow_raw[!, col_sym]), digits=3)
    end
    # SN9 = soil surface node — rename to SSOIL
    snowtemp_df[!, :SSOIL] = round.(Float64.(sunsnow_raw[!, :SN9]), digits=3)
    CSV.write(joinpath(outdir, "nmr_snowtemp.csv"), snowtemp_df)
    println("  Wrote nmr_snowtemp.csv ($(nrow(snowtemp_df)) rows × $(ncol(snowtemp_df)) cols)")
    println("  Columns: $(propertynames(snowtemp_df))")
else
    println("\nWARNING: sunsnow.csv not found at $sunsnow_path — skipping snow node temperatures")
    println("  Add  write.csv(micro_out\$sunsnow, file.path(path, 'sunsnow.csv'))  to snow_test.R")
end

# TIME is in minutes (0,60,...,1380); convert to hour (0-23) for alignment
nmr = DataFrame(
    DOY      = Int.(metout.DOY),
    HOUR     = Int.(round.(metout.TIME ./ 60.0)),
    SNOWDEP  = round.(metout.SNOWDEP,  digits=1),
    SNOWDENS = round.(metout.SNOWDENS, digits=3),
    D5cm     = round.(soil.D5cm,       digits=1),
    D10cm    = round.(soil.D10cm,      digits=1),
    D20cm    = round.(soil.D20cm,      digits=1),
    D50cm    = round.(soil.D50cm,      digits=1),
    D200cm   = round.(soil.D200cm,     digits=1),
    WC5cm    = round.(smoist.WC5cm,    digits=3),
    WC10cm   = round.(smoist.WC10cm,   digits=3),
    WC20cm   = round.(smoist.WC20cm,   digits=3),
    WC50cm   = round.(smoist.WC50cm,   digits=3),
)
CSV.write(joinpath(outdir, "nmr_hourly.csv"), nmr)
println("  Wrote nmr_hourly.csv ($(nrow(nmr)) rows × $(ncol(nmr)) cols)")

# ── 4. SNOTEL 329 observations ────────────────────────────────────────────────
# Subset to simulation period 2013-01-01 to 2016-12-31.
# Raw units: SNWD.I = snow depth in inches, WTEQ.I = SWE in inches (both convert ×2.54 → cm)
# Soil temp STO.I_* in °C, soil moisture SMS.I_* in % volumetric
println("\nReading SNOTEL 329 observations...")
# normalizenames=true converts "STO.I_2" → "STO_I_2" etc.
obs_raw = DataFrame(CSV.File(joinpath(DATADIR, "climate/329/329.csv"),
    normalizenames=true))
obs_raw[!, :DateTime] = DateTime.(string.(obs_raw[!, :DateTime]),
    dateformat"yyyy-mm-dd HH:MM:SS")

t_start = DateTime(2013, 1, 1)
t_end   = DateTime(2016, 12, 31, 23, 0, 0)
mask    = (obs_raw[!, :DateTime] .>= t_start) .& (obs_raw[!, :DateTime] .<= t_end)

# After normalizenames: "SNWD.I"→"SNWD_I", "STO.I_2"→"STO_I_2", etc.
obs_sub = obs_raw[mask, [:DateTime, :SNWD_I, :WTEQ_I,
    :STO_I_2, :STO_I_8, :STO_I_20,
    :SMS_I_2, :SMS_I_8, :SMS_I_20]]

# Rename to depth-labelled identifiers
rename!(obs_sub,
    :SNWD_I   => :SNWD_mm,
    :WTEQ_I   => :WTEQ_in,     # inches; convert ×2.54 → cm in test script
    :STO_I_2  => :STO_5cm,
    :STO_I_8  => :STO_20cm,
    :STO_I_20 => :STO_50cm,
    :SMS_I_2  => :SMS_5cm,
    :SMS_I_8  => :SMS_20cm,
    :SMS_I_20 => :SMS_50cm,
)
CSV.write(joinpath(outdir, "obs.csv"), obs_sub)
println("  Wrote obs.csv ($(nrow(obs_sub)) rows)")

println("\nDone. Data files written to $outdir")
println("Hard-wire initial_soil_moisture in snow_obs_comparison.jl:")
println("  $initial_sm")
