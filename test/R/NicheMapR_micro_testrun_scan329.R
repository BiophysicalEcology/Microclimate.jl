# NicheMapR_micro_testrun_scan329.R
#
# Called by prepare_scan329_data.jl as:
#   Rscript NicheMapR_micro_testrun_scan329.R <outdir>
#
# Reads five CSV files written by Julia:
#   nmr_params.csv    — scalar model parameters
#   nmr_forcing.csv   — daily forcing (same as Julia run)
#   nmr_soil.csv      — 10-node thermal soil properties
#   nmr_soil19.csv    — 19-node soil moisture properties
#   nmr_initial.csv   — initial soil T and SM at 10 NMR nodes
#
# Calls NicheMapR::microclimate() directly (the Fortran solver) with Julia's
# identical forcing, bypassing micro_usa's internal OPeNDAP fetch.
#
# Writes to outdir:
#   metout.csv, soil.csv, soilmoist.csv, sunsnow.csv
#   nmr_timing.csv  — solver elapsed time (seconds)

library(NicheMapR)

args   <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
if (is.na(outdir) || !dir.exists(outdir)) stop("outdir not found: ", outdir)

cat("run_nmr.R: reading inputs from", outdir, "\n")

# ── 1. Read input CSVs ───────────────────────────────────────────────────────
p       <- read.csv(file.path(outdir, "nmr_params.csv"),  stringsAsFactors = FALSE)
forcing <- read.csv(file.path(outdir, "nmr_forcing.csv"), stringsAsFactors = FALSE)
soil10  <- read.csv(file.path(outdir, "nmr_soil.csv"),    stringsAsFactors = FALSE)
soil19  <- read.csv(file.path(outdir, "nmr_soil19.csv"),  stringsAsFactors = FALSE)
init    <- read.csv(file.path(outdir, "nmr_initial.csv"), stringsAsFactors = FALSE)

ndays <- nrow(forcing)
cat("  ndays =", ndays, "\n")

# ── 2. Location-derived scalars (mirroring micro_usa) ────────────────────────
lat    <- p$latitude[1]
lon    <- p$longitude[1]
ALTITUDES <- p$elevation_m[1]

x <- c(lon, lat)
HEMIS  <- ifelse(lat < 0, 2, 1)
ALAT   <- abs(trunc(lat))
AMINUT <- (abs(lat) - ALAT) * 60
ALONG  <- abs(trunc(lon))
ALMINT <- (abs(lon) - ALONG) * 60
ALREF  <- abs(trunc(lon))

# ── 3. Forcing vectors (already corrected by Julia) ──────────────────────────
TMAXX    <- as.matrix(forcing$Tmax_C)
TMINN    <- as.matrix(forcing$Tmin_C)
RHMAXX   <- forcing$RHmax_pct       # max RH (occurs at Tmin / night)
RHMINN   <- forcing$RHmin_pct       # min RH (occurs at Tmax / day)
CCMAXX   <- forcing$CCmax_pct
CCMINN   <- forcing$CCmin_pct
WNMAXX   <- forcing$Wind_ms         # already at 2 m — no height correction
WNMINN   <- forcing$Wind_ms
RAINFALL <- forcing$Rain_mm
tannulrun <- forcing$tannul_C
doy       <- forcing$doy

# ── 4. Fixed microclimate settings ───────────────────────────────────────────
RUF      <- p$RUF[1]
SLE      <- p$SLE[1]
REFL     <- p$REFL[1]
ERR      <- 1.5
Usrhyt   <- 0.01
Refhyt   <- 2.0
microdaily <- 1
EC       <- 0.0167238
CMH2O    <- 1
VIEWF    <- 1
runshade <- 0
runmoist <- 1
maxpool  <- 10000
evenrain <- 0
snowmodel <- 1
writecsv <- 0
hourly   <- 0
rainhourly <- 0
lamb     <- 0
IUV      <- 0
RW       <- 2.5e10; PC <- -1500; RL <- 2e6; SP <- 10; R1 <- 0.001
IM       <- 1e-6;   MAXCOUNT <- 500; IR <- 0; message_flag <- 0
fail     <- ndays * 24
ZH       <- 0; D0 <- 0; Z01 <- 0; Z02 <- 0; ZH1 <- 0; ZH2 <- 0
solonly  <- 0; grasshade <- 0; ndmax <- 3; maxsurf <- 85
TIMAXS   <- c(1.0, 1.0, 0.0, 0.0)
TIMINS   <- c(0, 0, 1, 1)
slope    <- 0; azmuth <- 0
idayst   <- 1; ida <- ndays
Numtyps  <- 10
LAI      <- rep(0.1, ndays)
hori     <- rep(0, 24)

# ── 5. Snow parameters ───────────────────────────────────────────────────────
snowtemp   <- p$snowtemp[1]
snowdens   <- p$snowdens[1]
snowmelt   <- p$snowmelt[1]
undercatch <- p$undercatch[1]
rainmelt   <- p$rainmelt[1]
densfun    <- c(p$densfun1[1], p$densfun2[1], p$densfun3[1], p$densfun4[1])
snowcond   <- p$snowcond[1]
intercept  <- p$intercept[1]
rainmult   <- 1

# ── 6. Soil properties ───────────────────────────────────────────────────────
DEP <- c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200)

soilprops <- matrix(data = 0, nrow = 10, ncol = 5)
soilprops[, 1] <- soil10$BulkDensity
soilprops[, 2] <- 1 - soil10$BulkDensity / soil10$MineralDensity
soilprops[soilprops[, 2] < 0.26, 2] <- 0.26
soilprops[, 3] <- soil10$Thcond       # cap already applied by Julia
soilprops[, 4] <- soil10$SpecHeat     # cap already applied by Julia
soilprops[, 5] <- soil10$MineralDensity

PE <- soil19$PE
KS <- soil19$KS
BB <- soil19$BB
BD <- soil19$BD
DD <- soil19$DD
L  <- soil19$L

# ── 7. Initial conditions ────────────────────────────────────────────────────
ST_init <- as.numeric(init[1, paste0("ST", 1:10)])
SM_init <- as.numeric(init[1, paste0("SM", 1:10)])

avetemp <- mean(c(TMAXX, TMINN))
# snowmodel=1 layout: 8 padding + 10 soil + 2 padding
soilinit <- c(ST_init, rep(avetemp, 10))
spinup   <- 0

moists <- matrix(nrow = 10, ncol = ndays, data = 0)
moists[1:10, ] <- SM_init

# ── 8. Derived daily arrays ───────────────────────────────────────────────────
REFLS   <- rep(REFL, ndays)
SLES    <- matrix(nrow = ndays, data = SLE)
MINSHADES <- rep(0, ndays)
MAXSHADES <- rep(90, ndays)

rainwet <- 1.5
soilwet <- RAINFALL
soilwet[soilwet <= rainwet] <- 0
soilwet[soilwet > 0]        <- 90
PCTWET <- pmax(soilwet, 0)

Nodes <- matrix(data = 0, nrow = 10, ncol = ndays)
Nodes[1:10, ] <- 1:10

# ── 9. GADS aerosol optical depth (run.gads=0 branch from micro_usa) ─────────
TAI <- c(0.0670358341290886,0.0662612704779235,0.065497075238002,0.0647431301168489,
  0.0639993178022531,0.0632655219571553,0.0625416272145492,0.0611230843885423,
  0.0597427855962549,0.0583998423063099,0.0570933810229656,0.0558225431259535,
  0.0545864847111214,0.0533843764318805,0.0522154033414562,0.0499736739981675,
  0.047855059159556,0.0458535417401334,0.0439633201842001,0.0421788036108921,
  0.0404946070106968,0.0389055464934382,0.0374066345877315,0.0359930755919066,
  0.0346602609764008,0.0334037648376212,0.0322193394032758,0.0311029105891739,
  0.0300505736074963,0.0290585886265337,0.0281233764818952,0.0272415144391857,
  0.0264097320081524,0.0256249068083005,0.0248840604859789,0.0241843546829336,
  0.0235230870563317,0.0228976873502544,0.0223057135186581,0.0217448478998064,
  0.0212128934421699,0.0207077699817964,0.0202275105711489,0.0197702578594144,
  0.0193342605242809,0.0189178697551836,0.0177713140039894,0.0174187914242432,
  0.0170790495503944,0.0167509836728154,0.0164335684174899,0.0161258546410128,
  0.0158269663770596,0.0155360978343254,0.0152525104459325,0.0149755299703076,
  0.0147045436435285,0.0144389973831391,0.0141783930434343,0.0134220329447663,
  0.0131772403830191,0.0129356456025128,0.0126970313213065,0.0124612184223418,
  0.0122280636204822,0.01199745718102,0.0115436048739351,0.0110993711778668,
  0.0108808815754663,0.0106648652077878,0.0104513876347606,0.0102405315676965,
  0.00982708969547694,0.00962473896278535,0.00903679230300494,0.00884767454432418,
  0.0083031278398166,0.00796072474935954,0.00755817587626185,0.00718610751850881,
  0.00704629977586921,0.00684663903049612,0.00654155580333479,0.00642947339729728,
  0.00627223096874308,0.00603955966866779,0.00580920937536261,0.00568506186880564,
  0.00563167068287251,0.00556222005081865,0.00550522989971023,0.00547395763028062,
  0.0054478983436216,0.00541823364504573,0.00539532163908382,0.00539239864119488,
  0.00541690124712384,0.00551525885358836,0.00564825853509463,0.00577220185074264,
  0.00584222986640171,0.00581645238345584,0.00566088137411449,0.00535516862329704,
  0.00489914757707667,0.00432017939770409,0.0036813032251836,0.00309019064543606,
  0.00270890436501562,0.00276446109239711,0.00356019862584603)

# ── 10. Empty hourly arrays (hourly=0) ───────────────────────────────────────
TAIRhr <- rep(0,   24 * ndays)
RHhr   <- rep(0,   24 * ndays)
WNhr   <- rep(0,   24 * ndays)
CLDhr  <- rep(0,   24 * ndays)
SOLRhr <- rep(0,   24 * ndays)
ZENhr  <- rep(-1,  24 * ndays)
IRDhr  <- rep(-1,  24 * ndays)
RAINhr <- rep(0,   24 * ndays)
tides  <- matrix(data = 0, nrow = 24 * ndays, ncol = 3)

tannul <- mean(tannulrun)

# ── 11. Assemble microinput vector (must match micro_usa line ~1296) ──────────
microinput <- c(
  ndays, RUF, ERR, Usrhyt, Refhyt, Numtyps,
  Z01, Z02, ZH1, ZH2, idayst, ida,
  HEMIS, ALAT, AMINUT, ALONG, ALMINT, ALREF,
  slope, azmuth, ALTITUDES, CMH2O, microdaily, tannul,
  EC, VIEWF, snowtemp, snowdens, snowmelt, undercatch, rainmult,
  runshade, runmoist, maxpool, evenrain, snowmodel, rainmelt, writecsv,
  densfun, hourly, rainhourly, lamb, IUV,
  RW, PC, RL, SP, R1, IM, MAXCOUNT, IR, message_flag, fail,
  snowcond, intercept, grasshade, solonly, ZH, D0,
  TIMAXS, TIMINS, spinup, 0, 360, maxsurf
)

# ── 12. Build micro list ─────────────────────────────────────────────────────
micro <- list(
  tides      = tides,
  microinput = microinput,
  doy        = doy,
  SLES       = SLES,
  DEP        = DEP,
  Nodes      = Nodes,
  MAXSHADES  = MAXSHADES,
  MINSHADES  = MINSHADES,
  TMAXX      = TMAXX,
  TMINN      = TMINN,
  RHMAXX     = RHMAXX,
  RHMINN     = RHMINN,
  CCMAXX     = CCMAXX,
  CCMINN     = CCMINN,
  WNMAXX     = WNMAXX,
  WNMINN     = WNMINN,
  TAIRhr     = TAIRhr,
  RHhr       = RHhr,
  WNhr       = WNhr,
  CLDhr      = CLDhr,
  SOLRhr     = SOLRhr,
  RAINhr     = RAINhr,
  ZENhr      = ZENhr,
  IRDhr      = IRDhr,
  REFLS      = REFLS,
  PCTWET     = PCTWET,
  soilinit   = soilinit,
  hori       = hori,
  TAI        = TAI,
  soilprops  = soilprops,
  moists     = moists,
  RAINFALL   = RAINFALL,
  tannulrun  = tannulrun,
  PE         = PE,
  KS         = KS,
  BB         = BB,
  BD         = BD,
  DD         = DD,
  L          = L,
  LAI        = LAI
)

# ── 13. Run the microclimate solver ─────────────────────────────────────────
cat("  Running NicheMapR microclimate solver (", ndays, "days)...\n")
ptm    <- proc.time()
microut <- NicheMapR::microclimate(micro)
elapsed <- (proc.time() - ptm)[["elapsed"]]
cat(sprintf("  NMR solver: %.2f s\n", elapsed))

# ── 14. Save outputs ─────────────────────────────────────────────────────────
metout   <- as.data.frame(microut$metout)
soil     <- as.data.frame(microut$soil)
soilmoist <- as.data.frame(microut$soilmoist)
sunsnow  <- as.data.frame(microut$sunsnow)

write.csv(metout,    file = file.path(outdir, "metout.csv"),    row.names = FALSE)
write.csv(soil,      file = file.path(outdir, "soil.csv"),      row.names = FALSE)
write.csv(soilmoist, file = file.path(outdir, "soilmoist.csv"), row.names = FALSE)
write.csv(sunsnow,   file = file.path(outdir, "sunsnow.csv"),   row.names = FALSE)

write.csv(data.frame(elapsed_s = elapsed),
          file = file.path(outdir, "nmr_timing.csv"), row.names = FALSE)

cat("  Outputs written to", outdir, "\n")
