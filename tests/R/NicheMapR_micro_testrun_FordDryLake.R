library(NicheMapR)
library(zoo)
head(SCANsites)
sitenum <- '2184' # Ford Dry Lake
site <- subset(SCANsites, id == sitenum) # subset the SCANsites dataset for Ford Dry Lake
name <- site$name # the name of the sites
Latitude <- site$lat # the latitude in decimal degrees
Longitude <- site$lon # the longitude in decimal degrees
Elevation <- site$elev / 3.28084 # elevation, converted from feet to metres
TZoffset <- site$`GMT offset` # the offset from Greenwich Mean Time, in hours
ystart <- 2015 # start year
yfinish <- 2015 # end year
nyears <- yfinish - ystart + 1 # number of years to run
weather <- SCAN_FordDryLake_2015 # make SCAN_FordDrylake_2015 supplied package data the weather input variable

write_input <- 1
writecsv <- 0 # make Fortran code write output as csv files
runshade <- 0 # run the model twice, once for each shade level (1) or just for the first shade level (0)?
runmoist <- 1 # run soil moisture model (0 = no, 1 = yes)?
snowmodel <- 0 # run the snow model (0 = no, 1 = yes)? - note that this runs slower
hourly <- 1 # run the model with hourly input data
rainhourly <- 0 # run the model with hourly rainfall input data (irrelevant if hourly = 1)
microdaily <- 1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day
IR <- 0 # compute clear-sky longwave radiation using Campbell and Norman (1998) eq. 10.10 (includes humidity)
solonly <- 0 # Only run SOLRAD to get solar radiation? 1=yes, 0=no
message <- 0 # do not allow the Fortran integrator to output warnings
fail <- 24*365 # how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)

longlat <- c(Longitude, Latitude) # decimal degrees longitude and latitude from the SCAN site data table
doynum <- floor(nrow(weather) / 24) # number of days to run, determined by counting the number of rows in the weather dataset and dividing by 24 to get days, but keeping it as a whole number
idayst <- 1 # start day
ida <- doynum # end day
HEMIS <- ifelse(longlat[2] < 0, 2, 1) # chose hemisphere based on latitude
ALAT <- abs(trunc(longlat[2])) # degrees latitude
AMINUT <- (abs(longlat[2]) - ALAT) * 60 # minutes latitude
ALONG <- abs(trunc(longlat[1])) # degrees longitude
ALMINT <- (abs(longlat[1]) - ALONG) * 60 # minutes latitude
ALREF <- ALONG # reference longitude for time zone

EC <- 0.0167238 # Eccenricity of the earth's orbit (current value 0.0167238, ranges between 0.0034 to 0.058)
RUF <- 0.004 # Roughness height (m), , e.g. sand is 0.0005, grass may be 0.02, current allowed range: 0.00001 (snow) - 0.02 cm.
ZH <- 0 # heat transfer roughness height (m) for Campbell and Norman air temperature/wind speed profile (invoked if greater than 1, 0.02 * canopy height in m if unknown)
D0 <- 0 # zero plane displacement correction factor (m) for Campbell and Norman air temperature/wind speed profile (0.6 * canopy height in m if unknown)
# Next four parameters are segmented velocity profiles due to bushes, rocks etc. on the surface
#IF NO EXPERIMENTAL WIND PROFILE DATA SET ALL THESE TO ZERO! (then roughness height is based on the parameter RUF)
Z01 <- 0 # Top (1st) segment roughness height(m)
Z02 <- 0 # 2nd segment roughness height(m)
ZH1 <- 0 # Top of (1st) segment, height above surface(m)
ZH2 <- 0 # 2nd segment, height above surface(m)
SLE <- 0.96 # Substrate longwave IR emissivity (decimal %), typically close to 1
ERR <- 1 # Integrator error for soil temperature calculations
Refhyt <- 2 # Reference height (m), reference height at which air temperature, wind speed and relative humidity input data are measured
DEP <- c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200) # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
Thcond <- 2.5 # soil minerals thermal conductivity (W/mC)
Density <- 2.56 # soil minerals density (Mg/m3)
SpecHeat <- 870 # soil minerals specific heat (J/kg-K)
BulkDensity <- 1.3 # soil bulk density (Mg/m3)
SatWater <- 0.26 # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
REFL <- 0.20 # soil reflectance (decimal %)
ALTT <- Elevation # altitude (m)
slope <- 0 # slope (degrees, range 0-90)
azmuth <- 180 # aspect (degrees, 0 = North, range 0-360)
hori <- rep(0, 24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
VIEWF <- 1 - sum(sin(hori * pi / 180)) / length(hori) # convert horizon angles to radians and calc view factor(s)
lamb <- 0 # Return wavelength-specific solar radiation output?
IUV <- 0 # Use gamma function for scattered solar radiation? (computationally intensive)
PCTWET <- 0 # percentage of surface area acting as a free water surface (%)
CMH2O <- 1 # precipitable cm H2O in air column, 0.1 = VERY DRY; 1.0 = MOIST AIR CONDITIONS; 2.0 = HUMID, TROPICAL CONDITIONS (note this is for the whole atmospheric profile, not just near the ground)
TIMAXS <- c(1, 1, 0, 0)   # Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise
TIMINS <- c(0, 0, 1, 1)   # Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon
minshade <- 0 # minimum available shade (%)
maxshade <- 90 # maximum available shade (%)
Usrhyt <- 0.01# local height (m) at which air temperature, relative humidity and wind speed calculatinos will be made
# Aerosol profile using GADS
relhum <- 1
optdep.summer <- as.data.frame(rungads(longlat[2],longlat[1],relhum, 0))
optdep.winter <- as.data.frame(rungads(longlat[2],longlat[1],relhum, 1))
optdep <- cbind(optdep.winter[,1],rowMeans(cbind(optdep.summer[,2],optdep.winter[,2])))
optdep <- as.data.frame(optdep)
colnames(optdep)<-c("LAMBDA","OPTDEPTH")
a <- lm(OPTDEPTH~poly(LAMBDA, 6, raw = TRUE),data = optdep)
LAMBDA <- c(290, 295, 300, 305, 310, 315, 320, 330, 340, 350, 360, 370, 380, 390, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020, 1080, 1100, 1120, 1140, 1160, 1180, 1200, 1220, 1240, 1260, 1280, 1300, 1320, 1380, 1400, 1420, 1440, 1460, 1480, 1500, 1540, 1580, 1600, 1620, 1640, 1660, 1700, 1720, 1780, 1800, 1860, 1900, 1950, 2000, 2020, 2050, 2100, 2120, 2150, 2200, 2260, 2300, 2320, 2350, 2380, 2400, 2420, 2450, 2490, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000)
TAI <- predict(a, data.frame(LAMBDA))


if(is.na(weather$TAVG.H[1])==TRUE){ # mean hourly air temperature
  weather$TAVG.H[1]<-weather$TAVG.H[which(!is.na(weather$TAVG.H))[1]]
}
if(is.na(weather$WSPDV.H[1])==TRUE){ # mean hourly wind speed
  weather$WSPDV.H[1]<-weather$WSPDV.H[which(!is.na(weather$WSPDV.H))[1]]
}
if(is.na(weather$RHUM[1])==TRUE){ # mean hourly relative humidity
  weather$RHUM[1]<-weather$RHUM[which(!is.na(weather$RHUM))[1]]
}
if(is.na(weather$SRADV.H[1])==TRUE){ # mean hourly solar radiation
  weather$SRADV.H[1]<-weather$SRADV.H[which(!is.na(weather$SRADV.H))[1]]
}
if(is.na(weather$PRCP.H[1])==TRUE){ # hourly precipitation
  weather$PRCP.H[1]<-weather$PRCP.H[which(!is.na(weather$PRCP.H))[1]]
}
weather$TAVG.H <- weather$TAVG.H
# use na.approx function from zoo package to fill in missing data
TAIRhr <- weather$TAVG.H <- na.approx(weather$TAVG.H)
RHhr <- weather$RHUM.I <- na.approx(weather$RHUM.I)
SOLRhr <- weather$SRADV.H <- na.approx(weather$SRADV.H)
RAINhr <- weather$PRCP.H <- na.approx(weather$PRCP.H * 25.4) # convert rainfall from inches to mm
WNhr <- weather$WSPDV.H <- na.approx(weather$WSPDV.H * 0.44704) # convert wind speed from miles/hour to m/s
ZENhr <- TAIRhr * 0 - 1 # negative zenith angles to force model to compute them
IRDhr <- TAIRhr * 0 - 1 #

micro <- micro_global(loc = c(Longitude, Latitude), timeinterval = 365, clearsky = 1)
# append dates
metout <- as.data.frame(micro$metout)
tzone <- paste0("Etc/GMT", TZoffset)
dates <- seq(ISOdate(ystart, 1, 1, tz = tzone)-3600 * 12, ISOdate((ystart+nyears),1, 1, tz = tzone)-3600 * 13, by="hours")
clear <- as.data.frame(cbind(dates, as.data.frame(rep(micro$metout[1:(365 * 24),13],nyears))),stringsAsFactors = FALSE)
doy <- rep(seq(1, 365),nyears)[1:floor(nrow(weather)/24)] # days of year to run
clear <- as.data.frame(clear, stringsAsFactors = FALSE)
colnames(clear)=c("datetime", "sol")

# find the maximum observed solar and adjust the clear sky prediction to this value
maxsol <- max(SOLRhr)
clear2 <- clear[, 2]*(maxsol / max(clear[, 2])) # get ratio of max observed to predicted max clear sky solar

# compute cloud cover from ratio of max to observed solar
sol <- SOLRhr
clr <- clear2
zenthresh <- 85
sol[metout$ZEN > zenthresh] <- NA # remove values for nighttime (and very early morning/late afternoon)
clr[metout$ZEN > zenthresh] <- NA # remove values for nighttime (and very early morning/late afternoon)
a <- ((clr - sol) / clr) * 100 # get ratio of observed to predicted solar, convert to %
a[a > 100] <- 100 # cap max 100%
a[a < 0] <- 0 # cap min at 0%
a[is.na(a)] <- 0 # replace NA with zero
a[is.infinite(a)]=0 # replace infinity with zero
a[a == 0] <- NA # change all zeros to NA for na.approx
a <- na.approx(a, na.rm = FALSE) # apply na.approx, but leave any trailing NAs
a[is.na(a)] <- 0 # make trailing NAs zero
CLDhr <- a # now we have hourly cloud cover

CCMAXX <- aggregate(CLDhr, by = list(weather$Date), FUN = max)[,2]#c(100, 100) # max cloud cover (%)
CCMINN <- aggregate(CLDhr, by = list(weather$Date), FUN = min)[,2]#c(0, 15.62) # min cloud cover (%)
TMAXX <- aggregate(TAIRhr, by = list(weather$Date), FUN = max)[,2]#c(40.1, 31.6) # maximum air temperatures (°C)
TMINN <- aggregate(TAIRhr, by = list(weather$Date), FUN = min)[,2]#c(19.01, 19.57) # minimum air temperatures (°C)
RAINFALL <- aggregate(RAINhr, by = list(weather$Date), FUN = sum)[,2]#c(19.01, 19.57) # minimum air temperatures (°C)
RHMAXX <- aggregate(RHhr, by = list(weather$Date), FUN = max)[,2]#c(90.16, 80.92) # max relative humidity (%)
RHMINN <- aggregate(RHhr, by = list(weather$Date), FUN = min)[,2]#c(11.05, 27.9) # min relative humidity (%)
WNMAXX <- aggregate(WNhr, by = list(weather$Date), FUN = max)[,2]#c(1.35, 2.0) # max wind speed (m/s)
WNMINN <- aggregate(WNhr, by = list(weather$Date), FUN = min)[,2]#c(0.485, 0.610) # min wind speed (m/s)

#tannul <- mean(c(TMAXX, TMINN)) # annual mean temperature for getting monthly deep soil temperature (°C)
tannul <- mean(TAIRhr) # annual mean temperature for getting monthly deep soil temperature (°C)
tannulrun <- rep(tannul, doynum) # monthly deep soil temperature (2m) (°C)
# creating the arrays of environmental variables that are assumed not to change with month for this simulation
MAXSHADES <- rep(maxshade, doynum) # daily max shade (%)
MINSHADES <- rep(minshade, doynum) # daily min shade (%)
SLES <- rep(SLE, doynum) # set up vector of ground emissivities for each day
REFLS <- rep(REFL, doynum) # set up vector of soil reflectances for each day
PCTWET <- rep(PCTWET, doynum) # set up vector of soil wetness for each day

# set up a profile of soil properites with depth for each day to be run
Numtyps <- 10 # number of soil types
Nodes <- matrix(data = 0, nrow = 10, ncol = doynum) # array of all possible soil nodes for max time span of 20 years
Nodes[1, 1:doynum]<-10 # deepest node for first substrate type

# now make the soil properties matrix
# columns are:
#1) bulk density (Mg/m3)
#2) volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
#3) thermal conductivity (W/mK)
#4) specific heat capacity (J/kg-K)
#5) mineral density (Mg/m3)
soilprops <- matrix(data = 0, nrow = 10, ncol = 5) # create an empty soil properties matrix
soilprops[, 1]<-BulkDensity # insert soil bulk density to profile 1
soilprops[, 2]<-SatWater # insert saturated water content to profile 1
soilprops[, 3]<-Thcond # insert thermal conductivity to profile 1
soilprops[, 4]<-SpecHeat # insert specific heat to profile 1
soilprops[, 5]<-Density # insert mineral density to profile 1
soilinit <- rep(tannul, 20) # make initial soil temps equal to mean annual

#use Campbell and Norman Table 9.1 soil moisture properties
soiltype <- 3 # 3 = sandy loam
PE <- rep(CampNormTbl9_1[soiltype, 4],19) #air entry potential J/kg
KS <- rep(CampNormTbl9_1[soiltype, 6],19) #saturated conductivity, kg s/m3
BB <- rep(CampNormTbl9_1[soiltype, 5],19) #soil 'b' parameter
BD <- rep(BulkDensity, 19) # soil bulk density, Mg/m3
DD <- rep(Density, 19) # soil mineral density, Mg/m3
soiltype <- 5 # change deeper nodes to 5 = a silt loam
PE[10:19]<-CampNormTbl9_1[soiltype, 4] #air entry potential J/kg
KS[10:19]<-CampNormTbl9_1[soiltype, 6] #saturated conductivity, kg s/m3
BB[10:19]<-CampNormTbl9_1[soiltype, 5] #soil 'b' parameter

L <- c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0) * 10000 # root density at each node, mm/m3 (from Campell 1985 Soil Physics with Basic, p. 131)
R1 <- 0.001 #root radius, m}\cr\cr
RW <- 2.5e+10 #resistance per unit length of root, m3 kg-1 s-1
RL <- 2e+6 #resistance per unit length of leaf, m3 kg-1 s-1
PC <- -1500 #critical leaf water potential for stomatal closure, J kg-1
SP <- 10 #stability parameter for stomatal closure equation, -
IM <- 1e-06 #maximum allowable mass balance error, kg
MAXCOUNT <- 500 #maximum iterations for mass balance, -
LAI <- rep(0.1, doynum) # leaf area index, used to partition traspiration/evaporation from PET
rainmult <- 1 # rainfall multiplier to impose catchment
maxpool <- 10 # max depth for water pooling on the surface, mm (to account for runoff)
evenrain <- 0 # spread daily rainfall evenly across 24hrs (1) or one event at midnight (0)
SoilMoist_Init <- rep(0.2, 10) # initial soil water content for each node, m3/m3
moists <- matrix(nrow = 10, ncol = doynum, data = 0) # set up an empty vector for soil moisture values through time
moists[1:10,] <- SoilMoist_Init # insert initial soil moisture
spinup <- 0 # repeat first day 3 times for steady state
dewrain <- 0 # don't feed dew back into soil as rain
moiststep <- 360 # how many steps within the hour is soil moisture solved over
maxsurf <- 85 # what is the maximum allowable soil surface temp (for stability purposes), deg C

snowtemp <- 1.5 # temperature at which precipitation falls as snow (used for snow model)
snowdens <- 0.375 # snow density (Mg/m3)
densfun <- c(0.5979, 0.2178, 0.001, 0.0038) # slope and intercept of linear model of snow density as a function of day of year - if it is c(0, 0) then fixed density used
snowmelt <- 1 # proportion of calculated snowmelt that doesn't refreeze
undercatch <- 1 # undercatch multipier for converting rainfall to snow
rainmelt <- 0.0125 # parameter in equation from Anderson's SNOW-17 model that melts snow with rainfall as a function of air temp
snowcond <- 0 # effective snow thermal conductivity W/mC (if zero, uses inbuilt function of density)
intercept <- 0 # snow interception fraction for when there's shade (0-1)
grasshade <- 0 # if 1, means shade is removed when snow is present, because shade is cast by grass/low veg

# intertidal simulation input vector (col 1 = tide in(1)/out(0), col 2 = sea water temperature in °C, col 3 = % wet from wave splash)
tides <- matrix(data = 0, nrow = 24 * doynum, ncol = 3) # matrix for tides

microinput <- c(doynum, RUF, ERR, Usrhyt, Refhyt, Numtyps, Z01, Z02, ZH1, ZH2, idayst, ida, HEMIS, ALAT, AMINUT, ALONG, ALMINT, ALREF, slope, azmuth, ALTT, CMH2O, microdaily, tannul, EC, VIEWF, snowtemp, snowdens, snowmelt, undercatch, rainmult, runshade, runmoist, maxpool, evenrain, snowmodel, rainmelt, writecsv, densfun, hourly, rainhourly, lamb, IUV, RW, PC, RL, SP, R1, IM, MAXCOUNT, IR, message, fail, snowcond, intercept, grasshade, solonly, ZH, D0, TIMAXS, TIMINS, spinup, dewrain, moiststep, maxsurf)

# all microclimate data input list - all these variables are expected by the input argument of the fortran micro2014 subroutine
microin <- list(microinput = microinput, tides = tides, doy = doy, SLES = SLES, DEP = DEP, Nodes = Nodes, MAXSHADES = MAXSHADES, MINSHADES = MINSHADES, TMAXX = TMAXX, TMINN = TMINN, RHMAXX = RHMAXX, RHMINN = RHMINN, CCMAXX = CCMAXX, CCMINN = CCMINN, WNMAXX = WNMAXX, WNMINN = WNMINN, TAIRhr = TAIRhr, RHhr = RHhr, WNhr = WNhr, CLDhr = CLDhr, SOLRhr = SOLRhr, RAINhr = RAINhr, ZENhr = ZENhr, IRDhr = IRDhr, REFLS = REFLS, PCTWET = PCTWET, soilinit = soilinit, hori = hori, TAI = TAI, soilprops = soilprops, moists = moists, RAINFALL = RAINFALL, tannulrun = tannulrun, PE = PE, KS = KS, BB = BB, BD = BD, DD = DD, L = L, LAI = LAI)


if(write_input){
  if(dir.exists("data/init_FDL") == FALSE){
    dir.create("data/init_FDL")
  }
  write.table(as.matrix(microinput), file = "data/init_FDL/microinput.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(doy, file = "data/init_FDL/doy.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(SLES, file = "data/init_FDL/SLES.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(DEP, file = "data/init_FDL/DEP.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(Nodes, file = "data/init_FDL/Nodes.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(MAXSHADES, file = "data/init_FDL/Maxshades.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(MINSHADES, file = "data/init_FDL/Minshades.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(TIMAXS, file = "data/init_FDL/TIMAXS.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(TIMINS, file = "data/init_FDL/TIMINS.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(TMAXX, file = "data/init_FDL/TMAXX.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(TMINN, file = "data/init_FDL/TMINN.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(RHMAXX, file = "data/init_FDL/RHMAXX.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(RHMINN, file = "data/init_FDL/RHMINN.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(CCMAXX, file = "data/init_FDL/CCMAXX.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(CCMINN, file = "data/init_FDL/CCMINN.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(WNMAXX, file = "data/init_FDL/WNMAXX.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(WNMINN, file = "data/init_FDL/WNMINN.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(REFLS, file = "data/init_FDL/REFLS.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(PCTWET, file = "data/init_FDL/PCTWET.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(soilinit, file = "data/init_FDL/soilinit.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(hori, file = "data/init_FDL/hori.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(TAI, file = "data/init_FDL/TAI.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(soilprops, file="data/init_FDL/soilprop.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(moists,file="data/init_FDL/moists.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(RAINFALL,file="data/init_FDL/rain.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(tannulrun,file="data/init_FDL/tannulrun.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(PE,file="data/init_FDL/PE.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(BD,file="data/init_FDL/BD.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(DD,file="data/init_FDL/DD.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(BB,file="data/init_FDL/BB.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(KS,file="data/init_FDL/KS.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(L,file="data/init_FDL/L.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(LAI,file="data/init_FDL/LAI.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(tides,file="data/init_FDL/tides.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(TAIRhr,file="data/init_FDL/TAIRhr.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(RHhr,file="data/init_FDL/RHhr.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(WNhr,file="data/init_FDL/WNhr.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(CLDhr,file="data/init_FDL/CLDhr.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(SOLRhr,file="data/init_FDL/SOLRhr.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(RAINhr,file="data/init_FDL/RAINhr.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(ZENhr,file="data/init_FDL/ZENhr.csv", sep = ",", col.names = NA, qmethod = "double")
  write.table(IRDhr,file="data/init_FDL/IRDhr.csv", sep = ",", col.names = NA, qmethod = "double")
}

micro <- microclimate(microin) # run the model in Fortran

dates <- weather$datetime[1:nrow(micro$metout)]
metout <- as.data.frame(micro$metout) # retrieve above ground microclimatic conditions, min shade
shadmet <- as.data.frame(micro$shadmet) # retrieve above ground microclimatic conditions, max shade
soil <- as.data.frame(micro$soil) # retrieve soil temperatures, minimum shade
shadsoil <- as.data.frame(micro$shadsoil) # retrieve soil temperatures, maximum shade
soilmoist <- as.data.frame(micro$soilmoist) # retrieve soil moisture, minimum shade
shadmoist <- as.data.frame(micro$shadmoist) # retrieve soil moisture, maximum shade
humid <- as.data.frame(micro$humid) # retrieve soil humidity, minimum shade
shadhumid <- as.data.frame(micro$shadhumid) # retrieve soil humidity, maximum shade
soilpot <- as.data.frame(micro$soilpot) # retrieve soil water potential, minimum shade
shadpot <- as.data.frame(micro$shadpot) # retrieve soil water potential, maximum shade
tcond <- as.data.frame(micro$tcond)
specheat <- as.data.frame(micro$specheat)
densit <- as.data.frame(micro$densit)

# append dates
metout <- cbind(dates, metout)
shadmet <- cbind(dates, shadmet)
soil <- cbind(dates, soil)
shadsoil <- cbind(dates, shadsoil)
soilmoist <- cbind(dates, soilmoist)
shadmoist <- cbind(dates, shadmoist)
humid <- cbind(dates, humid)
shadhumid <- cbind(dates, shadhumid)
soilpot <- cbind(dates, soilpot)
shadpot <- cbind(dates, shadpot)
tcond <- cbind(dates, tcond)
specheat <- cbind(dates, specheat)
densit <- cbind(dates, densit)

tstart <- as.POSIXct("2015-01-01",format="%Y-%m-%d")
tfinish <- as.POSIXct("2015-01-31",format="%Y-%m-%d")

# set up plot parameters
par(mfrow = c(1, 1)) # set up for 6 plots in 1 columns
par(oma = c(2, 1, 2, 2) + 0.1) # margin spacing stuff
par(mar = c(3, 3, 1, 1) + 0.1) # margin spacing stuff
par(mgp = c(2, 1, 0) ) # margin spacing stuff

# plot the soil temperatures
plot(dates, soil$D0cm, type='l',ylim = c(-10, 70),xlim = c(tstart, tfinish),xaxt = "n",ylab = expression("soil temperature (" * degree * C *")"),xlab="")
points(dates, soil$D2.5cm, type='l',col = 2)
#points(weather$datetime, weather$STO.I_2, type='l',col="red")
axis.POSIXct(side = 1, x = micro$dates, at = seq(tstart, tfinish, "weeks"), format = "%d-%m",  las = 2)
abline(0, 0, lty = 2, col='light blue')

plot(dates, tcond$TC0cm, type='l',xlim = c(tstart, tfinish),xaxt = "n",ylab = expression("soil temperature (" * degree * C *")"),xlab="")
points(dates, tcond$TC2.5cm, type='l',xlim = c(tstart, tfinish), col = 2, xaxt = "n",ylab = expression("soil temperature (" * degree * C *")"),xlab="")
points(dates, tcond$TC5cm, type='l',xlim = c(tstart, tfinish), col = 3, xaxt = "n",ylab = expression("soil temperature (" * degree * C *")"),xlab="")

# tstart <- as.POSIXct("2015-01-01",format="%Y-%m-%d")
# tfinish <- as.POSIXct("2015-01-31",format="%Y-%m-%d")
# 
# # set up plot parameters
# par(mfrow = c(5, 1)) # set up for 6 plots in 1 columns
# par(oma = c(2, 1, 2, 2) + 0.1) # margin spacing stuff
# par(mar = c(3, 3, 1, 1) + 0.1) # margin spacing stuff
# par(mgp = c(2, 1, 0) ) # margin spacing stuff
# 
# # plot the soil temperatures
# plot(dates, soil$D5cm, type='l',ylim = c(-10, 70),xlim = c(tstart, tfinish),xaxt = "n",ylab = expression("soil temperature (" * degree * C *")"),xlab="")
# points(weather$datetime, weather$STO.I_2, type='l',col="red")
# axis.POSIXct(side = 1, x = micro$dates, at = seq(tstart, tfinish, "weeks"), format = "%d-%m",  las = 2)
# text(tstart, 10,"5cm",col="black",pos = 4, cex = 1.5)
# abline(0, 0, lty = 2, col='light blue')
# #points(dates, metout$SNOWDEP, type='h',col='light blue')
# plot(dates, soil$D10cm, type='l',ylim = c(-10, 70),xlim = c(tstart, tfinish),xaxt = "n",ylab = expression("soil temperature (" * degree * C *")"),xlab="")
# points(weather$datetime, weather$STO.I_4, type='l',col="red")
# axis.POSIXct(side = 1, x = micro$dates, at = seq(tstart, tfinish, "weeks"), format = "%d-%m",  las = 2)
# text(tstart, 10,"10cm",col="black",pos = 4, cex = 1.5)
# abline(0, 0, lty = 2, col='light blue')
# plot(dates, soil$D20cm, type='l',ylim = c(-10, 70),xlim = c(tstart, tfinish),xaxt = "n",ylab = expression("soil temperature (" * degree * C *")"),xlab="")
# points(weather$datetime, weather$STO.I_8, type='l',col="red")
# axis.POSIXct(side = 1, x = micro$dates, at = seq(tstart, tfinish, "weeks"), format = "%d-%m",  las = 2)
# text(tstart, 10,"20cm",col="black",pos = 4, cex = 1.5)
# abline(0, 0, lty = 2, col='light blue')
# plot(dates, soil$D50cm, type='l',ylim = c(-10, 70),xlim = c(tstart, tfinish),xaxt = "n",ylab = expression("soil temperature (" * degree * C *")"),xlab="")
# points(weather$datetime, weather$STO.I_20, type='l',col="red")
# axis.POSIXct(side = 1, x = micro$dates, at = seq(tstart, tfinish, "weeks"), format = "%d-%m",  las = 2)
# text(tstart, 10,"50cm",col="black",pos = 4, cex = 1.5)
# abline(0, 0, lty = 2, col='light blue')
# plot(dates, soil$D100cm, type='l',ylim = c(-10, 70),xlim = c(tstart, tfinish),xaxt = "n",ylab = expression("soil temperature (" * degree * C *")"),xlab="")
# points(weather$datetime, weather$STO.I_40, type='l',col="red")
# axis.POSIXct(side = 1, x = micro$dates, at = seq(tstart, tfinish, "weeks"), format = "%d-%m",  las = 2)
# abline(0, 0, lty = 2, col='light blue')
# text(tstart, 10,"100cm",col="black",pos = 4, cex = 1.5)
# mtext(site$name, outer = TRUE)
# 
# # plot the soil moisture
# plot(dates, soilmoist$WC5cm * 100, type='l',ylim = c(0, 60),xaxt = "n",xlim = c(tstart, tfinish),col="blue",ylab="soil moisture (%)",xlab="")
# axis.POSIXct(side = 1, x = micro$dates, at = seq(tstart, tfinish, "weeks"), format = "%d-%m",  las = 2)
# points(weather$datetime, weather$SMS.I_2, type='l',col="red")
# text(tstart, 40,"5cm",col="black",pos = 4, cex = 1.5)
# plot(dates, soilmoist$WC10cm * 100, type='l',ylim = c(0, 60),xaxt = "n",xlim = c(tstart, tfinish),col="blue",ylab="soil moisture (%)",xlab="")
# axis.POSIXct(side = 1, x = micro$dates, at = seq(tstart, tfinish, "weeks"), format = "%d-%m",  las = 2)
# points(weather$datetime, weather$SMS.I_4, type='l',col="red")
# text(tstart, 40,"10cm",col="black",pos = 4, cex = 1.5)
# plot(dates, soilmoist$WC20cm * 100, type='l',ylim = c(0, 60),xaxt = "n",xlim = c(tstart, tfinish),col="blue",ylab="soil moisture (%)",xlab="")
# axis.POSIXct(side = 1, x = micro$dates, at = seq(tstart, tfinish, "weeks"), format = "%d-%m",  las = 2)
# points(weather$datetime, weather$SMS.I_8, type='l',col="red")
# text(tstart, 40,"20cm",col="black",pos = 4, cex = 1.5)
# plot(dates, soilmoist$WC50cm * 100, type='l',ylim = c(0, 60),xaxt = "n",xlim = c(tstart, tfinish),col="blue",ylab="soil moisture (%)",xlab="")
# axis.POSIXct(side = 1, x = micro$dates, at = seq(tstart, tfinish, "weeks"), format = "%d-%m",  las = 2)
# points(weather$datetime, weather$SMS.I_20, type='l',col="red")
# text(tstart, 40,"50cm",col="black",pos = 4, cex = 1.5)
# plot(dates, soilmoist$WC100cm * 100, type='l',ylim = c(0, 100),xaxt = "n",xlim = c(tstart, tfinish),col="blue",ylab="soil moisture (%)",xlab="")
# axis.POSIXct(side = 1, x = micro$dates, at = seq(tstart, tfinish, "weeks"), format = "%d-%m",  las = 2)
# points(weather$datetime, weather$SMS.I_40, type='l',col="red")
# text(tstart, 40,"100cm",col="black",pos = 4, cex = 1.5)
# mtext(site$name, outer = TRUE)

write.csv(metout, file = 'c:/git/Microclimate.jl/test/data/metout_FordDryLake.csv')
write.csv(soil, file = 'c:/git/Microclimate.jl/test/data/soil_FordDryLake.csv')
write.csv(soilmoist, file = 'c:/git/Microclimate.jl/test/data/soilmoist_FordDryLake.csv')
write.csv(soilpot, file = 'c:/git/Microclimate.jl/test/data/soilpot_FordDryLake.csv')
write.csv(tcond, file = 'c:/git/Microclimate.jl/test/data/tcond_FordDryLake.csv')
write.csv(specheat, file = 'c:/git/Microclimate.jl/test/data/specheat_FordDryLake.csv')
write.csv(densit, file = 'c:/git/Microclimate.jl/test/data/densit_FordDryLake.csv')
# write.csv(drlam, file = 'c:/git/Microclimate.jl/test/data/drlam_FordDryLake.csv')
# write.csv(drrlam, file = 'c:/git/Microclimate.jl/test/data/drrlam_FordDryLake.csv')
# write.csv(srlam, file = 'c:/git/Microclimate.jl/test/data/srlam_FordDryLake.csv')

write.csv(TAI, file = 'c:/git/Microclimate.jl/test/data/TAI_FordDryLake.csv', row.names = FALSE)
write.csv(TAIRhr, file = 'c:/git/Microclimate.jl/test/data/TAIRhr_FordDryLake.csv', row.names = FALSE)
write.csv(RHhr, file = 'c:/git/Microclimate.jl/test/data/RHhr_FordDryLake.csv', row.names = FALSE)
write.csv(SOLRhr, file = 'c:/git/Microclimate.jl/test/data/SOLRhr_FordDryLake.csv', row.names = FALSE)
write.csv(RAINhr, file = 'c:/git/Microclimate.jl/test/data/RAINhr_FordDryLake.csv', row.names = FALSE)
write.csv(WNhr, file = 'c:/git/Microclimate.jl/test/data/WNhr_FordDryLake.csv', row.names = FALSE)
write.csv(CLDhr, file = 'c:/git/Microclimate.jl/test/data/CLDhr_FordDryLake.csv', row.names = FALSE)

write.csv(CampNormTbl9_1, file = 'c:/git/Microclimate.jl/test/data/CampNormTbl9_1.csv', row.names = FALSE)
