using DifferentialEquations
using Interpolations


days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
#days = [196]
hours = collect(0.:1:24.) # hour of day
lat = -40° # latitude
elev = 0u"m" # elevation (m)
hori = fill(0.0°, 24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
slope = 90.0° # slope (degrees, range 0-90)
aspect = 0.0° # aspect (degrees, 0 = North, range 0-360)
refl = 0.10 # substrate solar reflectivity (decimal %)
iuv = false
# Time varying environmental data
TIMINS = [0, 0, 1, 1] # time of minima for air temp, wind, humidity and cloud cover (h), air & wind mins relative to sunrise, humidity and cloud cover mins relative to solar noon
TIMAXS = [1, 1, 0, 0] # time of maxima for air temp, wind, humidity and cloud cover (h), air temp & wind maxs relative to solar noon, humidity and cloud cover maxs relative to sunrise
TMINN = [-14.3, -12.1, -5.1, 1.2, 6.9, 12.3, 15.2, 13.6, 8.9, 3, -3.2, -10.6]u"°C" # minimum air temperatures (°C)
TMAXX = [-3.2, 0.1, 6.8, 14.6, 21.3, 26.4, 29, 27.7, 23.3, 16.6, 7.8, -0.4]u"°C" # maximum air temperatures (°C)
RHMINN = [50.2, 48.4, 48.7, 40.8, 40, 42.1, 45.5, 47.3, 47.6, 45, 51.3, 52.8] # min relative humidity (%)
RHMAXX = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100] # max relative humidity (%)
WNMINN = [4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8]u"m/s" # min wind speed (m/s)
WNMAXX = [4.9, 4.8, 5.2, 5.3, 4.6, 4.3, 3.8, 3.7, 4, 4.6, 4.9, 4.8]u"m/s" # max wind speed (m/s)
CCMINN = [50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1] # min cloud cover (%)
CCMAXX = [50.3, 47, 48.2, 47.5, 40.9, 35.7, 34.1, 36.6, 42.6, 48.4, 61.1, 60.1] # max cloud cover (%)
RAINFALL = ([28, 28.2, 54.6, 79.7, 81.3, 100.1, 101.3, 102.5, 89.7, 62.4, 54.9, 41.2]) / 1000u"m" # monthly mean rainfall (mm)
daily = false

# compute solar radiation
solrad_out = solrad(days = days, hours = hours, lat = lat, elev = elev, hori = hori, slope = slope, aspect = aspect, refl = refl, iuv = iuv)

# interpolate air temperature to hourly
TAIRs, WNs, RHs, CLDs = hourly_vars(TMINN=TMINN, TMAXX=TMAXX, WNMINN=WNMINN, WNMAXX=WNMAXX, RHMINN=RHMINN, RHMAXX=RHMAXX, CCMINN=CCMINN, CCMAXX=CCMAXX,solrad_out=solrad_out, TIMINS=TIMINS, TIMAXS=TIMAXS, daily=daily)
plot(TAIRs)
plot(WNs)
plot(RHs)
plot(CLDs)

# simulate a day
iday = 2
ndays = length(days)
sub = (iday*25-24):(iday*25)

# get today's weather
SOLR = solrad_out.Global[sub1]
ZENR = solrad_out.Zenith[sub1]
ZSL = solrad_out.ZenithSlope[sub1]
TAIR = TAIRs[sub]
VEL = WNs[sub]
RH = RHs[sub]
CLD = CLDs[sub]

# create forcing weather variable splines
tspan = 0.:60:1440
tmin = tspan .* u"minute"
interpSOLR = interpolate(SOLR, BSpline(Cubic(Line(OnGrid()))))
interpZENR = interpolate(ZENR, BSpline(Cubic(Line(OnGrid()))))
interpZSL = interpolate(ZSL, BSpline(Cubic(Line(OnGrid()))))
interpTAIR = interpolate(u"K".(TAIR), BSpline(Cubic(Line(OnGrid()))))
interpVEL = interpolate(VEL, BSpline(Cubic(Line(OnGrid()))))
interpRH = interpolate(RH, BSpline(Cubic(Line(OnGrid()))))
interpCLD = interpolate(CLD, BSpline(Cubic(Line(OnGrid()))))
SOLRt = scale(interpSOLR, tspan)
ZENRt = scale(interpZENR, tspan)
ZSLt = scale(interpZSL, tspan)
TAIRt = scale(interpTAIR, tspan)
VELt = scale(interpVEL, tspan)
RHt = scale(interpRH, tspan)
CLDt = scale(interpCLD, tspan)

plot(tspan, SOLRt(tspan))
plot(tspan, ZENRt(tspan))
plot(tspan, TAIRt(tspan))
plot(tspan, VELt(tspan))
plot(tspan, RHt(tspan))
plot(tspan, CLDt(tspan))


# Initial conditions
T0 = [15.0, 14.0].*u"°C"

# Parameters
params = SoilParams(
    density = [1300.0, 1400.0],
    spheat = [1000.0, 1000.0],
    thconduct = [0.3, 0.4],
    depp = [0.0, 0.05, 0.20],
    n = 2,
    slope = 10.0,
    shayd = 30.0,
    viewf = 0.95,
    altt = 500.0,
    refls = [0.15],
    sles = [0.95],
    sigp = 5.67e-8,
    irmode = 0,
    runsnow = 0
)

tspan = (0.0, 48.0)  # 2 days
prob = ODEProblem(soil_energy_balance!, T0, tspan, params)
sol = solve(prob, Tsit5())

plot(sol, xlabel="Time (hr)", ylabel="Soil Temperature (°C)", label=["Top layer" "Bottom layer"], lw=2)










# Parameters for the simplified model
β = 1.0    # absorption coefficient
γ = 0.5    # scattering coefficient
S = 0.1    # source term

function radiative_transfer!(du, u, τ, p)
    I⁺, I⁻ = u
    du[1] = -β * I⁺ + γ * I⁻ + S      # dI⁺/dτ
    du[2] = β * I⁻ - γ * I⁺ + S       # dI⁻/dτ
end

# Initial conditions at τ = 0
u0 = [1.0, 0.0]  # e.g., I⁺(0)=1 (incoming), I⁻(0)=0 (no reflection initially)
τspan = (0.0, 1.0)  # Optical depth from 0 to 1

# Solve the system
prob = ODEProblem(radiative_transfer!, u0, τspan)
sol = solve(prob, Tsit5())

# Plot the result
using Plots
plot(sol, labels=["I⁺(τ)" "I⁻(τ)"], xlabel="Optical Depth τ", ylabel="Intensity", lw=2)


Base.@kwdef struct SoilParams
    density::Vector{Float64}
    spheat::Vector{Float64}
    thconduct::Vector{Float64}
    depp::Vector{Float64}
    n::Int
    slope::Float64
    shayd::Float64
    viewf::Float64
    altt::Float64
    refl::Float64
    sle::Float64
    sigp::Float64
end

dT = zeros(Float64, N)
function soil_energy_balance!(dT, T, p::SoilParams, t)
    N = p.n
    dT .= 0.0

    # Get environmental data at time t
    tair = TAIRt(t)
    zenr = ZENRt(t)
    solr = SOLRt(t)
    cloud = CLDt(t)
    rh = RHt(t)
    zslr = ZSLt(t)

    # Surface properties
    sabnew = 1.0 - p.refls[1]#[Int(floor(t))]  # simplistic time index
    sle = p.sles[1]#[Int(floor(t))]
    refl = p.refls[1]#[Int(floor(t))]

    # Compute soil layer properties
    wc = zeros(N)
    c = zeros(N)
    for i in 1:N
        rcsp = p.density[i] * p.spheat[i]
        if i == 1
            wc[i] = rcsp * p.depp[2] / 2
            sok = p.thconduct[1]
            c[i] = sok / p.depp[2]
        else
            wc[i] = rcsp * (p.depp[i+1] - p.depp[i-1]) / 2
            sok = p.thconduct[i]
            c[i] = sok / (p.depp[i+1] - p.depp[i])
        end
    end

    # Solar radiation
    qsolar = sabnew * solr * ((100 - p.shayd) / 100)
    if p.slope > 0 && zenr < 90
        cz = cosd(zenr)
        czsl = cosd(zslr)
        qsolar = (qsolar / cz) * czsl
    end

    # Longwave radiation (handle both IR modes)
# Constants
SIGP = 5.67e-8  # Stefan-Boltzmann constant (W/m²/K⁴)

# Atmospheric radiation (based on empirical model)
arad = (0.0000092 * (tair + 273.16)^2) * SIGP * (tair + 273.16)^4 * 60 / (4.185 * 10000)

# Cloud radiation temperature (shade approximation, TAIR - 2°C)
crad = SIGP * SLEP * (tair + 271)^4

# Hillshade radiation temperature (approximated as air temperature)
hrad = SIGP * SLEP * (tair + 273)^4

# Ground surface radiation temperature
srad = SIGP * SLE * (T[1] + 273)^4

# Clear sky fraction
clr = 1.0 - cloud / 100
    clod = crad * (cloud / 100)
    qradsk = (arad * clr + clod) * ((100 - p.shayd) / 100)
    qradvg = (p.shayd / 100) * hrad
    qradgr = ((100 - p.shayd) / 100) * srad + (p.shayd / 100) * hrad
    qradhl = hrad
    qrad = (qradsk + qradvg) * p.viewf + qradhl * (1.0 - p.viewf) - qradgr

    # get convection
    profile_out = get_profile(D0cm=T[1], TAREF=u"°C"(tair), ZEN=zenr, heights=[0.01].*u"m", RH = rh) # Stub
    qconv = profile_out.QCONV
    # Air properties
    bp = get_pressure(p.altt)
    dry_air_out = dry_air(u"K"(TAREF), elev = elev)
    wet_air_out = wet_air(u"K"(TAREF), rh=RH)
    e, cp, denair = wetair(tair, rh, bp) # Stub

    # Convection and evaporation
    hc = max(abs((qconv * 4.184 / 60 * 10000) / (T[1] - tair)), 0.5)
    hd = (hc / (cp * denair)) * (0.71 / 0.60)^0.666
    qevap = evap(T[1], tair, rh, hd) # Stub

    # Energy balance at surface
    dT[1] = (qsolar + qrad + qconv - qevap) / wc[1]

    # Soil conduction for internal nodes
    for i in 2:N-1
        dT[i] = (c[i-1] * (T[i-1] - T[i]) + c[i] * (T[i+1] - T[i])) / wc[i]
    end

    # Lower boundary condition
    dT[N] = 0.0  # or set T[N] = T_surface from data
end