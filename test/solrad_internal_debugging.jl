hours = [12.]
days = [15.]
year=2001.0 # needed to determine if a leap year
 
    lonc=0.0 # longitude correction, hours
    slope=0u"°"
    aspect=0u"°"
    #hori::Vector{typeof(0.0u"°")}=fill(0.0, 24) .* u"°"
    refl=0.1 # substrate solar reflectivity (decimal %)
    cmH2O=1 # precipitable cm H2O in air column, 0.1 = VERY DRY; 1 = MOIST AIR CONDITIONS; 2 = HUMID, TROPICAL CONDITIONS (note this is for the whole atmospheric profile, not just near the ground)
    ϵ=0.0167238
    ω=2π / 365
    se=0.39784993 #0.39779,
    d0=80.0
    iuv=false # Use gamma function for scattered solar radiation? (computationally intensive)
    noscat::Bool=true
    amr=25.0u"km"
    nmax=111 # maximum number of wavelengths
    Iλ=float.([ # wavelengths across which to integrate
        290, 295, 300, 305, 310, 315, 320, 330, 340, 350, 360, 370, 380, 390,
        400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700,
        720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020,
        1080, 1100, 1120, 1140, 1160, 1180, 1200, 1220, 1240, 1260, 1280, 1300, 1320,
        1380, 1400, 1420, 1440, 1460, 1480, 1500, 1540, 1580, 1600, 1620, 1640, 1660,
        1700, 1720, 1780, 1800, 1860, 1900, 1950, 2000, 2020, 2050, 2100, 2120, 2150,
        2200, 2260, 2300, 2320, 2350, 2380, 2400, 2420, 2450, 2490, 2500, 2600, 2700,
        2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000
    ]) * u"nm"
    OZ=reshape([ # seasonal variation of atmospheric ozone in cm, Robinson 1966 Table 4.2
            0.31, 0.31, 0.32, 0.32, 0.31, 0.3, 0.27, 0.24, 0.23, 0.22, 0.23, 0.24,
            0.27, 0.3, 0.32, 0.33, 0.34, 0.34, 0.33, 0.3, 0.31, 0.31, 0.31, 0.3, 0.29, 0.28,
            0.25, 0.24, 0.22, 0.24, 0.26, 0.28, 0.32, 0.36, 0.39, 0.4, 0.4, 0.39, 0.3, 0.31,
            0.31, 0.3, 0.29, 0.28, 0.26, 0.24, 0.24, 0.23, 0.24, 0.26, 0.29, 0.33, 0.38, 0.42,
            0.45, 0.46, 0.46, 0.27, 0.28, 0.29, 0.3, 0.3, 0.29, 0.27, 0.25, 0.24, 0.23, 0.25,
            0.27, 0.3, 0.34, 0.38, 0.4, 0.42, 0.43, 0.42, 0.34, 0.35, 0.34, 0.33, 0.32, 0.31,
            0.28, 0.25, 0.24, 0.24, 0.26, 0.28, 0.3, 0.34, 0.37, 0.39, 0.4, 0.4, 0.39, 0.38,
            0.4, 0.39, 0.38, 0.36, 0.33, 0.28, 0.25, 0.24, 0.24, 0.25, 0.27, 0.3, 0.33, 0.35,
            0.36, 0.36, 0.36, 0.34, 0.43, 0.44, 0.43, 0.41, 0.39, 0.35, 0.29, 0.25, 0.24, 0.24,
            0.25, 0.26, 0.29, 0.31, 0.33, 0.34, 0.34, 0.33, 0.32, 0.45, 0.46, 0.45, 0.42, 0.4,
            0.37, 0.31, 0.26, 0.24, 0.23, 0.24, 0.26, 0.28, 0.3, 0.31, 0.32, 0.31, 0.3, 0.3,
            0.41, 0.42, 0.43, 0.42, 0.4, 0.38, 0.32, 0.26, 0.24, 0.23, 0.24, 0.26, 0.27, 0.28,
            0.3, 0.3, 0.29, 0.28, 0.27, 0.37, 0.38, 0.4, 0.4, 0.39, 0.37, 0.32, 0.26, 0.24,
            0.22, 0.23, 0.25, 0.26, 0.27, 0.28, 0.28, 0.28, 0.27, 0.26, 0.34, 0.36, 0.38, 0.39,
            0.37, 0.34, 0.29, 0.26, 0.24, 0.22, 0.23, 0.25, 0.26, 0.28, 0.29, 0.3, 0.29, 0.29,
            0.28, 0.31, 0.32, 0.34, 0.35, 0.35, 0.32, 0.29, 0.25, 0.23, 0.22, 0.23, 0.25, 0.27,
            0.29, 0.3, 0.31, 0.31, 0.31, 0.3
        ], (19, 12))
    τR=[
        1.41, 1.31, 1.23, 1.14, 1.05, 0.99, 0.92, 0.81, 0.72, 0.63, 0.56, 0.5,
        0.45, 0.4, 0.36, 0.3, 0.25, 0.2, 0.17, 0.15, 0.12, 0.1, 0.09, 0.08, 0.07, 0.06,
        0.06, 0.05, 0.04, 0.04, 0.03, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01,
        0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ]
    τO=[
        11.5, 6.3, 3.2, 1.62, 0.83, 0.44, 0.26, 0.03, 0.02, 0.01, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ]
    τA=[0.269904738, 0.266147825, 0.262442906, 0.258789404, 0.255186744, 0.251634356, 0.248131676, 0.2412732,
        0.234606887, 0.228128378, 0.221833385, 0.215717692, 0.20977715, 0.204007681, 0.198405272, 0.187685927,
        0.177588357, 0.168082846, 0.159140695, 0.150734206, 0.142836655, 0.135422274, 0.128466227, 0.12194459,
        0.115834329, 0.110113284, 0.104760141, 0.099754417, 0.09507644, 0.090707328, 0.086628967, 0.082823998,
        0.07927579, 0.075968428, 0.072886691, 0.070016034, 0.067342571, 0.064853053, 0.062534858, 0.060375964,
        0.058364941, 0.056490925, 0.054743609, 0.053113222, 0.051590514, 0.050166738, 0.046408775, 0.045302803,
        0.044259051, 0.043271471, 0.042334415, 0.041442618, 0.040591184, 0.039775572, 0.038991583, 0.038235345,
        0.037503301, 0.036792197, 0.036099067, 0.034101935, 0.033456388, 0.032817888, 0.032184949, 0.031556287,
        0.030930816, 0.030307633, 0.029065372, 0.027825562, 0.027205981, 0.026586556, 0.025967391, 0.025348692,
        0.024114005, 0.023498886, 0.021669152, 0.021066668, 0.019292088, 0.018144698, 0.016762709, 0.015451481,
        0.014949794, 0.014224263, 0.013093462, 0.012670686, 0.012070223, 0.011164062, 0.010241734, 0.009731103,
        0.009507687, 0.009212683, 0.008965785, 0.008827751, 0.008710756, 0.008574128, 0.008462605, 0.008446967,
        0.008539475, 0.009015237, 0.009748444, 0.010586023, 0.011359647, 0.011901268, 0.012062153, 0.011735443,
        0.010882215, 0.009561062, 0.007961182, 0.006438984, 0.005558204, 0.006133532, 0.009277754
    ]
    τW=[
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0.123, 0.117, 0.1, 0.23, 0.174, 0.058, 0, 0.024, 0.027,
        0.036, 0.215, 0.25, 0.136, 0.058, 0.047, 0.036, 0.042, 0.098, 0.044, 0, 0.038,
        0.83, 0, 0, 0.38, 0.289, 0.258, 0.173, 0.008, 0, 0, 0, 0, 0, 0, 0, 0, 0.57, 0.76, 0,
        0.185, 0.291, 0.178, 0.196, 0.112, 0.075, 0.074, 0.07, 0.007, 0, 0, 0, 0.086,
        0.122, 0.132, 0.14, 0.207, 0.259, 0, 0, 0, 0.549, 0.297, 0.462, 0.52, 0.374, 0.222,
        0.614, 0.058, 0.038, 0.03, 0.04, 0.16
    ]
    Sλ=[
        48.2, 58.4, 51.4, 60.2, 68.6, 75.7, 81.9, 103.7, 105, 107.4, 105.5,
        117.3, 111.7, 109.9, 143.3, 175.8, 182.3, 208, 208.5, 194.6, 183.3, 178.3,
        169.5, 170.5, 164.6, 157.6, 151.7, 146.8, 141.8, 136.9, 131.4, 126, 120, 115,
        110.7, 105, 100, 95, 91, 88, 85, 83, 80, 77, 75, 70, 61, 59, 56, 54, 51, 49, 48, 45,
        42, 41, 40, 39, 38, 34, 33, 32, 31, 30, 29, 28, 26, 25, 24, 24, 23, 22, 21, 19, 16, 15,
        12, 11, 10.7, 10.3, 10, 9.7, 9, 8.8, 8.5, 7.9, 7.4, 6.8, 6.7, 6.6, 6.5, 6.4, 6.2,
        5.9, 5.5, 5.4, 4.8, 4.3, 3.9, 3.5, 3.1, 2.6, 2.3, 1.9, 1.7, 1.5, 1.4, 1.2, 1.1, 1, 1
    ] * u"mW * cm^-2 * nm^-1"

    FD=reshape([
            8.00e-5, 6.50e-5, 4.00e-5, 2.30e-5, 1.00e-5, 4.50e-6, 1.00e-6, 1.00e-7, 5.50e-9,
            1.00e-9, 3.50e-10, 1.60e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10, 1.00e-10,
            1.00e-10, 1.00e-10, 1.00e-3, 9.50e-4, 9.00e-4, 8.00e-4, 7.00e-4, 6.00e-4,
            4.50e-4, 3.00e-4, 1.70e-4, 8.00e-5, 3.30e-5, 1.80e-5, 1.00e-5, 8.00e-6,
            6.30e-6, 5.00e-6, 4.10e-6, 3.50e-6, 3.00e-6, 3.50e-2, 3.30e-2, 3.10e-2,
            2.90e-2, 2.50e-2, 2.20e-2, 1.75e-2, 1.30e-2, 7.50e-3, 4.50e-3, 2.50e-3,
            1.30e-3, 5.60e-4, 2.70e-4, 1.40e-4, 7.10e-5, 4.20e-5, 2.60e-5, 1.70e-5,
            1.55e-1, 1.40e-1, 1.40e-1, 1.30e-1, 1.20e-1, 1.10e-1, 1.00e-1, 9.00e-2,
            8.20e-2, 7.00e-2, 5.20e-2, 3.50e-2, 2.20e-2, 1.00e-2, 4.00e-3, 1.20e-3,
            4.50e-4, 1.90e-4, 8.00e-5, 3.70e-1, 3.70e-1, 3.60e-1, 3.50e-1, 3.30e-1,
            3.10e-1, 2.90e-1, 2.60e-1, 2.30e-1, 2.00e-1, 1.70e-1, 1.50e-1, 1.00e-1,
            7.00e-2, 3.80e-2, 1.60e-2, 5.00e-3, 1.20e-3, 3.00e-4, 5.50e-1, 5.50e-1,
            5.50e-1, 5.30e-1, 5.10e-1, 5.00e-1, 4.70e-1, 4.40e-1, 4.00e-1, 3.60e-1,
            3.10e-1, 2.70e-1, 2.25e-1, 1.80e-1, 1.20e-1, 6.00e-2, 2.50e-2, 4.50e-3,
            7.00e-4, 6.50e-1, 6.50e-1, 6.50e-1, 6.50e-1, 6.20e-1, 6.00e-1, 5.70e-1,
            5.50e-1, 5.00e-1, 4.50e-1, 4.20e-1, 3.80e-1, 3.25e-1, 2.70e-1, 1.90e-1,
            1.15e-1, 5.00e-2, 1.10e-2, 1.20e-3, 7.88e-1, 7.80e-1, 8.00e-1, 8.00e-1,
            8.00e-1, 7.60e-1, 7.35e-1, 7.10e-1, 7.00e-1, 6.50e-1, 6.00e-1, 5.50e-1,
            5.14e-1, 4.50e-1, 3.50e-1, 2.68e-1, 1.50e-1, 6.27e-2, 1.20e-2, 7.48e-1,
            7.40e-1, 7.40e-1, 7.30e-1, 7.20e-1, 7.10e-1, 7.04e-1, 6.90e-1, 6.70e-1,
            6.20e-1, 5.70e-1, 5.30e-1, 5.16e-1, 4.80e-1, 3.90e-1, 2.90e-1, 1.70e-1,
            7.62e-2, 3.00e-2, 7.00e-1, 7.00e-1, 7.00e-1, 6.90e-1, 6.80e-1, 6.80e-1,
            6.60e-1, 6.50e-1, 6.30e-1, 6.00e-1, 5.60e-1, 5.21e-1, 5.00e-1, 4.70e-1,
            3.90e-1, 3.00e-1, 1.85e-1, 9.00e-2, 2.60e-2, 6.51e-1, 6.50e-1, 6.50e-1,
            6.40e-1, 6.30e-1, 6.25e-1, 6.22e-1, 6.00e-1, 5.90e-1, 5.70e-1, 5.50e-1,
            5.20e-1, 4.89e-1, 4.60e-1, 3.90e-1, 3.08e-1, 2.00e-1, 9.55e-2, 2.20e-2
        ], (11, 19))
    FDQ=reshape([
            8.00e-6, 7.00e-6, 5.20e-6, 3.50e-6, 1.70e-6, 5.50e-7, 1.00e-7, 2.50e-8, 6.00e-9,
            1.50e-9, 3.00e-10, 6.00e-11, 1.00e-11, 1.00e-11, 1.00e-11, 1.00e-11, 1.00e-11,
            1.00e-11, 1.00e-11, 6.10e-4, 6.00e-4, 5.50e-4, 4.50e-4, 3.40e-4, 2.30e-4,
            1.20e-4, 5.50e-5, 2.60e-5, 1.20e-5, 6.00e-6, 3.50e-6, 2.00e-6, 1.50e-6,
            1.00e-6, 6.00e-7, 4.00e-7, 2.30e-7, 1.00e-7, 2.40e-2, 2.30e-2, 2.20e-2,
            2.10e-2, 1.80e-2, 1.50e-2, 1.20e-2, 7.50e-3, 4.80e-3, 2.50e-3, 1.20e-3,
            5.00e-4, 2.50e-4, 1.00e-4, 4.50e-5, 2.00e-5, 1.00e-5, 5.50e-6, 2.50e-6,
            1.30e-1, 1.20e-1, 1.10e-1, 1.00e-1, 9.10e-2, 8.50e-2, 7.20e-2, 6.50e-2,
            5.20e-2, 4.00e-2, 2.80e-2, 1.70e-2, 1.00e-2, 3.00e-3, 1.00e-3, 4.00e-4,
            1.70e-4, 6.70e-5, 2.50e-5, 3.40e-1, 3.30e-1, 3.20e-1, 3.00e-1, 2.90e-1,
            2.70e-1, 2.50e-1, 2.00e-1, 1.70e-1, 1.50e-1, 1.20e-1, 8.00e-2, 5.80e-2,
            3.00e-2, 1.30e-2, 6.20e-3, 2.00e-3, 6.00e-4, 1.90e-4, 5.40e-1, 5.30e-1,
            5.10e-1, 5.00e-1, 4.80e-1, 4.50e-1, 4.20e-1, 3.70e-1, 3.20e-1, 2.70e-1,
            2.20e-1, 1.70e-1, 1.20e-1, 8.00e-2, 5.00e-2, 2.30e-2, 8.00e-3, 4.70e-3,
            9.00e-4, 6.50e-1, 6.40e-1, 6.20e-1, 6.10e-1, 6.00e-1, 5.60e-1, 5.10e-1,
            4.70e-1, 4.20e-1, 3.60e-1, 3.00e-1, 2.50e-1, 1.80e-1, 1.30e-1, 8.00e-2,
            5.00e-2, 2.00e-2, 4.60e-3, 1.50e-3, 8.20e-1, 8.10e-1, 8.10e-1, 8.00e-1,
            7.90e-1, 7.50e-1, 7.00e-1, 6.50e-1, 6.00e-1, 5.40e-1, 4.60e-1, 3.80e-1,
            3.20e-1, 2.50e-1, 1.80e-1, 1.20e-1, 6.00e-2, 2.50e-2, 8.00e-3, 8.60e-1,
            8.40e-1, 8.20e-1, 8.00e-1, 7.50e-1, 7.00e-1, 6.80e-1, 6.30e-1, 5.80e-1,
            5.00e-1, 4.40e-1, 3.70e-1, 3.20e-1, 2.40e-1, 1.70e-1, 1.20e-1, 6.00e-2,
            2.70e-2, 1.00e-2, 8.00e-1, 7.90e-1, 7.70e-1, 7.60e-1, 7.30e-1, 6.85e-1,
            6.40e-1, 5.90e-1, 5.40e-1, 4.75e-1, 4.20e-1, 3.65e-1, 3.20e-1, 2.50e-1,
            1.80e-1, 1.30e-1, 7.00e-2, 2.90e-2, 1.10e-2, 7.50e-1, 7.40e-1, 7.30e-1,
            7.20e-1, 7.00e-1, 6.70e-1, 6.10e-1, 5.50e-1, 5.00e-1, 4.50e-1, 4.00e-1,
            3.60e-1, 3.20e-1, 2.60e-1, 1.90e-1, 1.40e-1, 8.00e-2, 3.10e-2, 1.20e-2
        ], (11, 19))
    S=[0.2, 0.255, 0.315, 0.365, 0.394, 0.405, 0.405, 0.395, 0.37, 0.343, 0.32]

     ndays = length(days)    # number of days
    ntimes = length(hours)  # number of times
    nsteps = ndays * ntimes # total time steps

    # arrays to hold every time step's radiation between 300 and 320 nm in 2 nm steps
    GRλs = fill(0.0, nsteps, nmax)u"mW/nm/cm^2" # wavelength-specific global radiation
    DRRλs = fill(0.0, nsteps, nmax)u"mW/nm/cm^2"# wavelength-specific direct Rayleigh radiation
    DRλs = fill(0.0, nsteps, nmax)u"mW/nm/cm^2" # wavelength-specific direct radiation
    SRλs = fill(0.0, nsteps, nmax)u"mW/nm/cm^2" # wavelength-specific scattered radiation
    GRs = fill(0.0, nsteps)u"mW/cm^2"           # total global radiation
    DRRs = fill(0.0, nsteps)u"mW/cm^2"          # total direct Rayleigh radiation
    DRs = fill(0.0, nsteps)u"mW/cm^2"           # total direct radiation
    SRs = fill(0.0, nsteps)u"mW/cm^2"           # total scattered radiation

    # arrays to hold zenith angles each step
    Zs = fill(90.0, nsteps)u"°"                  # zenith angles
    ZSLs = fill(90.0, nsteps)u"°"                # slope zenith angles
    HHs = fill(0.0, ndays)                      # hour angles
    tsns = fill(0.0, ndays)                     # hour at solar noon
    DOYs = Vector{Int}(undef, nsteps)           # day of year
    times = Vector{Real}(undef, nsteps)         # time

    step = 1
    HH = 0.0 # initialise sunrise hour angle
    tsn = 12.0 # initialise time of solar noon

    i=1
    GRINT = fill(0.0, nmax)u"mW/cm^2"   # integrated global radiation component (direct + scattered)
        DRRINT = fill(0.0, nmax)u"mW/cm^2"  # integrated direct Rayleigh radiation component
        DRINT = fill(0.0, nmax)u"mW/cm^2"   # integrated direct radiation component
        SRINT = fill(0.0, nmax)u"mW/cm^2"   # integrated scattered radiation component
        AIλ = fill(0.0, nmax)u"nm"
        GRλ = GRINT * u"1/nm"               # wavelength-specific global radiation component (direct + scattered)
        DRRλ = GRINT * u"1/nm"              # wavelength-specific direct Rayleigh radiation component
        DRλ = GRINT * u"1/nm"               # wavelength-specific direct radiation component
        SRλ = GRINT * u"1/nm"               # wavelength-specific scattered radiation component
j=1
 d = days[i]
            t = hours[j]
            h, tsn = hour_angle(t, lonc) # hour angle (radians)
            ζ, δ, z, AR2 = solar_geometry(d=d, lat=lat, h=h, d0=d0, ω=ω, ϵ=ϵ, se=se) # compute ecliptic, declination, zenith angle and (a/r)^2
            Z = uconvert(u"°", z)
            Zsl = Z
           #check_skylight!(z, nmax, SRINT, GRINT) # checking zenith angle for possible skylight before sunrise or after sunset
            if Z < 107.0u"°"
                if Z > 88.0u"°"
                    # Compute skylight based on G.V. Rozenberg. 1966. Twilight. Plenum Press.
                    # p. 18,19.  First computing lumens: y = b - mx. Regression of data is:
                    Elog = 41.34615384 - 0.423076923 * ustrip(u"°", Z)
                    # Converting lux (lumen/m2) to W/m2 on horizontal surface -
                    # Twilight - scattered skylight before sunrise or after sunset
                    # From p. 239 Documenta Geigy Scientific Tables. 1966. 6th ed. K. Diem, ed.
                    # Mech./elect equiv. of light = 1.46*10^-3 kW/lumen
                    Skylum = (10.0^Elog) * 1.46E-03u"mW * cm^-2"
                    SRINT[nmax] = Skylum
                    GRINT[nmax] = SRINT[nmax]
                end
            end
            
            # testing cos(h) to see if it exceeds +1 or -1
            TDTL = -tan(δ) * tan(lat) # from eq.7 McCullough & Porter 1971
            if abs(TDTL) >= 1 # long day or night
                H = π
            else
                H = abs(acos(TDTL))
            end
            # check if sunrise
            HH = 12.0 * H / π
            ts = t - tsn

            sun_up = true

            if ts <= 0.0 && abs(ts) > HH
                sun_up = false
            elseif ts > 0.0 && ts >= HH
                sun_up = false
            end

                h, tsn = hour_angle(t, lonc) # hour angle (radians)
                ζ, δ, z, AR2 = solar_geometry(d=d, lat=lat, h=h, d0=d0, ω=ω, ϵ=ϵ, se=se) # compute ecliptic, declination, zenith angle and (a/r)^2 - redundant?
                alt = (π / 2 - z)u"rad"
                altdeg = uconvert(u"°", alt).val
                cazsun = (sin(δ) - sin(lat) * sin(alt)) / (cos(lat) * cos(alt)) # cos(solar azimuth)
                #      Error checking for instability in trig function
                if cazsun < -0.9999999
                    azsun = π
                else
                    if cazsun > 0.9999999
                        azsun = 0
                    else
                        azsun = acos(cazsun)
                    end
                end
                if h <= 0
                    # Morning
                    dazsun = uconvert(u"°", azsun).val
                else
                    # Afternoon
                    if sign(lat) < 0
                        dazsun = 360u"°" - uconvert(u"°", azsun).val
                    else
                        dazsun = 180u"°" + (180u"°" - uconvert(u"°", azsun).val)
                    end
                end

                cz = cos(z)
                intcz = Int(floor(100.0 * cz + 1.0))
                Z = uconvert(u"°", z)  # zenith angle in degrees

                # horizon angle - check this works when starting at 0 rather than e.g. 15 deg
                azi = range(0u"°", stop=360u"°" - 360u"°" / length(hori), length=length(hori))
                ahoriz = hori[argmin(abs.(dazsun .- azi))]
                                
                # slope zenith angle calculation (Eq. 3.15 in Sellers 1965. Physical Climatology. U. Chicago Press)
                if slope > 0u"°"
                    czsl = cos(z) * cos(slope) + sin(z) * sin(slope) * cos(dazsun - aspect)
                    zsl = acos(czsl)
                    Zsl = min(uconvert(u"°", zsl), 90u"°") # cap at 90 degrees if sun is below slope horizon
                    intczsl = Int(floor(100.0 * czsl + 1.0))
                else
                    czsl = cz
                    zsl = z
                    Zsl = Z
                    intczsl = intcz
                end

                # refraction correction check
                if z < 1.5358896
                    # skip refraction correction
                else
                    refr = 16.0 + ((z - 1.53589) * 15) / (π / 90)
                    refr = (refr / 60) * (π / 180)
                    z -= refr
                end

                # optical air mass (Rozenberg 1966 formula p.159 in book 'Twilight') ---
                airms = 1.0 / (cos(z) + (0.025 * exp(-11.0 * cos(z))))
                cz = cos(z)
                intcz = Int(floor(100.0 * cz + 1.0))
                Z = uconvert(u"°", z)  # zenith angle in degrees

                # atmospheric ozone lookup
                # convert latitude in degrees to nearest 10-degree index
                tlat = (lat + 100.0u"°") / 10.0u"°"
                llat = Int(floor(tlat))
                allat = llat
                ala = allat + 0.5
                if tlat > ala
                    llat += 1
                end
                # clamp llat index to valid range
                mon = month(Date(year, 1, 1) + Day(d - 1)) # month from day of year
                llat = clamp(llat, 1, size(OZ, 1))
                ozone = OZ[llat, mon]  # ozone thickness (cm) from lookup table
                ELEVFCT1, ELEVFCT2, ELEVFCT3, ELEVFCT4 = elev_corr(elev)

                # deal with this:
                # c     mutliplier to correct hourly solar data for horizon angle
                #     if(altdeg.lt.ahoriz)then
                # c	   diffuse only - cut down to diffuse fraction      
                #     TDD(111+IT)=TDD(111+IT)* (0.12 + 0.83 * ((CCMINN(IDAY) + 
                #    &  CCMAXX(IDAY))/ 2. / 100.)) ! from Butt et al. 2010 Agricultural and Forest Meteorology 150 (2010) 361–368
                #     endif
                P = Microclimate.get_pressure(elev) # pressure from elevation

                for N in 1:nmax
                    τλ1 = (P / 101300u"Pa") * τR[N] * ELEVFCT1
                    τλ2 = (25.0u"km" / amr) * τA[N] * ELEVFCT2
                    τλ3 = (ozone / 0.34) * τO[N] * ELEVFCT3
                    τλ4 = τW[N] * sqrt(airms * cmH2O * ELEVFCT4)
                    τλ = ((float(τλ1) + τλ2 + τλ3) * airms) + τλ4

                    if τλ > 80.0 # making sure that at low sun angles air mass doesn't make τλ too large
                        τλ = 80.0
                    end

                    part1 = Sλ[N] * AR2 * cz
                    part2 = τλ > 0.0 ? exp(-τλ) : 0.0
                    if part2 < 1.0e-24
                        DRλ[N] = 0.0u"mW / cm^2 / nm"
                    else
                        DRλ[N] = ((ustrip(part1) * part2) / 1000.0) * u"mW / cm^2 / nm"
                    end

                    # so the integrator doesn't get confused at very low sun angles
                    if DRλ[N] < 1.0e-25u"mW / cm^2 / nm"
                        DRλ[N] = 1.0e-25u"mW / cm^2 / nm"
                    end

                    DRRλ[N] = (Sλ[N] * AR2 * cz) * exp(-float(τλ1) * airms) / 1000.0

                    if altdeg < ahoriz
                        DRλ[N] = 1.0e-25u"mW / cm^2 / nm"
                        DRRλ[N] = 1.0e-25u"mW / cm^2 / nm"
                    end

                    # Sky (SRλ) and Global Radiation (GRλ)
                    if noscat == false
                        SRλ[N] = 0.0u"mW / cm^2 / nm"
                    elseif iuv
                        if τλ1 >= 0.03
                            GAMR, GAML, SBAR = GAMMA(τλ1)
                            SRλ[N] = (
                                         ((float(GAML[intcz]) + float(GAMR[intcz])) / (2.0 * (1.0 - refl * float(SBAR))))
                                         -
                                         exp(-float(τλ1) * airms)
                                     ) * cz * Sλ[N] * AR2 / 1000.0
                        else
                            SRλ[N] = 0.0u"mW / cm^2 / nm"
                        end
                    else
                        if N > 11
                            SRλ[N] = 0.0u"mW / cm^2 / nm"
                        else
                            #I = round(Int, (ustrip(Z) + 5) / 5)
                            B = ustrip(Z) / 5
                            IA = trunc(Int, B)
                            C = B - IA
                            if C > 0.5
                                I = IA + 2
                            else
                                I = IA + 1
                            end
                            FDAV = FD[N, I]
                            FDQDAV = FDQ[N, I]
                            SRλ[N] = (Sλ[N] / π) * (FDAV + FDQDAV * (refl / (1.0 - (refl * S[N])))) / 1000.0
                            SRλ[N] *= AR2
                        end
                    end

                    GRλ[N] = SRλ[N] + DRλ[N]

                    if N == 1
                        SRINT[1] = 0.0u"mW / cm^2"
                        DRRINT[1] = 0.0u"mW / cm^2"
                        DRINT[1] = 0.0u"mW / cm^2"
                        GRINT[1] = 0.0u"mW / cm^2"
                    else
                        AIλ[N] = Iλ[N]
                        AIλ[N-1] = Iλ[N-1]

                        Δλ = AIλ[N] - AIλ[N-1]

                        DRINT[N] = DRINT[N-1] + (Δλ * DRλ[N-1]) + (0.5 * Δλ * (DRλ[N] - DRλ[N-1]))
                        DRRINT[N] = DRRINT[N-1] + (Δλ * DRRλ[N-1]) + (0.5 * Δλ * (DRRλ[N] - DRRλ[N-1]))
                        SRINT[N] = SRINT[N-1] + (Δλ * SRλ[N-1]) + (0.5 * Δλ * (SRλ[N] - SRλ[N-1]))
                        GRINT[N] = GRINT[N-1] + (Δλ * GRλ[N-1]) + (0.5 * Δλ * (GRλ[N] - GRλ[N-1]))
                    end
                end
                GRλs[step, :] .= GRλ
                DRRλs[step, :] .= DRRλ
                DRλs[step, :] .= DRλ
                SRλs[step, :] .= SRλ
                GRs[step] = GRINT[nmax]
                DRRs[step] = DRRINT[nmax]
                DRs[step] = DRINT[nmax]
                SRs[step] = SRINT[nmax]
            else # sunrise, sunset or long day

            end
            # Store into row `step`
            Zs[step] = Z
            ZSLs[step] = Zsl
            DOYs[step] = d
            times[step] = t
            #Zs[Zs.>90u"°"] .= 90u"°"
            step += 1
        end
        HHs[i] = HH     # save today's sunrise hour angle
        tsns[i] = tsn   # save today's time of sunrise