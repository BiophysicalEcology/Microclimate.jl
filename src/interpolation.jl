function sinec!(
    ITAIR, TIMARY, TAIRRY, TMIN, TMAX, TMIN2, TMAX2, TIMSR, TIMSS, TIMTMX, daily, iday
)
    TMIN = u"K"(TMIN)
    TMAX = u"K"(TMAX)
    TMIN2 = u"K"(TMIN2)
    TMAX2 = u"K"(TMAX2)
    ITAIR = u"K"(ITAIR)
    T=TMIN[1] # initialise T
    A = (TMAX - TMIN) / 2
    TSR = TMIN
    TREF = (TIMTMX - TIMSR) / 2 + TIMSR
    SS = 360.0 * (TIMSS - TREF) / (2.0 * (TIMTMX - TIMSR))
    SY = SS / 57.29577
    ZS = sin(SY)
    TSS = A * ZS + TMIN + A
    TAU = 3.0 / ((2400.0 - TIMSS) + TIMSR)

    for I in 1:24
        J = I + 1
        TIME = I * 100.0
        if TIME == TIMSR # sunrise
            T = TMIN
        end
        if TIME < TIMSR
            TI = (2400.0 - TIMSS) + TIME
            if daily && iday > 1
                T = ((TMIN - ITAIR) / TIMSR) * TIME + ITAIR   
            else
                E = TI * TAU
                T = (TSS - TSR) / exp(E) + TSR
            end
        end
        if TIME > TIMSR
            if TIME <= TIMSS # before or at sunset
                X = 360.0 * (TIME - TREF) / (2.0 * (TIMTMX - TIMSR))
                Y = X / 57.29577
                Z = sin(Y)
                T = A * Z + TMIN + A
            else # after suset
                TI = (2400.0 - TIMSS) - (2400. - TIME)
                E = TI * TAU
                if daily
                    T = (TSS - TMIN2) / exp(E) + TMIN2
                else
                    T = (TSS - TSR) / exp(E) + TSR
                end
            end
        end

        T = u"°C"(T)
        ITIME = Int(floor(TIME / 100.0))
        FRMIN = TIME / 100.0 - ITIME
        ITIME *= 60
        FRMIN *= 100.0
        TIMEC = ITIME + FRMIN

        TIMARY[J] = TIMEC
        TAIRRY[J] = T
    end

    # Set first time step
    TIMARY[1] = 0.0
    if daily == 1 && iday > 1
        TAIRRY[1] = ITAIR
    else
        TAIRRY[1] = T
    end
    ITAIR = TAIRRY[25]
end

function vsine(VMIN, VMAX, TIMSR, TIMSS, TIMIN, TIMAX, daily, iday, IVAR)
    vinit = 0.0 * VMIN[1]
    XA = fill(0.0, 25)
    YA = fill(vinit, 25)
    TIMARY = fill(0.0, 25)  
    vave = (VMAX + VMIN) / 2.0

    if daily == 1 && iday > 1
        vinit = IVAR
    end

    vsmIN = VMIN + 0.01 * vave
    vsmAX = VMAX - 0.01 * vave

    if TIMIN < TIMAX
        # morning min, afternoon max
        ITEST1 = Int(floor(TIMIN / 100))
        ITEST2 = Int(floor(TIMAX / 100))
        slope1 = (vave - vsmIN) / (100.0 - TIMIN)
        slope2 = (vsmIN - vsmAX) / (TIMIN - TIMAX)
        slope3 = (vsmAX - vave) / (TIMAX - 2400.0)
    else
        # morning max, afternoon min
        ITEST1 = Int(floor(TIMAX / 100))
        ITEST2 = Int(floor(TIMIN / 100))
        slope1 = (vave - vsmAX) / (100.0 - TIMAX)
        slope2 = (vsmAX - vsmIN) / (TIMAX - TIMIN)
        slope3 = (vsmIN - vave) / (TIMIN - 2400.0)
    end


    for I in 1:25
        XA[I] = I * 100.0 - 100.0
        TIME = XA[I]

        ITIME = Int(floor(TIME / 100.0))
        FRMIN = (TIME / 100.0) - ITIME
        ITIME *= 60
        FRMIN *= 100.0
        TIMEC = ITIME + FRMIN
        TIMARY[I] = TIMEC

        if I == 1
            YA[I] = (daily == 1 && iday > 1) ? vinit : vave
            continue
        end

        if I < ITEST1
            YA[I] = YA[1] - slope1 * (XA[1] - XA[I])
            continue
        end

        if I == ITEST1
            YA[I] = (TIMIN < TIMAX) ? vsmIN : vsmAX
            continue
        end

        if I < ITEST2
            if TIMIN < TIMAX
                YA[I] = vsmIN - slope2 * (TIMIN - XA[I])
            else
                YA[I] = vsmAX - slope2 * (TIMAX - XA[I])
                if YA[I] > vsmAX
                    YA[I] = vsmAX
                end
            end
            continue
        end

        if I == ITEST2
            YA[I] = (TIMIN < TIMAX) ? vsmAX : vsmIN
            continue
        end

        if I > ITEST2
            if TIMIN < TIMAX
                YA[I] = vsmAX - slope3 * (TIMAX - XA[I])
            else
                YA[I] = vsmIN - slope3 * (TIMIN - XA[I])
            end
        end
    end

    # Fix any negative values
    for JCT in 1:25
        if YA[JCT] < 0.0 * VMIN[1]
            if JCT < 25 && YA[JCT + 1] > 0.0 * VMIN[1]
                YA[JCT] = (YA[JCT - 1] + YA[JCT + 1]) / 2.0
            else
                YA[JCT] = abs(YA[JCT])
            end
        end
    end

    return YA
end

# this version does all variables
function hourly_vars(
    air_temperature_min::Vector,
    air_temperature_max::Vector,
    wind_min::Vector,
    wind_max::Vector,
    humidity_min::Vector,
    humidity_max::Vector,
    cloud_min::Vector,
    cloud_max::Vector,
    solrad_out::Any,
    minima_times::Vector=[0, 0, 1, 1],
    maxima_times::Vector=[1, 1, 0, 0],
    daily::Bool = false)

    ndays = length(air_temperature_min)
    air_temperatures = fill(air_temperature_min[1], 25 * ndays)
    cloud_covers = fill(cloud_min[1], 25 * ndays)
    humidites = fill(humidity_min[1], 25 * ndays)
    wind_speeds = fill(wind_min[1], 25 * ndays)

    for iday in 1:ndays
        ITAIR = air_temperature_min[iday] # initial air temperature for daily
        IWN = wind_min[iday]
        IRH = humidity_max[iday]
        ICLD = cloud_min[iday]
        TIMARY = fill(0.0, 25)
        TAIRRY = fill(0.0u"°C", 25)
        WNARRY = fill(IWN, 25)
        RHARRY = fill(IRH, 25)
        CLDARRY = fill(ICLD, 25)

        HH = solrad_out.hour_angle_sunrise[iday]
        tsn = solrad_out.hour_solar_noon[iday]

        #     Air temperature calculations
        #     SUNSET IN MILITARY TIME
        HHINT = trunc(HH)
        FRACT = (HH - HHINT) * 60.0
        HSINT = trunc(tsn)
        FRACTS = (tsn - HSINT) * 60.0
        TIMSS = (HSINT * 100.0 + FRACTS) + (HHINT * 100.0 + FRACT)
        #     SUNRISE IN MILITARY TIME
        DELT = tsn - HH
        HRINT = trunc(DELT)
        FRACTR = (DELT - HRINT) * 60.0
        #     TIME OF AIR TEMPERATURE MINIMUM (NOTE: A KLUGE OF 200, I.E. 2
        #     HOURS IS BEING USED TO MAKE SINEC WORK WHEN MORNING MINIMUM
        #     IS NOT AT TRUE SUNRISE, THE ALGORITHM IS 2 HOURS OFF FOR SOME
        #     UNKNOWN (6/7/89) REASON)
        if minima_times[1] > 0
            TSRHR = minima_times[1] + 2.0
        else
            TSRHR = minima_times[1]
        end
        TIMSR = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        #     TIME OF AIR TEMPERATURE MAXIMUM
        TSNHR = maxima_times[1]# + TIMCOR (TIMCOR never used in Fortran version - uninitialised?)
        TIMTMX = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        #     SETTING TMIN, TMAX FROM ARRAYS OBTAINED FROM SUB. IOSOLR
        TMIN = air_temperature_min[iday]
        TMAX = air_temperature_max[iday]
        if iday < ndays & ndays > 1
            TMIN2 = air_temperature_min[iday+1]
            TMAX2 = air_temperature_max[iday+1]
        else
            TMIN2 = TMIN
            TMAX2 = TMAX
        end

        #     SETTING TIME OF MIN & MAX (HOURS BEFORE SUNRISE OR
        #     AFTER SOLAR NOON)
        TIMIN = TIMSR
        TIMAX = TIMTMX
        sinec!(ITAIR, TIMARY, TAIRRY, TMIN, TMAX, TMIN2, TMAX2, TIMSR, TIMSS, TIMTMX, daily, iday)

        # wind speed
        VMIN = wind_min[iday]
        VMAX = wind_max[iday]
        #      SETTING MAX & MIN TIMES RELATIVE TO SUNRISE & SOLAR NOON
        #     TIME OF MINIMUM
        TSRHR = minima_times[2]
        TIMSR = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        #      TIME OF MAXIMUM at sunrise for relative humidites
        TSNHR = maxima_times[2] #+ TIMCOR
        TIMTMX = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        TIMIN = TIMSR
        TIMAX = TIMTMX
        IVAR = IWN
        WNARRY = vsine(VMIN, VMAX, TIMSR, TIMSS, TIMIN, TIMAX, daily, iday, IVAR)

        # relative humidites
        VMIN = humidity_min[iday]
        VMAX = humidity_max[iday]
        #      SETTING MAX & MIN TIMES RELATIVE TO SUNRISE & SOLAR NOON
        #     TIME OF MINIMUM
        TSRHR = minima_times[3]
        TIMSR = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        #      TIME OF MAXIMUM at sunrise for relative humidites
        TSNHR = maxima_times[3] #+ TIMCOR
        TIMTMX = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        TIMIN = TIMTMX
        TIMAX = TIMSR
        IVAR = IRH
        RHARRY = vsine(VMIN, VMAX, TIMSR, TIMSS, TIMIN, TIMAX, daily, iday, IVAR)

        # cloud_covers cover
        VMIN = cloud_min[iday]
        VMAX = cloud_max[iday]
        #      SETTING MAX & MIN TIMES RELATIVE TO SUNRISE & SOLAR NOON
        #     TIME OF MINIMUM
        TSRHR = minima_times[4]
        TIMSR = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        #      TIME OF MAXIMUM at sunrise for relative humidites
        TSNHR = maxima_times[4] #+ TIMCOR
        TIMTMX = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        TIMIN = TIMSR
        TIMAX = TIMTMX
        IVAR = ICLD
        CLDARRY = vsine(VMIN, VMAX, TIMSR, TIMSS, TIMIN, TIMAX, daily, iday, IVAR)

        air_temperatures[(iday*25-24):(iday*25)] = TAIRRY[1:25]
        wind_speeds[(iday*25-24):(iday*25)] = WNARRY[1:25]
        humidites[(iday*25-24):(iday*25)] = RHARRY[1:25]
        cloud_covers[(iday*25-24):(iday*25)] = CLDARRY[1:25]
    end
    return air_temperatures, wind_speeds, humidites, cloud_covers
end

# TODO this does just cloud_covers but should generalise first version better down the track
function hourly_vars(
    cloud_min::Vector,
    cloud_max::Vector,
    solrad_out::Any,
    minima_times::Vector=[0, 0, 1, 1],
    maxima_times::Vector=[1, 1, 0, 0],
    daily::Bool=false)

    ndays = length(cloud_min)
    cloud_covers = fill(cloud_min[1], 25 * ndays)

    for iday in 1:ndays
        ICLD = cloud_min[iday]
        TIMARY = fill(0.0, 25)
        CLDARRY = fill(ICLD, 25)

        HH = solrad_out.hour_angle_sunrise[iday]
        tsn = solrad_out.hour_solar_noon[iday]

        #     Air temperature calculations
        #     SUNSET IN MILITARY TIME
        HHINT = trunc(HH)
        FRACT = (HH - HHINT) * 60.0
        HSINT = trunc(tsn)
        FRACTS = (tsn - HSINT) * 60.0
        TIMSS = (HSINT * 100.0 + FRACTS) + (HHINT * 100.0 + FRACT)
        #     SUNRISE IN MILITARY TIME
        DELT = tsn - HH
        HRINT = trunc(DELT)
        FRACTR = (DELT - HRINT) * 60.0

        # cloud_covers cover
        VMIN = cloud_min[iday]
        VMAX = cloud_max[iday]
        #      SETTING MAX & MIN TIMES RELATIVE TO SUNRISE & SOLAR NOON
        #     TIME OF MINIMUM
        TSRHR = minima_times[4]
        TIMSR = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        #      TIME OF MAXIMUM at sunrise for relative humidites
        TSNHR = maxima_times[4] #+ TIMCOR
        TIMTMX = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        TIMIN = TIMSR
        TIMAX = TIMTMX
        IVAR = ICLD
        CLDARRY = vsine(VMIN, VMAX, TIMSR, TIMSS, TIMIN, TIMAX, daily, iday, IVAR)

        cloud_covers[(iday*25-24):(iday*25)] = CLDARRY[1:25]
    end
    return cloud_covers
end