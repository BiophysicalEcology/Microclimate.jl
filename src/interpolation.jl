function sine_exponential!(
    initial_temperature,
    times,
    temperatures,
    minimum_temperature,
    maximum_temperature,
    next_minimum_temperature,
    next_maximum_temperature,
    time_sunrise,
    time_sunset,
    time_maximum_temperature,
    daily,
    iday
)
    nhours = 24
    minimum_temperature = u"K"(minimum_temperature)
    maximum_temperature = u"K"(maximum_temperature)
    next_minimum_temperature = u"K"(next_minimum_temperature)
    next_maximum_temperature = u"K"(next_maximum_temperature)
    initial_temperature = u"K"(initial_temperature)

    temperature = minimum_temperature[1] # initialise temperature

    # constants
    amplitude = (maximum_temperature - minimum_temperature) / 2
    sunrise_temperature = minimum_temperature
    reference_time = (time_maximum_temperature - time_sunrise) / 2 + time_sunrise
    
    # Temperature at sunset
    sunset_angle = 360 * (time_sunset - reference_time) / (2 * (time_maximum_temperature - time_sunrise))
    sunset_temperature = amplitude * sin(sunset_angle / 57.29577) + minimum_temperature + amplitude
    
    # Nighttime exponential decay constant
    decay_rate = 3 / ((2400 - time_sunset) + time_sunrise)

    for i in 1:nhours
        j = i + 1
        time = i * 100
        if time == time_sunrise # sunrise
            temperature = minimum_temperature
        end
        if time < time_sunrise
            TI = (2400 - time_sunset) + time
            if daily && iday > 1
                temperature = ((minimum_temperature - initial_temperature) / time_sunrise) * time + initial_temperature   
            else
                E = TI * decay_rate
                temperature = (sunset_temperature - sunrise_temperature) / exp(E) + sunrise_temperature
            end
        end
        if time > time_sunrise
            if time <= time_sunset # before or at sunset, sine wave
                angle = 360 * (time - reference_time) / (2 * (time_maximum_temperature - time_sunrise))
                temperature = amplitude * sin(angle / 57.29577) + minimum_temperature + amplitude
            else # after suset, exponential decay
                E = ((2400 - time_sunset) - (2400 - time)) * decay_rate
                if daily
                    temperature = (sunset_temperature - next_minimum_temperature) / exp(E) + next_minimum_temperature
                else
                    temperature = (sunset_temperature - sunrise_temperature) / exp(E) + sunrise_temperature
                end
            end
        end

        temperature = u"°C"(temperature)
        ITIME = Int(floor(time / 100))
        FRMIN = time / 100 - ITIME
        ITIME *= 60
        FRMIN *= 100
        TIMEC = ITIME + FRMIN

        if j <= nhours 
            times[j] = TIMEC
            temperatures[j] = temperature
        end
    end

    # Set first time step
    times[1] = 0
    if daily == 1 && iday > 1
        temperatures[1] = initial_temperature
    else
        temperatures[1] = temperature
    end

    return nothing
end

function vsine(VMIN, VMAX, time_sunrise, time_sunset, TIMIN, TIMAX, daily, iday, IVAR)
    nhours = 24
    vinit = 0.0 * VMIN[1]
    XA = fill(0.0, nhours)
    YA = fill(vinit, nhours)
    times = fill(0.0, nhours)  
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


    for i in 1:nhours
        XA[i] = i * 100.0 - 100.0
        time = XA[i]

        ITIME = Int(floor(time / 100.0))
        FRMIN = (time / 100.0) - ITIME
        ITIME *= 60
        FRMIN *= 100.0
        TIMEC = ITIME + FRMIN
        times[i] = TIMEC

        if i == 1
            YA[i] = (daily == 1 && iday > 1) ? vinit : vave
            continue
        end

        if i < ITEST1
            YA[i] = YA[1] - slope1 * (XA[1] - XA[i])
            continue
        end

        if i == ITEST1
            YA[i] = (TIMIN < TIMAX) ? vsmIN : vsmAX
            continue
        end

        if i < ITEST2
            if TIMIN < TIMAX
                YA[i] = vsmIN - slope2 * (TIMIN - XA[i])
            else
                YA[i] = vsmAX - slope2 * (TIMAX - XA[i])
                if YA[i] > vsmAX
                    YA[i] = vsmAX
                end
            end
            continue
        end

        if i == ITEST2
            YA[i] = (TIMIN < TIMAX) ? vsmAX : vsmIN
            continue
        end

        if i > ITEST2
            if TIMIN < TIMAX
                YA[i] = vsmAX - slope3 * (TIMAX - XA[i])
            else
                YA[i] = vsmIN - slope3 * (TIMIN - XA[i])
            end
        end
    end

    # Fix any negative values
    for JCT in 1:nhours
        if YA[JCT] < 0.0 * VMIN[1]
            if JCT < nhours && YA[JCT + 1] > 0.0 * VMIN[1]
                # FIXME: this can index out of bounds
                YA[JCT] = (YA[JCT - 1] + YA[JCT + 1]) / 2.0
            else
                YA[JCT] = abs(YA[JCT])
            end
        end
    end

    return YA
end

function hourly_vars(minmax, solrad_out, daily::Bool=false)
    (; reference_temperature_min, reference_temperature_max, reference_wind_min, reference_wind_max, reference_humidity_min, 
        reference_humidity_max, cloud_min, cloud_max, minima_times, maxima_times) = minmax

    ndays = length(reference_temperature_min)
    nhours = 24
    all_hours = nhours * ndays
    temperature = fill(reference_temperature_min[1], all_hours)
    cloud_cover = fill(cloud_min[1], all_hours)
    humidity = fill(reference_humidity_min[1], all_hours)
    wind = fill(reference_wind_min[1], all_hours)

    for iday in 1:ndays
        initial_temperature = reference_temperature_min[iday] # initial air temperature for daily
        initial_wind = reference_wind_min[iday]
        initial_humidity = reference_humidity_max[iday]
        initial_cloud = cloud_min[iday]
        times = fill(0.0, nhours)
        temperature = fill(0.0u"°C", nhours)
        wind = fill(initial_wind, nhours)
        humidity = fill(initial_humidity, nhours)
        clouds = fill(initial_cloud, nhours)

        HH = solrad_out.hour_angle_sunrise[iday]
        tsn = solrad_out.hour_solar_noon[iday]

        #     Air temperature calculations
        # sunset in military time
        HHINT = trunc(HH)
        FRACT = (HH - HHINT) * 60.0
        HSINT = trunc(tsn)
        FRACTS = (tsn - HSINT) * 60.0
        time_sunset = (HSINT * 100.0 + FRACTS) + (HHINT * 100.0 + FRACT)
        # sunrise in military time
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
        time_sunrise = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        # time of air temperature maximum
        TSNHR = maxima_times[1]# + TIMCOR (TIMCOR never used in Fortran version - uninitialised?)
        time_maximum_temperature = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        #     SETTING minimum_temperature, maximum_temperature FROM ARRAYS OBTAINED FROM SUB. IOSOLR
        minimum_temperature = reference_temperature_min[iday]
        maximum_temperature = reference_temperature_max[iday]
        if iday < ndays & ndays > 1
            next_minimum_temperature = reference_temperature_min[iday+1]
            next_maximum_temperature = reference_temperature_max[iday+1]
        else
            next_minimum_temperature = minimum_temperature
            next_maximum_temperature = maximum_temperature
        end

        # setting time of minimum and maximum (hours before sunrise or after solar noon)     #     AFTER SOLAR NOON)
        TIMIN = time_sunrise
        TIMAX = time_maximum_temperature
        sine_exponential!(initial_temperature, times, temperature, minimum_temperature, maximum_temperature, next_minimum_temperature, next_maximum_temperature, time_sunrise, time_sunset, time_maximum_temperature, daily, iday)

        # wind speed
        VMIN = wind_min[iday]
        VMAX = wind_max[iday]
        #      SETTING MAX & MIN TIMES RELATIVE TO SUNRISE & SOLAR NOON
        #     TIME OF MINIMUM
        TSRHR = minima_times[2]
        time_sunrise = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        #      TIME OF MAXIMUM at sunrise for relative humidities
        TSNHR = maxima_times[2] #+ TIMCOR
        time_maximum_temperature = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        TIMIN = time_sunrise
        TIMAX = time_maximum_temperature
        IVAR = initial_wind
        wind = vsine(VMIN, VMAX, time_sunrise, time_sunset, TIMIN, TIMAX, daily, iday, IVAR)

        # relative humidities
        VMIN = reference_humidity_min[iday]
        VMAX = reference_humidity_max[iday]
        #      SETTING MAX & MIN TIMES RELATIVE TO SUNRISE & SOLAR NOON
        #     TIME OF MINIMUM
        TSRHR = minima_times[3]
        time_sunrise = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        #      TIME OF MAXIMUM at sunrise for relative humidities
        TSNHR = maxima_times[3] #+ TIMCOR
        time_maximum_temperature = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        TIMIN = time_maximum_temperature
        TIMAX = time_sunrise
        IVAR = initial_humidity
        humidity = vsine(VMIN, VMAX, time_sunrise, time_sunset, TIMIN, TIMAX, daily, iday, IVAR)

        # cloud_covers cover
        VMIN = cloud_min[iday]
        VMAX = cloud_max[iday]
        #      SETTING MAX & MIN TIMES RELATIVE TO SUNRISE & SOLAR NOON
        #     TIME OF MINIMUM
        TSRHR = minima_times[4]
        time_sunrise = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        #      TIME OF MAXIMUM at sunrise for relative humidities
        TSNHR = maxima_times[4] #+ TIMCOR
        time_maximum_temperature = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        TIMIN = time_sunrise
        TIMAX = time_maximum_temperature
        IVAR = initial_cloud
        clouds = vsine(VMIN, VMAX, time_sunrise, time_sunset, TIMIN, TIMAX, daily, iday, IVAR)

        dayrange = (iday*nhours-nhours+1):(iday*nhours)
        reference_temperature[dayrange] .= temperature
        reference_wind_speed[dayrange] .= wind
        reference_humidity[dayrange] .= humidity
        cloud_cover[dayrange] .= clouds
    end

    return (; reference_temperature, reference_wind_speed, reference_humidity, cloud_cover)
end

# TODO this does just cloud_cover but should generalise first version better down the track
function hourly_vars(
    cloud_min::Vector,
    cloud_max::Vector,
    solrad_out::Any,
    minima_times::Vector,
    maxima_times::Vector,
    daily::Bool=false
)

    ndays = length(cloud_min)
    nhours = 24
    cloud_cover = fill(cloud_min[1], nhours * ndays)

    for iday in 1:ndays
        initial_cloud = cloud_min[iday]
        times = fill(0.0, 25)
        cloud = fill(initial_cloud,  nhours)

        HH = solrad_out.hour_angle_sunrise[iday]
        tsn = solrad_out.hour_solar_noon[iday]

        #     Air temperature calculations
        #     SUNSET IN MILITARY TIME
        HHINT = trunc(HH)
        FRACT = (HH - HHINT) * 60.0
        HSINT = trunc(tsn)
        FRACTS = (tsn - HSINT) * 60.0
        time_sunset = (HSINT * 100.0 + FRACTS) + (HHINT * 100.0 + FRACT)
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
        time_sunrise = (HRINT * 100.0 + FRACTR) + (TSRHR * 100.0)
        #      TIME OF MAXIMUM at sunrise for relative humidities
        TSNHR = maxima_times[4] #+ TIMCOR
        time_maximum_temperature = (HSINT * 100.0 + FRACTS) + (TSNHR * 100.0)
        TIMIN = time_sunrise
        TIMAX = time_maximum_temperature
        IVAR = initial_cloud
        cloud = vsine(VMIN, VMAX, time_sunrise, time_sunset, TIMIN, TIMAX, daily, iday, IVAR)

        cloud_cover[(iday*nhours-nhours+1):(iday*nhours)] = cloud[1:nhours]
    end
    return cloud_cover
end
