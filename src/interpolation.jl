function daily_cycle_sine_exponential(
    initial_temperature,
    times,
    temperatures,
    min_temperature,
    max_temperature,
    next_min_temperature,
    next_max_temperature,
    time_sunrise,
    time_sunset,
    time_max_temperature,
    daily,
    iday
)
    nhours = 24
    min_temperature = u"K"(min_temperature)
    max_temperature = u"K"(max_temperature)
    next_min_temperature = u"K"(next_min_temperature)
    next_max_temperature = u"K"(next_max_temperature)
    initial_temperature = u"K"(initial_temperature)

    temperature = min_temperature[1] # initialise temperature

    # constants
    amplitude = (max_temperature - min_temperature) / 2
    sunrise_temperature = min_temperature
    reference_time = (time_max_temperature - time_sunrise) / 2 + time_sunrise
    
    # Temperature at sunset
    sunset_angle = 360 * (time_sunset - reference_time) / (2 * (time_max_temperature - time_sunrise))
    sunset_temperature = amplitude * sin(sunset_angle / 57.29577) + min_temperature + amplitude
    
    # Nighttime exponential decay constant
    decay_rate = 3 / ((2400 - time_sunset) + time_sunrise)

    for i in 1:nhours
        j = i + 1
        time = i * 100
        if time == time_sunrise
            temperature = min_temperature
        end
        if time < time_sunrise
            hours_since_sunset = (2400 - time_sunset) + time
            if daily && iday > 1
                temperature = ((min_temperature - initial_temperature) / time_sunrise) * time + initial_temperature   
            else
                decay_exponent = hours_since_sunset * decay_rate
                temperature = (sunset_temperature - sunrise_temperature) / exp(decay_exponent) + sunrise_temperature
            end
        end
        if time > time_sunrise
            if time <= time_sunset # before or at sunset, sine wave
                angle = 360 * (time - reference_time) / (2 * (time_max_temperature - time_sunrise))
                temperature = amplitude * sin(angle / 57.29577) + min_temperature + amplitude
            else # after suset, exponential decay
                decay_exponent = ((2400 - time_sunset) - (2400 - time)) * decay_rate
                if daily
                    temperature = (sunset_temperature - next_min_temperature) / exp(decay_exponent) + next_min_temperature
                else
                    temperature = (sunset_temperature - sunrise_temperature) / exp(decay_exponent) + sunrise_temperature
                end
            end
        end

        temperature = u"°C"(temperature)
        hours = Int(time ÷ 100)
        minutes = time % 100
        minutes_since_midnight = hours * 60 + minutes

        if j <= nhours 
            times[j] = minutes_since_midnight
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

function daily_cycle_linear(daily_min, daily_max, time_sunrise, time_sunset, time_of_min, time_of_max, daily, iday, previous_value)
    nhours = 24
    initial_value = 0.0 * daily_min[1]
    time_vector = fill(0.0, nhours)
    result = fill(initial_value, nhours)
    times = fill(0.0, nhours)  
    daily_mean = (daily_max + daily_min) / 2.0

    if daily == 1 && iday > 1
        initial_value = previous_value
    end

    min_step = daily_min + 0.01 * daily_mean
    max_step = daily_max - 0.01 * daily_mean

    if time_of_min < time_of_max
        # morning min, afternoon max
        hour_index_min = Int(floor(time_of_min / 100))
        hour_index_max = Int(floor(time_of_max / 100))
        morning_slope = (daily_mean - min_step) / (100.0 - time_of_min)
        midday_slope = (min_step - max_step) / (time_of_min - time_of_max)
        evening_slope = (max_step - daily_mean) / (time_of_max - 2400.0)
    else
        # morning max, afternoon min
        hour_index_min = Int(floor(time_of_max / 100))
        hour_index_max = Int(floor(time_of_min / 100))
        morning_slope = (daily_mean - max_step) / (100.0 - time_of_max)
        midday_slope = (max_step - min_step) / (time_of_max - time_of_min)
        evening_slope = (min_step - daily_mean) / (time_of_min - 2400.0)
    end


    for i in 1:nhours
        time_vector[i] = i * 100.0 - 100.0
        time = time_vector[i]

        hours = Int(time ÷ 100)
        minutes = time % 100
        minutes_since_midnight = hours * 60 + minutes
        times[i] = minutes_since_midnight

        if i == 1
            result[i] = (daily == 1 && iday > 1) ? initial_value : daily_mean
            continue
        end

        if i < hour_index_min
            result[i] = result[1] - morning_slope * (time_vector[1] - time_vector[i])
            continue
        end

        if i == hour_index_min
            result[i] = (time_of_min < time_of_max) ? min_step : max_step
            continue
        end

        if i < hour_index_max
            if time_of_min < time_of_max
                result[i] = min_step - midday_slope * (time_of_min - time_vector[i])
            else
                result[i] = max_step - midday_slope * (time_of_max - time_vector[i])
                if result[i] > max_step
                    result[i] = max_step
                end
            end
            continue
        end

        if i == hour_index_max
            result[i] = (time_of_min < time_of_max) ? max_step : min_step
            continue
        end

        if i > hour_index_max
            if time_of_min < time_of_max
                result[i] = max_step - evening_slope * (time_of_max - time_vector[i])
            else
                result[i] = min_step - evening_slope * (time_of_min - time_vector[i])
            end
        end
    end

    # Fix any negative values
    for ihour in 1:nhours
        if result[ihour] < 0.0 * daily_min[1]
            if ihour < nhours && result[ihour + 1] > 0.0 * daily_min[1]
                # FIXME: this can index out of bounds
                result[ihour] = (result[ihour - 1] + result[ihour + 1]) / 2.0
            else
                result[ihour] = abs(result[ihour])
            end
        end
    end

    return result
end

# Convert solar times (in decimal hours) to HHMM format
_hour_to_hhmm(hour::Real) = (floor(hour) * 100.0) + (hour % 1) * 60.0

function _interpolate_variable(daily_min, daily_max, initial_value, minima_time, maxima_time, time_params, swap_times::Bool=false)
    (; half_day_length, hour_solar_noon, time_sunset, daily, iday) = time_params

    sunrise_offset = minima_time
    time_sunrise = _hour_to_hhmm(half_day_length + sunrise_offset)
    
    solar_noon_offset = maxima_time
    time_max_temperature = _hour_to_hhmm(hour_solar_noon + solar_noon_offset)

    if swap_times
        time_of_min = time_max_temperature
        time_of_max = time_sunrise
    else
        time_of_min = time_sunrise
        time_of_max = time_max_temperature
    end
    
    return daily_cycle_linear(daily_min, daily_max, time_sunrise, time_sunset, time_of_min, time_of_max, daily, iday, initial_value)
end

function _calculate_solar_times(hour_angle, hour_solar_noon)
    hour_angle_int = trunc(hour_angle)
    hour_solar_noon_minutes_fraction = (hour_angle - hour_angle_int) * 60.0
    hour_solar_noon_int = trunc(hour_solar_noon)
    hour_solar_noon_fraction = (hour_solar_noon - hour_solar_noon_int) * 60.0
    time_sunset = (hour_solar_noon_int * 100.0 + hour_solar_noon_fraction) + (hour_angle_int * 100.0 + hour_solar_noon_minutes_fraction)
    
    half_day_length = hour_solar_noon - hour_angle
    half_day_hours_int = trunc(half_day_length)
    half_day_minutes_fraction = (half_day_length - half_day_hours_int) * 60.0

    return (; time_sunset, half_day_length, half_day_hours_int, half_day_minutes_fraction, hour_solar_noon_int, hour_solar_noon_fraction)
end


function interpolate_hourly_temperature(reference_temperature_min, reference_temperature_max, minima_times, maxima_times, solar_radiation_out, daily)
    ndays = length(reference_temperature_min)
    nhours = 24
    all_hours = nhours * ndays
    air_temperatures = fill(reference_temperature_min[1], all_hours)

    for iday in 1:ndays
        initial_temperature = reference_temperature_min[iday]
        times = fill(0.0, nhours)
        temperatures = fill(0.0u"°C", nhours)

        hour_angle = solar_radiation_out.hour_angle_sunrise[iday]
        hour_solar_noon = solar_radiation_out.hour_solar_noon[iday]

        (; time_sunset, half_day_length) = _calculate_solar_times(hour_angle, hour_solar_noon)
        
        sunrise_offset = minima_times[1]
        time_sunrise = _hour_to_hhmm(half_day_length + sunrise_offset)
        
        solar_noon_offset = maxima_times[1]
        time_max_temperature = _hour_to_hhmm(hour_solar_noon + solar_noon_offset)
        
        min_temperature = reference_temperature_min[iday]
        max_temperature = reference_temperature_max[iday]
        
        next_min_temperature = (iday < ndays) ? reference_temperature_min[iday+1] : min_temperature
        next_max_temperature = (iday < ndays) ? reference_temperature_max[iday+1] : max_temperature

        daily_cycle_sine_exponential(
            initial_temperature, times, temperatures, min_temperature, max_temperature, 
            next_min_temperature, next_max_temperature, time_sunrise, time_sunset, 
            time_max_temperature, daily, iday
        )
        
        dayrange = (iday*nhours-nhours+1):(iday*nhours)
        air_temperatures[dayrange] .= temperatures
    end
    
    return air_temperatures
end

function interpolate_hourly_variable(var_min, var_max, minima_times, maxima_times, solar_radiation_out, daily, swap_times)
    ndays = length(var_min)
    nhours = 24
    all_hours = nhours * ndays
    result = fill(var_min[1], all_hours)

    for iday in 1:ndays
        initial_value = var_min[iday]
        
        hour_angle = solar_radiation_out.hour_angle_sunrise[iday]
        hour_solar_noon = solar_radiation_out.hour_solar_noon[iday]

        (; time_sunset, half_day_length) = _calculate_solar_times(hour_angle, hour_solar_noon)
        
        time_params = (; half_day_length, hour_solar_noon, time_sunset, daily, iday)
        
        interpolated_values = _interpolate_variable(var_min[iday], var_max[iday], initial_value, minima_times, maxima_times, time_params, swap_times)
        
        dayrange = (iday*nhours-nhours+1):(iday*nhours)
        result[dayrange] .= interpolated_values
    end
    
    return result
end

function hourly_from_min_max(minmax, solar_radiation_out, daily::Bool=false)
    (; reference_temperature_min, reference_temperature_max, reference_wind_min, reference_wind_max, reference_humidity_min, 
        reference_humidity_max, cloud_min, cloud_max, minima_times, maxima_times) = minmax

    air_temperatures = interpolate_hourly_temperature(reference_temperature_min, reference_temperature_max, minima_times, maxima_times, solar_radiation_out, daily)
    
    wind_speeds = interpolate_hourly_variable(reference_wind_min, reference_wind_max, minima_times[2], maxima_times[2], solar_radiation_out, daily, false)
    
    humidities = interpolate_hourly_variable(reference_humidity_min, reference_humidity_max, minima_times[3], maxima_times[3], solar_radiation_out, daily, true)
    
    cloud_covers = interpolate_hourly_variable(cloud_min, cloud_max, minima_times[4], maxima_times[4], solar_radiation_out, daily, false)

    return (; reference_temperature=air_temperatures,
            reference_wind_speed=wind_speeds,
            reference_humidity=humidities,
            cloud_cover=cloud_covers)
end
