function get_profile(;
    z0=0.004u"m",
    zh=0.0u"m",
    d0=0.0u"m",
    κ=0.4,
    heights::Vector{typeof(1.0u"m")}=[0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.2] .* u"m",
    reference_temperature=27.77818u"°C",
    reference_wind_speed=2.749575u"m/s",
    relative_humidity=49.0415,
    surface_temperature=48.58942u"°C",
    maximum_surface_temperature=40.0u"°C",
    zenith_angle=21.50564u"°",
    elevation=0.0u"m",
    P_atmos=atmospheric_pressure(elevation),
)
    reference_height = last(heights)
    if minimum(heights) < z0
        error("ERROR: the minimum height is not greater than the roughness height (z0).")
    end

    T_ref_height = u"K"(reference_temperature)
    T_surface = u"K"(surface_temperature)

    # Units: m to cm
    z = u"cm"(reference_height)
    z0 = u"cm"(z0)
    zh_cm = u"cm"(zh)
    d0_cm = u"cm"(d0)
    refence_wind_cm_min = u"cm/minute"(reference_wind_speed)
    # define air heights
    N_heights = length(heights)
    height_array = u"cm".(reverse(heights))
    wind_speeds = zeros(Float64, N_heights) .* 1u"cm/minute" # output wind speeds
    air_temperatures = Vector{typeof(0.0u"K")}(undef, N_heights) # output temperatures, need to do this otherwise get InexactError
    humidities = zeros(Float64, N_heights) # output relative humidities
    wind_speeds[1] = refence_wind_cm_min
    air_temperatures[1] = T_ref_height

    # compute rcptkg (was a constant in original Fortran version)
    dry_air_out = dry_air_properties(u"K"(reference_temperature); elevation, P_atmos)
    # TODO: why do we not pass P_atmos here
    wet_air_out = wet_air_properties(u"K"(reference_temperature); rh=relative_humidity)
    ρ = dry_air_out.ρ_air
    c_p = wet_air_out.c_p
    # TODO make this work with SI units
    rcptkg = u"cal*minute^2/cm^4"(ρ * c_p * T_ref_height / (κ * g_n))
    #rcptkg = 6.003e-8u"cal*minute^2/cm^4"
    gam = 16.0
    z_ratio = z / z0 + 1.0
    dum = log(z_ratio)
    u_star = κ * refence_wind_cm_min / dum
    ΔT = T_ref_height - T_surface
    T_mean = (T_surface + T_ref_height) / 2
    # TODO call rho_cp method specific to elevation and RH in final version but do it this way for NicheMapR comparison
    ρ_cp = rho_cp(T_mean)#, elevation, relative_humidity)
    amol = -30.0u"cm"
    sublayer_stanton_number = 0.62 / (ustrip(u"cm", z0) * ustrip(u"cm/minute", u_star) / 12)^(9//20)
    bulk_stanton_number = 0.64 / dum
    Q_convection = ρ_cp * ΔT * u_star * bulk_stanton_number / (1.0 + bulk_stanton_number / sublayer_stanton_number)
    if zh > 0.0u"m"
        for i in 2:N_heights
            if T_ref_height ≥ T_surface || T_surface ≤ u"K"(maximum_surface_temperature) || zenith_angle ≥ 90°
                wind_speeds[i] = (u_star / κ) * log(height_array[i] / z0 + 1)
            else
                x1 = phi(height_array[i], gam, amol)
                y1 = psi1(x1)
                adum = height_array[i] / z0 - y1
                wind_speeds[i] = (u_star / κ) * log(adum)
            end
            A = (T_ref_height - T_surface) / (1 - log((z - d0_cm) / zh_cm))
            T0 = T_ref_height + A * log((z - d0_cm) / zh_cm)
            air_temperatures[i] = T0 - A * log((height_array[i] - d0_cm) / zh_cm)
        end
    else
        if T_ref_height ≥ T_surface || T_surface ≤ u"K"(maximum_surface_temperature) || zenith_angle ≥ 90°
            for i in 2:N_heights
                wind_speeds[i] = (u_star / κ) * log(height_array[i] / z0 + 1.0)
                T_z0 = (T_ref_height * bulk_stanton_number + T_surface * sublayer_stanton_number) / (bulk_stanton_number + sublayer_stanton_number)
                air_temperatures[i] = T_z0 + (T_ref_height - T_z0) * log(height_array[i] / z0 + 1.0) / dum
            end
        else
            for i in 2:N_heights
                x1 = phi(height_array[i], gam, amol)
                y1 = psi1(x1)
                yy2 = psi2(x1)
                x = phi(z, gam, amol)
                yy = psi2(x)
                adum = height_array[i] / z0 - y1
                wind_speeds[i] = (u_star / κ) * log(adum)
                Obukhov_out = get_Obukhov(T_ref_height, T_surface, refence_wind_cm_min, height_array[i], z0, rcptkg, κ)
                T_z0 = (T_ref_height * Obukhov_out.bulk_stanton_number + T_surface * Obukhov_out.sublayer_stanton_number) / (Obukhov_out.bulk_stanton_number + Obukhov_out.sublayer_stanton_number)
                Q_convection = ρ_cp * ΔT * u_star * Obukhov_out.bulk_stanton_number / (1.0 + Obukhov_out.bulk_stanton_number / Obukhov_out.sublayer_stanton_number)
                air_temperatures[i] = T_z0 + (T_ref_height - T_z0) * log(height_array[i] / z0 - yy2) / log(z / z0 - yy)
            end
        end
    end
    wind_speeds = reverse(wind_speeds)
    air_temperatures = reverse(air_temperatures)
    
    e = wet_air_properties(T_ref_height; rh = relative_humidity).P_vap
    humidities .= clamp.(e ./ vapour_pressure.(air_temperatures) .* 100.0, 0.0, 100.0)

    return (;
        heights,
        wind_speeds,
        air_temperatures,
        humidities,
        qconv=u"W/m^2"(Q_convection),
        ustar=u"m/s"(u_star)
    )
end

function rho_cp(T_mean)
    return u"(cal*g)/(g*cm^3*K)" * (0.08472 / ustrip(u"K", T_mean))
end
function rho_cp(T_mean, elevation, relative_humidity, P_atmos)
    dry_air_out = dry_air_properties(u"K"(T_mean); elevation, P_atmos)
    # TODO: why do we not pass P_atmos here
    wet_air_out = wet_air_properties(u"K"(T_mean); rh=relative_humidity)
    ρ = dry_air_out.ρ_air
    c_p = wet_air_out.c_p
    return u"(cal*g)/(g*cm^3*K)"(ρ * c_p)
end

function phi(z, gam, amol)
    # TODO ustrip to what
    return (1.0 - min(1.0, gam * ustrip(z / amol)))^(1//4)
end

function psi1(x)
    return 2.0 * log((1.0 + x) / 2.0) + log((1.0 + x^2) / 2.0) - 2.0 * atan(x) + π / 2.0
end

function psi2(x)
    return 2.0 * log((1 + x^2.0) / 2.0)
end

function get_Obukhov(TA, TS, refence_wind_cm_min, z, z0, rcptkg, κ)
    amol = -30.0u"cm" # initial Monin-Obukhov length cm
    gam = 16.0 # -
    #RCPTKG = 6.003e-8u"cal/minute/cm/K" #CAL-MIN-CM-C
    z = u"cm"(z)
    z0 = u"cm"(z0)
    z_ratio = z / z0 + 1.0
    dum = log(z_ratio)
    TA = u"K"(TA)
    TS = u"K"(TS)
    ΔT = TA - TS
    T_mean = (TA + TS) / 2.0
    ρ_cp = rho_cp(T_mean)
    del = 1.0
    count = 0
    u_star = (κ * refence_wind_cm_min / dum)u"cm/minute"
    Q_convection = 0.0u"cal/minute/cm^2"
    STO = 0.0
    bulk_stanton_number = 0.0
    sublayer_stanton_number = 0.0

    while del > 1e-2 && count < 500
        count += 1
        x = phi(z, gam, amol)
        y = psi1(x)
        u_star = κ * refence_wind_cm_min / (log(z / z0) - y)

        if amol > 0.0u"cm"
            sublayer_stanton_number = 0.62 / (ustrip(u"cm", z0) * ustrip(u"cm/minute", u_star) / 12)^(9//20)
            bulk_stanton_number = 0.64 / dum
            Q_convection = ρ_cp * ΔT * u_star * bulk_stanton_number / (1.0 + bulk_stanton_number / sublayer_stanton_number)
        else
            sublayer_stanton_number = 0.62 / (ustrip(u"cm", z0) * ustrip(u"cm/minute", u_star) / 12)^(9//20)
            bulk_stanton_number = (0.64 / dum) * (1 - 0.1 * z / amol)
            STO = bulk_stanton_number / (1 + bulk_stanton_number / sublayer_stanton_number)
            Q_convection = ρ_cp * ΔT * u_star * STO
        end

        AMOLN = rcptkg * u_star^3 / Q_convection
        del = abs((AMOLN - amol) / amol)
        amol = AMOLN
    end

    return (; amol=u"m"(amol), sublayer_stanton_number, STO, bulk_stanton_number, u_star, Q_convection)
end

