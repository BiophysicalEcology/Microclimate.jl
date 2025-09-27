"""
    get_profile(; kwargs...)

Compute vertical profiles of wind speed, air temperature, and relative humidity 
in the atmospheric surface layer, using Monin–Obukhov similarity theory (MOST).

This function reproduces the `get_profile` routine from **NicheMapR**, ported to Julia.  
It calculates the microclimate profiles above the ground (or canopy) at specified heights,
based on reference conditions and surface parameters.

# Keyword Arguments
- `z0::Quantity=0.004u"m"`: roughness length (surface aerodynamic roughness).
- `zh::Quantity=0.0u"m"`: heat transfer roughness height
- `d0::Quantity=0.0u"m"`: zero plane displacement correction factor.
- `κ::Float64=0.4`: von Kármán constant.
- `heights::Vector{Quantity}`: Measurement heights above the surface (default: 0.01–1.2 m).
- `reference_temperature::Quantity=27.78u"°C"`: Air temperature at the reference height.
- `reference_wind_speed::Quantity=2.75u"m/s"`: Wind speed at the reference height.
- `relative_humidity::Float64=49.0`: Relative humidity at the reference height (%).
- `surface_temperature::Quantity=48.59u"°C"`: Soil or surface temperature.
- `maximum_surface_temperature::Quantity=40.0u"°C"`: Maximum allowed surface temperature.
- `zenith_angle::Quantity=21.5u"°"`: Solar zenith angle.
- `elevation::Quantity=0.0u"m"`: Elevation above sea level.

# Returns
Named tuple with fields:
- `heights`: Input heights (reversed to increasing order).
- `wind_speeds`: Wind speed profile at each height (`cm/min` internally, returned in SI units).
- `air_temperatures`: Air temperature profile at each height (`K`).
- `humidities`: Relative humidity profile (%) at each height.
- `qconv`: Convective heat flux (`W/m²`).
- `ustar`: Friction velocity (`m/s`).

# Notes
- Stability corrections use the **Businger–Dyer** formulations for unstable conditions.
- The Monin–Obukhov length is estimated iteratively through `calc_Obukhov_length`.
- Two broad options for aerodynamic roughness calculations are available: Campbell & Norman's (1998) approach
that handles canopy displacement, invoked if `zh > 0` and otherwise  
- When `zh > 0`, canopy displacement is considered in the profile calculation.
- zh and d0 for Campbell and Norman air temperature/wind speed profile (0.6 * canopy height in m if unknown
| Condition                   | Wind profile                   | Temperature profile                          |
| --------------------------- | ------------------------------ | -------------------------------------------- |
| `zh > 0` + neutral/hot      | log-law                        | log between `z` and `zh`                     |
| `zh > 0` + unstable/stable  | log-law with `calc_ψ_m` correction | log with displacement/`zh`                   |
| `zh == 0` + neutral/hot     | log-law                        | weighted by bulk/sublayer Stanton numbers    |
| `zh == 0` + unstable/stable | log-law with `calc_ψ_m` correction | full Monin–Obukhov profile via `calc_Obukhov_length` |

- Relative humidity profiles are estimated from vapor pressure at each height.

# References
- Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971).
  Flux–profile relationships in the atmospheric surface layer.
  *Journal of the Atmospheric Sciences*, 28(2), 181–189.
- Dyer, A. J. (1974). A review of flux–profile relationships.
  *Boundary-Layer Meteorology*, 7(3), 363–372.
- Kearney, M. R., et al. (2020). NicheMapR: an R package for microclimate and 
  biophysical modeling. *Ecography*, 43, 1–14.

# Example
```julia
profile = get_profile(
    reference_temperature = 25u"°C",
    reference_wind_speed = 2.0u"m/s",
    relative_humidity = 60.0,
    surface_temperature = 35u"°C",
    zenith_angle = 45u"°"
)

profile.air_temperatures  # vertical profile of air temperatures
profile.wind_speeds       # vertical profile of wind speeds
"""
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
    elevation=0.0u"m", # to be used if getting ρ_cp based on humidity and pressure
)
    # z0=0.004u"m"
    # zh=0.004u"m"
    # d0=0.12u"m"
    # κ=0.4
    # heights=[0.33, 2] .* u"m"
    # reference_temperature=27.77818u"°C"
    # reference_wind_speed=2.749575u"m/s"
    # relative_humidity=49.0415
    # surface_temperature=48.58942u"°C"
    # maximum_surface_temperature=40.0u"°C"
    # zenith_angle=21.50564u"°"
    # elevation=0.0u"m"
    
    reference_height = last(heights)
    if minimum(heights) < z0
        error("ERROR: the minimum height is not greater than the roughness height (z0).")
    end

    T_ref_height = u"K"(reference_temperature)
    T_surface = u"K"(surface_temperature)

    # Units: m to cm
    z = u"cm"(reference_height)
    z0 = u"cm"(z0) # roughness height
    zh_cm = u"cm"(zh)
    d0_cm = u"cm"(d0)
    Ū_ref_height = u"cm/minute"(reference_wind_speed)
    # define air heights
    N_heights = length(heights)
    height_array = u"cm".(reverse(heights))
    wind_speeds = zeros(Float64, N_heights) .* 1u"cm/minute" # output wind speeds
    air_temperatures = Vector{typeof(0.0u"K")}(undef, N_heights) # output temperatures, need to do this otherwise get InexactError
    humidities = zeros(Float64, N_heights) # output relative humidities
    wind_speeds[1] = Ū_ref_height
    air_temperatures[1] = T_ref_height

    # compute ρcpTκg (was a constant in original Fortran version)
    #dry_air_out = dry_air_properties(u"K"(reference_temperature), elevation=elevation)
    #wet_air_out = wet_air_properties(u"K"(reference_temperature), rh = relative_humidity)
    #ρ = dry_air_out.ρ_air
    #c_p = wet_air_out.c_p
    # TODO make this work with SI units
    #ρcpTκg = u"cal*minute^2/cm^4"(ρ * c_p * T_ref_height / (κ * g_n))
    ρcpTκg = 6.003e-8u"cal*minute^2/cm^4"
    γ = 16.0 # coefficient from Dyer and Hicks for Φ_m (momentum), TODO make it available as a user param?
    log_z_ratio = log(z / z0 + 1.0) # save compute time by precalculating, + 1.0 to avoid singularities
    u_star = κ * Ū_ref_height / log_z_ratio
    ΔT = T_ref_height - T_surface
    T_mean = (T_surface + T_ref_height) / 2
    # TODO call calc_ρ_cp method specific to elevation and RH in final version but do it this way for NicheMapR comparison
    ρ_cp = calc_ρ_cp(T_mean)#, elevation, relative_humidity)
    L_Obukhov = -30.0u"cm" # initialise Obukhov length
    sublayer_stanton_number = 0.62 / (ustrip(u"cm", z0) * ustrip(u"cm/minute", u_star) / 12)^(9//20)
    bulk_stanton_number = 0.64 / log_z_ratio
    Q_convection = ρ_cp * ΔT * u_star * bulk_stanton_number / (1.0 + bulk_stanton_number / sublayer_stanton_number)
    if zh > 0.0u"m" # Campbell & Norman canopy displacement approach
        for i in 2:N_heights
            A = (T_ref_height - T_surface) / (1 - log((z - d0_cm) / zh_cm))
            T0 = T_ref_height + A * log((z - d0_cm) / zh_cm)
            air_temperatures[i] = T0 - A * log((height_array[i] - d0_cm) / zh_cm)
        end
    end
    if T_ref_height ≥ T_surface || T_surface ≤ u"K"(maximum_surface_temperature) || zenith_angle ≥ 90°
        for i in 2:N_heights
            wind_speeds[i] = calc_wind(height_array[i], z0, κ, u_star, 1.0)
            T_z0 = (T_ref_height * bulk_stanton_number + T_surface * sublayer_stanton_number) / (bulk_stanton_number + sublayer_stanton_number)
            if zh <= 0.0u"m"
                air_temperatures[i] = T_z0 + (T_ref_height - T_z0) * log(height_array[i] / z0 + 1.0) / log_z_ratio
            end
        end
    else
        for i in 2:N_heights
            Obukhov_out = calc_Obukhov_length(T_ref_height, T_surface, Ū_ref_height, height_array[i], z0, ρcpTκg, κ, log_z_ratio, ΔT, ρ_cp)
            L_Obukhov = Obukhov_out.L_Obukhov
            T_z0 = Obukhov_out.T_z0
            Q_convection = Obukhov_out.Q_convection
            φ_m1 = calc_φ_m(height_array[i], γ, L_Obukhov)
            ψ_m1 = calc_ψ_m(φ_m1)
            ψ_h2 = calc_ψ_h(ψ_m1)
            φ_m = calc_φ_m(z, γ, L_Obukhov)
            ψ_h = calc_ψ_h(φ_m)
            wind_speeds[i] = calc_wind(height_array[i], z0, κ, u_star, -ψ_m1)
            if zh <= 0.0u"m"
                air_temperatures[i] = T_z0 + (T_ref_height - T_z0) * log(height_array[i] / z0 - ψ_h2) / log(z / z0 - ψ_h)
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

function calc_wind(z, z0, κ, u_star, b)
    return (u_star / κ) * log(z / z0 + b)
end

"""
    calc_ρ_cp(T_mean)

Compute the volumetric heat capacity of air (ρ·cₚ) as a function of mean temperature.

# Arguments
- `T_mean`: Mean air temperature (`Unitful.Temperature`), in Kelvin.

# Returns
- Volumetric heat capacity (`cal / (cm³·K)`).

This is a simplified empirical regression based only on temperature,
without accounting for moisture or elevation effects.
"""
function calc_ρ_cp(T_mean)
    return u"(cal*g)/(g*cm^3*K)" * (0.08472 / ustrip(u"K", T_mean))
end


"""
    calc_ρ_cp(T_mean, elevation, relative_humidity)

Compute the volumetric heat capacity of moist air (ρ·cₚ) given temperature,
elevation, and relative humidity.

# Arguments
- `T_mean`: Mean air temperature (`Unitful.Temperature`), in Kelvin.
- `elevation`: Elevation above sea level (with units of length).
- `relative_humidity`: Relative humidity (fraction between 0 and 1).

# Returns
- Volumetric heat capacity (`cal / (cm³·K)`).

Uses `dry_air_properties` to compute air density (ρ) and 
`wet_air_properties` to compute specific heat capacity (cₚ).
"""
function calc_ρ_cp(T_mean, elevation, relative_humidity)
    dry_air_out = dry_air_properties(u"K"(T_mean); elevation)
    wet_air_out = wet_air_properties(u"K"(T_mean); rh = relative_humidity)
    ρ = dry_air_out.ρ_air
    c_p = wet_air_out.c_p
    return u"(cal*g)/(g*cm^3*K)"(ρ * c_p)
end


"""
    calc_φ_m(z, γ, L_Obukhov)

Stability correction function φ for momentum in Monin–Obukhov similarity theory (MOST).

# Arguments
- `z`: Height above surface (with units of length).
- `γ`: Empirical constant (dimensionless, often ≈16).
- `L_Obukhov`: Monin–Obukhov length (with units of length).

# Returns
- Dimensionless stability correction factor φ.

This corresponds to the Businger–Dyer formulation for unstable stratification:

φₘ = (1 - γ z / L)^(1/4)

# References
- Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971).
  Flux–profile relationships in the atmospheric surface layer.
  *Journal of the Atmospheric Sciences*, 28(2), 181–189.
- Dyer, A. J. (1974). A review of flux–profile relationships.
  *Boundary-Layer Meteorology*, 7(3), 363–372.
"""
function calc_φ_m(z, γ, L_Obukhov)
    return (1.0 - min(1.0, γ * (z / L_Obukhov)))^(1//4)
end


"""
    calc_ψ_m(x)

Stability correction function ψₘ for momentum under unstable atmospheric stratification,
used in Monin–Obukhov similarity theory.

# Arguments
- `x`: Dimensionless argument, typically `(1 - γ z / L)^(1/4)`.

# Returns
- Correction factor ψₘ (dimensionless).

This is the Businger–Dyer form for momentum:

ψₘ(x) = 2 ln((1 + x) / 2) + ln((1 + x²) / 2) - 2 atan(x) + π/2

# References
- Businger et al. (1971).
- Dyer (1974).
"""
function calc_ψ_m(x)
    return 2.0 * log((1.0 + x) / 2.0) + log((1.0 + x^2) / 2.0) - 2.0 * atan(x) + π / 2.0
end


"""
    calc_ψ_h(x)

Stability correction function ψ_h for heat and moisture under unstable conditions,
used in Monin–Obukhov similarity theory.

# Arguments
- `x`: Dimensionless argument, typically `(1 - γ z / L)^(1/4)`.

# Returns
- Correction factor ψ_h (dimensionless).

This is the Businger–Dyer form for scalars:

ψ_h(x) = 2 ln((1 + x²) / 2)

# References
- Businger et al. (1971).
- Dyer (1974).
"""
function calc_ψ_h(x)
    return 2.0 * log((1 + x^2.0) / 2.0)
end


"""
    sublayer_stanton(z0, u_star)

Compute the Stanton number for the viscous sublayer.
"""
function sublayer_stanton(z0, u_star)
    return 0.62 / (ustrip(u"cm", z0) * ustrip(u"cm/minute", u_star) / 12)^(9//20)
end

"""
    bulk_stanton(log_z_ratio, z, L_Obukhov)

Compute the bulk Stanton number depending on stability (Monin-Obukhov length).
"""
function bulk_stanton(log_z_ratio, z, L_Obukhov)
    if L_Obukhov > 0.0u"cm"  # stable
        return 0.64 / log_z_ratio
    else                 # unstable
        return (0.64 / log_z_ratio) * (1 - 0.1 * z / L_Obukhov)
    end
end

"""
    convective_flux(ρ_cp, ΔT, u_star, St_bulk, St_sublayer)

Compute convective heat flux given bulk and sublayer Stanton numbers.
"""
function convective_flux(ρ_cp, ΔT, u_star, bulk_stanton_number, sublayer_stanton_number)
        return ρ_cp * ΔT * u_star * bulk_stanton_number / (1 + bulk_stanton_number / sublayer_stanton_number)
end

"""
    calc_Obukhov_length(T_ref_height, T_surface, Ū_ref_height, z, z0, ρcpTκg, κ, log_z_ratio, ΔT, ρ_cp, 
                         max_iter=500, tol=1e-2)

Iteratively solve for Monin-Obukhov length and convective heat flux.
"""
function calc_Obukhov_length(T_ref_height, T_surface, Ū_ref_height, z, z0, ρcpTκg, κ, 
    log_z_ratio, ΔT, ρ_cp; max_iter=500, tol=1e-2)
    L_Obukhov = -30.0u"cm" # initial Monin-Obukhov length cm
    γ = 16.0 # -

    # conversions
    z = u"cm"(z)
    z0 = u"cm"(z0)
    T_ref_height = u"K"(T_ref_height)
    T_surface = u"K"(T_surface)

    # initialise
    Q_convection = nothing
    effective_stanton_number = nothing
    bulk_stanton_number = nothing
    sublayer_stanton_number = nothing
    u_star = nothing
    δ = 1.0
    count = 0
    while δ > tol && count < max_iter
        count += 1
        φ_m = calc_φ_m(z, γ, L_Obukhov)
        ψ_m = calc_ψ_m(φ_m)
        u_star = κ * Ū_ref_height / (log(z / z0) - ψ_m)
        sublayer_stanton_number = sublayer_stanton(z0, u_star)
        bulk_stanton_number = bulk_stanton(log_z_ratio, z, L_Obukhov)
        Q_convection = convective_flux(ρ_cp, ΔT, u_star, bulk_stanton_number, sublayer_stanton_number)
        L_Obukhov_new = ρcpTκg * u_star^3 / Q_convection
        δ = abs((L_Obukhov_new - L_Obukhov) / L_Obukhov)
        L_Obukhov = L_Obukhov_new
    end
    T_z0 = (T_ref_height * bulk_stanton_number + T_surface * sublayer_stanton_number) / (bulk_stanton_number + sublayer_stanton_number)
    return (; L_Obukhov=u"m"(L_Obukhov), sublayer_stanton_number, effective_stanton_number, bulk_stanton_number, u_star, Q_convection, T_z0)
end

