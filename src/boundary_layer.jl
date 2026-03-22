function allocate_profile(heights)
    wind_speed = similar(heights, typeof(0.0u"m/s")) # output wind speeds
    height_array = similar(heights, typeof(0.0u"m"))
    height_array[end:-1:begin] .= heights 
    air_temperature = similar(heights, typeof(0.0u"K")) # output temperatures, need to do this otherwise get InexactError
    relative_humidity = similar(heights, Float64) # output relative humidities
    return (; heights, height_array, air_temperature, wind_speed, relative_humidity)
end

"""
    atmospheric_surface_profile(; kwargs...)

Compute vertical profiles of wind speed, air temperature, and relative humidity 
in the atmospheric surface layer, using Monin–Obukhov similarity theory (MOST).

This function reproduces the subroutine in `MICRO.f/get_profile.R` from **NicheMapR**, ported to Julia.  
It calculates the microclimate profiles above the ground (or canopy) at specified heights,
based on measured values at a reference height and computed or measured soil surface temperature, together
with surface roughness parameters. Zenith angle and a maximum allowed surface temperature are used
to assess whether conditions are stable or unstable.

# Keyword Arguments

- `z0::Quantity=0.004u"m"`: roughness length (surface aerodynamic roughness).
- `karman_constant::Float64=0.4`: von Kármán constant.
- `dyer_constant::Float=16, coefficient from Dyer and Hicks for Φ_m (momentum), γ
- `heights::Vector{Quantity}`: Requested heights above the surface, the last being the reference height.
- `reference_temperature::Quantity=27.78u"°C"`: Air temperature at the reference height.
- `reference_wind_speed::Quantity=2.75u"m/s"`: Wind speed at the reference height.
- `relative_humidity::Float64=49.0`: Relative humidity at the reference height (fractional).
- `surface_temperature::Quantity=48.59u"°C"`: Soil or surface temperature.
- `zenith_angle::Quantity=21.5u"°"`: Solar zenith angle.
- `elevation::Quantity=0.0u"m"`: Elevation above sea level.

# Returns
Named tuple with fields:
- `wind_speed`: Wind speed profile at each height (`m/s`).
- `air_temperature`: Air temperature profile at each height (`K`).
- `relative_humidity`: Relative humidity (fractional) at each height.
- `convective_heat_flux`: Convective heat flux (`W/m²`).
- `u_star`: Friction velocity (`m/s`).

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
profile = atmospheric_surface_profile(
    reference_temperature = 25u"°C",
    reference_wind_speed = 2.0u"m/s",
    relative_humidity = 0.6,
    surface_temperature = 35u"°C",
    zenith_angle = 45u"°"
)

profile.air_temperature  # vertical profile of air temperatures
profile.wind_speed       # vertical profile of wind speeds
```
"""
atmospheric_surface_profile(; heights=DEFAULT_HEIGHTS, kw...) =
    atmospheric_surface_profile!(allocate_profile(heights); kw...)
function atmospheric_surface_profile!(buffers;
    micro_terrain,
    environment_instant,
    surface_temperature,
    vapour_pressure_equation=GoffGratch(),
)
    (; roughness_height, karman_constant, dyer_constant, elevation) = micro_terrain
    (; atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity, zenith_angle) = environment_instant

    (; heights, height_array, air_temperature, wind_speed, relative_humidity) = buffers
    N_heights = length(heights)
    if minimum(heights) < roughness_height
        throw(ArgumentError("The minimum height is not greater than the roughness height."))
    end
    reference_height = last(heights)

    reference_temp = u"K"(reference_temperature)
    surface_temp = u"K"(surface_temperature)
    κ = karman_constant
    γ = dyer_constant
    z = reference_height
    z0 = roughness_height
    v_ref_height = reference_wind_speed

    # define air heights
    N_heights = length(heights)
    relative_humidity = zeros(Float64, N_heights) # output relative humidities
    wind_speed[1] = v_ref_height
    air_temperature[1] = reference_temp

    # compute rcptkg (was a constant in original Fortran version)
    # dry_air_out = dry_air_properties(u"K"(reference_temperature), P_atmos)
    # wet_air_out = wet_air_properties(u"K"(reference_temperature), rh=reference_humidity, P_atmos)
    # ρ = dry_air_out.ρ_air
    # c_p = wet_air_out.specific_heat
    # TODO make this work with SI units
    #ρcpTκg = u"cal*minute^2/cm^4"(ρ * c_p * T_ref_height / (κ * g_n))
    ρcpTκg = u"J*s^2/m^4"(6.003e-8u"cal*minute^2/cm^4")
    
    log_z_ratio = log(z / z0 + 1)
    ΔT = reference_temp - surface_temp
    mean_temp = (surface_temp + reference_temp) / 2
    # TODO call calc_ρ_cp method specific to elevation and RH in final version but do it this way for NicheMapR comparison
    ρ_cp = calc_ρ_cp(mean_temp)#, elevation, reference_humidity)


    # TODO name and explain this check, why `|| zenith_angle`
    if reference_temp ≥ surface_temp || zenith_angle ≥ 90°
        u_star = calc_u_star(; reference_wind_speed, log_z_ratio, κ)
        convective_heat_flux = calc_convection(; u_star, log_z_ratio, ΔT, ρ_cp, z0)
        for i in 2:N_heights
            wind_speed[i] = calc_wind(height_array[i], z0, κ, u_star, 1.0)
            roughness_height_temp = (reference_temp * bulk_stanton(log_z_ratio) + surface_temp * sublayer_stanton(z0, u_star)) / (bulk_stanton(log_z_ratio) + sublayer_stanton(z0, u_star))
            air_temperature[i] = roughness_height_temp + (reference_temp - roughness_height_temp) * log(height_array[i] / z0 + 1.0) / log_z_ratio
        end
    else
        obukhov_length = -0.3u"m" # initialise Obukhov length
        # TODO just pass the environment_instant through here
        Obukhov_out = calc_Obukhov_length(reference_temp, surface_temp, v_ref_height, z0, z, ρcpTκg, κ, log_z_ratio, ΔT, ρ_cp; max_iter=30, tol=1e-2)
        obukhov_length = Obukhov_out.obukhov_length
        roughness_height_temp = Obukhov_out.roughness_height_temperature
        convective_heat_flux = Obukhov_out.convective_heat_flux
        u_star = Obukhov_out.u_star
        ψ_h = Obukhov_out.ψ_h
        for i in 2:N_heights
            φ_m1 = calc_φ_m(height_array[i], γ, obukhov_length)
            ψ_m1 = calc_ψ_m(φ_m1)
            ψ_h2 = calc_ψ_h(φ_m1)
            h_ratio = height_array[i] / z0  # dimensionless h/z0
            # Clamp log arguments to a small positive value to prevent NaN when
            # stability corrections exceed h/z0 at near-surface heights (e.g. 0.01 m).
            wind_log_arg = max(h_ratio - ψ_m1, 1e-6)
            wind_speed[i] = (u_star / κ) * log(wind_log_arg)
            temp_log_arg = max(h_ratio - ψ_h2, 1e-6)
            air_temperature[i] = roughness_height_temp + (reference_temp - roughness_height_temp) * log(temp_log_arg) / log(z / z0 - ψ_h)
        end
    end
    wind_speed = reverse(wind_speed)
    air_temperature = reverse(air_temperature)
    reference_vapor_pressure = wet_air_properties(reference_temp, reference_humidity, atmospheric_pressure; vapour_pressure_equation).vapour_pressure
    relative_humidity .= clamp.(reference_vapor_pressure ./ vapour_pressure.(Ref(vapour_pressure_equation), air_temperature) .* 1.0, 0.0, 1.0)

    return (;
        wind_speed,
        air_temperature,
        relative_humidity,
        convective_heat_flux=u"W/m^2"(convective_heat_flux),
        u_star,
    )
end

"""
    calc_ρ_cp(mean_temperature)

Compute the volumetric heat capacity of air (ρ·cₚ) as a function of mean temperature.

# Arguments
- `mean_temperature`: Mean air temperature (`Unitful.Temperature`), in Kelvin.

# Returns
- Volumetric heat capacity (`cal / (cm³·K)`).

This is a simplified empirical regression based only on temperature,
without accounting for moisture or elevation effects.
"""
function calc_ρ_cp(mean_temperature)
    return u"J/m^3/K"(u"(cal*g)/(g*cm^3*K)" * (0.08472 / ustrip(u"K", mean_temperature)))
end

"""
    calc_ρ_cp(mean_temperature, elevation, relative_humidity, atmospheric_pressure)

Compute the volumetric heat capacity of moist air (ρ·cₚ) given temperature,
elevation, and relative humidity.

# Arguments
- `mean_temperature`: Mean air temperature (`Unitful.Temperature`), in Kelvin.
- `elevation`: Elevation above sea level (with units of length).
- `relative_humidity`: Relative humidity (fraction between 0 and 1).
- `atmospheric_pressure`: Atmospheric pressure.

# Returns
- Volumetric heat capacity (`cal / (cm³·K)`).

Uses `dry_air_properties` to compute air density (ρ) and
`wet_air_properties` to compute specific heat capacity (cₚ).
"""
function calc_ρ_cp(mean_temperature, elevation, relative_humidity, atmospheric_pressure; vapour_pressure_equation=GoffGratch())
    dry_air_out = dry_air_properties(u"K"(mean_temperature), atmospheric_pressure)
    wet_air_out = wet_air_properties(u"K"(mean_temperature), relative_humidity, atmospheric_pressure; vapour_pressure_equation)
    air_density = dry_air_out.density
    air_heat_capacity = wet_air_out.specific_heat
    return air_density * air_heat_capacity
end

"""
    calc_u_star(; reference_wind_speed, log_z_ratio, κ=0.4)

Compute the friction velocity (u*) from a reference wind speed using the
logarithmic wind profile.

# Arguments
- `reference_wind_speed::Quantity{<:Real,𝐋/𝐓}`: Wind speed at the reference height (e.g. `m/s`, `cm/min`).
- `log_z_ratio::Real`: Precomputed log height ratio, typically `log(z/z0 + 1.0)`.
- `κ::Real`: von Kármán constant (default = 0.4).

# Returns
- Friction velocity `u_star::Quantity{<:Real,𝐋/𝐓}`.

# See also
[`calc_convection`](@ref), [`calc_wind`](@ref)
"""
function calc_u_star(; reference_wind_speed, log_z_ratio, κ=0.4)
    v_ref_height = reference_wind_speed
    return κ * v_ref_height / log_z_ratio
end

"""
    calc_wind(z, z0, κ, u_star, b)

Calculate wind speed at height `z` using the logarithmic wind profile.

# Arguments
- `z::Quantity{<:Real,𝐋}`: Height above the surface (e.g. `m`, `cm`).
- `z0::Quantity{<:Real,𝐋}`: Roughness length (e.g. `m`, `cm`).
- `κ::Real`: von Kármán constant.
- `u_star::Quantity{<:Real,𝐋/𝐓}`: Friction velocity.
- `b::Real`: Offset term (e.g. `1.0` for neutral stability, or stability correction).

# Returns
- Wind speed at height `z::Quantity{<:Real,𝐋/𝐓}`.

# See also
[`calc_u_star`](@ref), [`calc_convection`](@ref)
"""
function calc_wind(z, z0, κ, u_star, b)
    return (u_star / κ) * log(z / z0 + b)
end


"""
    calc_convection(; u_star, log_z_ratio, ΔT, ρ_cp, z0)

Calculate the convective heat flux (sensible heat exchange between surface and air).

# Arguments
- `u_star::Quantity{<:Real,𝐋/𝐓}`: Friction velocity (e.g. `m/s`, `cm/min`).
- `log_z_ratio::Real`: Precomputed logarithmic height ratio, typically `log(z/z0 + 1.0)`.
- `ΔT::Quantity{<:Real,Θ}`: Temperature difference between reference air and surface (Kelvin).
- `ρ_cp::Quantity{<:Real,(𝐌*𝐋^-1*𝐓^-2)}`: Volumetric heat capacity of air (e.g. `J/m³/K`, `cal/cm³/K`).
- `z0::Quantity{<:Real,𝐋}`: Surface roughness length (length).

# Returns
- Convective heat flux as `Quantity{<:Real,(𝐌*𝐓^-3)}` (e.g. `W/m²`, `cal/min/cm²`).

Uses bulk and sublayer Stanton numbers to account for turbulence near the surface.

# See also
[`calc_u_star`](@ref), [`calc_wind`](@ref), [`sublayer_stanton`](@ref), [`bulk_stanton`](@ref), [`convective_flux`](@ref)
"""
function calc_convection(; u_star, log_z_ratio, ΔT, ρ_cp, z0)
    sublayer_stanton_number = sublayer_stanton(z0, u_star)
    bulk_stanton_number = bulk_stanton(log_z_ratio)
    return convective_flux(ρ_cp, ΔT, u_star, bulk_stanton_number, sublayer_stanton_number)
end

"""
    convective_flux(ρ_cp, ΔT, u_star, St_bulk, St_sublayer)

Compute convective heat flux given bulk and sublayer Stanton numbers.
"""
function convective_flux(ρ_cp, ΔT, u_star, bulk_stanton_number, sublayer_stanton_number)
    return ρ_cp * ΔT * u_star * bulk_stanton_number / (1 + bulk_stanton_number / sublayer_stanton_number)
end


"""
    sublayer_stanton(z0, u_star)

Compute the Stanton number for the viscous sublayer.
"""
function sublayer_stanton(z0, u_star)
    return 0.62 / (ustrip(u"cm", z0) * ustrip(u"cm/minute", u_star) / 12)^(9//20)
end

"""
    bulk_stanton(log_z_ratio)

Compute the bulk Stanton number for stable conditions.
"""
function bulk_stanton(log_z_ratio)
    return 0.64 / log_z_ratio
end

"""
    bulk_stanton(log_z_ratio, z, obukhov_length)

Compute the bulk Stanton number for unstable conditions.
"""
function bulk_stanton(log_z_ratio, z, obukhov_length)
    return (0.64 / log_z_ratio) * (1 - 0.1 * z / obukhov_length)
end


"""
    calc_φ_m(z, γ, obukhov_length)

Stability correction function φ for momentum in Monin–Obukhov similarity theory (MOST).

# Arguments
- `z`: Height above surface (with units of length).
- `γ`: Empirical constant (dimensionless, often ≈16).
- `obukhov_length`: Monin–Obukhov length (with units of length).

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
function calc_φ_m(z, γ, obukhov_length)
    #return (1.0 - min(1.0, γ * (z / obukhov_length)))^(1//4)
    return (1.0 - γ * z / obukhov_length)^(1//4)
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
    return 2.0 * log((1.0 + x^2.0) / 2.0)
end


"""
    calc_Obukhov_length(reference_temp, surface_temp, v_ref_height, z, z0, ρcpTκg, κ, log_z_ratio, ΔT, ρ_cp,
                         γ=16.0, max_iter=500, tol=1e-2)

Iteratively solve for Monin-Obukhov length and convective heat flux.
"""
function calc_Obukhov_length(
    reference_temp, surface_temp, v_ref_height, z0, z, ρcpTκg, κ, log_z_ratio, ΔT, ρ_cp;
    γ=16.0, max_iter=30, tol=1e-2
)
    obukhov_length = -0.3u"m" # initial Monin-Obukhov length

    # initialise with zeros
    convective_heat_flux = 0.0u"W/m^2"
    bulk_stanton_number = 0.0
    sublayer_stanton_number = 0.0
    ψ_h = 0.0
    φ_m = 0.0
    u_star = 0.0u"m/s"
    obukhov_length_new = 0.0u"m"

    relative_error = 1.0
    count = 0
    just_above_zero = 1.0e-6
    while relative_error > tol && count < max_iter
        count += 1
        φ_m = calc_φ_m(z, γ, obukhov_length)
        ψ_m = calc_ψ_m(φ_m)
        ψ_h = calc_ψ_h(φ_m)
        log_ratio_corrected = log(z / z0) - ψ_m
        if log_ratio_corrected <= 0.0
            log_ratio_corrected = just_above_zero
        end
        u_star = κ * v_ref_height / log_ratio_corrected
        if u_star < just_above_zero * 1u"m/s"
            u_star = just_above_zero * 1u"m/s"
        end
        sublayer_stanton_number = sublayer_stanton(z0, u_star)
        bulk_stanton_number = bulk_stanton(log_ratio_corrected, z, obukhov_length)
        convective_heat_flux = convective_flux(ρ_cp, ΔT, u_star, bulk_stanton_number, sublayer_stanton_number)
        obukhov_length_new = ρcpTκg * u_star^3 / convective_heat_flux
        relative_error = abs((obukhov_length_new - obukhov_length) / obukhov_length)
        obukhov_length = obukhov_length_new
    end

    roughness_height_temperature = (reference_temp * bulk_stanton_number + surface_temp * sublayer_stanton_number) / (bulk_stanton_number + sublayer_stanton_number)

    return (;
        obukhov_length=u"m"(obukhov_length),
        sublayer_stanton_number,
        bulk_stanton_number,
        u_star,
        ψ_h,
        convective_heat_flux=u"W/m^2"(convective_heat_flux),
        roughness_height_temperature,
    )
end

