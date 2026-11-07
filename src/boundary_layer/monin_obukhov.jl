"""
    MoninObukhov(; karman_constant=0.4, stability_functions=ChengBrutsaert(),
                   min_friction_velocity=0.1u"m/s", thermal_roughness_model=SublayerStantonRoughness(),
                   stability_limit=NoStabilityLimit())

Monin–Obukhov similarity theory boundary-layer formulation, with the
stability corrections of `stability_functions` ([`ChengBrutsaert`](@ref) or
[`BusingerDyer`](@ref)).
`stability_limit` bounds the canopy-top stability ([`NoStabilityLimit`](@ref) or
[`MomentumCorrectionLimit`](@ref)).

The thermal_roughness_model is a z0-to-z0h correction for the sensible heat transfer
coefficient. Default is SublayerStantonRoughness(); canopy uses ScalarRoughnessRatio().
"""
@kwdef struct MoninObukhov{KC,SF,MINUF,TRM,SL} <: AbstractBoundaryLayerModel
    karman_constant::KC = 0.4
    stability_functions::SF = ChengBrutsaert()
    min_friction_velocity::MINUF = 0.1u"m/s" # stability clamp
    thermal_roughness_model::TRM = SublayerStantonRoughness()
    stability_limit::SL = NoStabilityLimit()
end

"""
    AbstractStabilityLimit

How far the atmospheric stability above a canopy is allowed to go.
"""
abstract type AbstractStabilityLimit end

"""
    NoStabilityLimit()

Atmospheric stability above a canopy has no limit.
"""
struct NoStabilityLimit <: AbstractStabilityLimit end

"""
    MomentumCorrectionLimit(; fraction=0.9)

Limits the atmospheric stability above a canopy, so that stability changes the friction
velocity by at most a factor `1 ± fraction`.

On calm nights a canopy can cool faster than stability reduces its heat flux, and then no
stability agrees with that flux. This limit gives one. Stability approaches the limit smoothly.
"""
@kwdef struct MomentumCorrectionLimit{F} <: AbstractStabilityLimit
    fraction::F = 0.9
end

"""
    AbstractStabilityFunctions

Monin–Obukhov stability corrections: variants implement [`momentum_stability_correction`](@ref),
[`heat_stability_correction`](@ref) and [`heat_gradient_stability_factor`](@ref) of the stability parameter `ζ = z/L`.
"""
abstract type AbstractStabilityFunctions end

"""
    BusingerDyer(; dyer_constant=16.0, stable_beta=4.7, turbulent_prandtl_number=0.74,
        stable_heat_gradient_coefficient=6.0, min_stable_heat_gradient=0.5, max_stable_heat_gradient=1.5)

NicheMapR's stability corrections. Unstable: Businger–Dyer/Paulson with
`dyer_constant`. Stable: the linear form `ψ_m = -stable_beta·ζ`,
`ψ_h = -(stable_beta/turbulent_prandtl_number)·ζ` (Businger et al. 1971),
clamped to `[-4, 3]`, and a separate saturating diffusivity multiplier
`Φ_h = clamp(1 + stable_heat_gradient_coefficient·ζ/(1+ζ), min_stable_heat_gradient, max_stable_heat_gradient)`.
`ζ` is bounded to ±10. The clamps put kinks in every function of stability;
[`ChengBrutsaert`](@ref) has none.

# References
Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971).
Flux–profile relationships in the atmospheric surface layer.
*Journal of the Atmospheric Sciences*, 28(2), 181–189.

Dyer, A. J. (1974). A review of flux–profile relationships.
*Boundary-Layer Meteorology*, 7(3), 363–372.
"""
@kwdef struct BusingerDyer{DC,SB,PR,PHC,MINPH,MAXPH} <: AbstractStabilityFunctions
    dyer_constant::DC = 16.0
    stable_beta::SB = 4.7
    turbulent_prandtl_number::PR = 0.74
    stable_heat_gradient_coefficient::PHC = 6.0
    min_stable_heat_gradient::MINPH = 0.5
    max_stable_heat_gradient::MAXPH = 1.5
end

"""
    ChengBrutsaert(; momentum_dyer_constant=16.0, heat_dyer_constant=16.0, momentum_scale=6.1, momentum_exponent=2.5, heat_scale=5.3, heat_exponent=1.1)

Smooth stability corrections, without clamps. Unstable: Businger–Dyer/Paulson
with `momentum_dyer_constant` in `ψ_m`, `heat_dyer_constant` in `ψ_h` and `Φ_h`;
the slopes at neutral are `-momentum_dyer_constant/4` and `-heat_dyer_constant/2`. Stable: Cheng & Brutsaert (2005),
`ψ_m = -a·ln(ζ + (1 + ζᵇ)^(1/b))` with `a = momentum_scale` and `b = momentum_exponent`,
`ψ_h = -c·ln(ζ + (1 + ζᵈ)^(1/d))` with `c = heat_scale` and `d = heat_exponent`, which
level off into the z-less limit rather than growing linearly, so the bulk
Richardson number has no critical value and every stable state has a
Monin–Obukhov solution.

# References
Cheng, Y., & Brutsaert, W. (2005). Flux-profile relationships for wind speed
and temperature in the stable atmospheric boundary layer.
*Boundary-Layer Meteorology*, 114, 519–538.
"""
@kwdef struct ChengBrutsaert{DM,DH,A,B,C,D} <: AbstractStabilityFunctions
    momentum_dyer_constant::DM = 16.0
    heat_dyer_constant::DH = 16.0
    momentum_scale::A = 6.1
    momentum_exponent::B = 2.5
    heat_scale::C = 5.3
    heat_exponent::D = 1.1
end

"""
    AbstractThermalRoughnessModel

How the momentum roughness length `z0` relates to the scalar (heat)
transfer resistance -- the `z0`-to-`z0h` correction, `kB⁻¹` in the
literature's usual notation.
"""
abstract type AbstractThermalRoughnessModel end

"""
    SublayerStantonRoughness()

[`sublayer_stanton`](@ref)'s `Re*`-dependent correction, combined in series
with the bulk (log-law) Stanton number via [`convective_flux`](@ref).
Default; appropriate at the bare-ground/short-roughness scale.
"""
struct SublayerStantonRoughness <: AbstractThermalRoughnessModel end

"""
    ScalarRoughnessRatio(; ratio=0.2)

Fixed scalar-to-momentum roughness ratio `z0h = ratio·z0` (`kB⁻¹ =
-log(ratio)`, constant, not `Re*`-dependent).
"""
@kwdef struct ScalarRoughnessRatio{R} <: AbstractThermalRoughnessModel
    ratio::R = 0.2
end

const _MIN_LOGLAW_RATIO = 5.0 # avoid the log-law inversion becoming ill-conditioned

function allocate_profile(::MoninObukhov, heights)
    wind_speed = similar(heights, typeof(0.0u"m/s")) # output wind speeds
    air_temperature = similar(heights, typeof(0.0u"K")) # output temperatures, need to do this otherwise get InexactError
    relative_humidity = similar(heights, Float64) # output relative humidities
    obukhov_length_prev = Ref(-0.3u"m")  # warm-start across timesteps
    warned_below_roughness = Ref(false)  # emit warning at most once per buffer
    warned_above_reference = Ref(false)
    return (; heights, air_temperature, wind_speed, relative_humidity,
              obukhov_length_prev, warned_below_roughness, warned_above_reference)
end

"""
    atmospheric_surface_profile(; kwargs...)

Compute vertical profiles of wind speed, air temperature, and relative humidity
using Monin–Obukhov similarity theory (MOST). Profiles are returned at all requested
`heights`; heights below the roughness length return zero wind and linearly interpolated
temperature; heights above `reference_height` are extrapolated via the log-law.

# Keyword Arguments

- `heights`: Heights at which profiles are returned (default: `DEFAULT_HEIGHTS`).
- `reference_height`: Measurement height for met inputs (default: `last(heights)`).
- `surface_temperature`: Soil/surface temperature.
- `site`: `Site` holding `roughness_height` and `elevation`.
- `environment_instant`: Named tuple with `reference_temperature`, `reference_wind_speed`,
  `reference_humidity`, `atmospheric_pressure`, and `zenith_angle`.

# Returns
Named tuple: `wind_speed`, `air_temperature`, `relative_humidity` (profiles at each height),
`convective_heat_flux` (`W/m²`), `friction_velocity` (`m/s`), `obukhov_length`, `inverse_obukhov_length`.

# References
- Businger et al. (1971). *J. Atmos. Sci.* 28, 181–189.
- Dyer (1974). *Bound.-Layer Meteor.* 7, 363–372.

# Example

```julia
using Microclimate, Unitful

boundary_layer_model = MoninObukhov()
site = example_site()
environment_instant = (;
    reference_temperature = 25u"°C",
    reference_wind_speed  = 2.0u"m/s",
    reference_humidity    = 0.6,
    atmospheric_pressure  = 101.325u"kPa",
    zenith_angle          = 45u"°",
)
profile = atmospheric_surface_profile(boundary_layer_model;
    heights          = [0.01, 0.5, 2.0, 5.0]u"m",  # 5 m is above the reference height
    reference_height = 2.0u"m",
    site,
    environment_instant,
    surface_temperature = 35u"°C",
)
profile.wind_speed
profile.air_temperature
```
"""
atmospheric_surface_profile(boundary_layer_model::MoninObukhov; heights=DEFAULT_HEIGHTS, reference_height=last(heights), kw...) =
    atmospheric_surface_profile!(boundary_layer_model, allocate_profile(boundary_layer_model, heights); reference_height, kw...)
"""
    atmospheric_surface_profile!(boundary_layer_model::MoninObukhov, buffers; site, environment_instant,
                                  surface_temperature, vapour_pressure_equation=GoffGratch(),
                                  reference_height=last(buffers.heights),
                                  displacement_height=0.0u"m", roughness_length=site.roughness_height)

`displacement_height`/`roughness_length` shift the effective log-profile
origin to a canopy-displaced height with a canopy-derived roughness length,
instead of bare ground. Defaults (`0.0u"m"` / `site.roughness_height`)
give the bare-ground profile. Heights are measured from the true ground (as
in `MicroModel.atmosphere.heights`); `displacement_height` is subtracted internally
wherever the log-law needs height above the displaced origin.

`temperature_anchor_height`: by default (`nothing`) the temperature profile
is anchored at `z0` using `roughness_height_temperature`. Pass an absolute height 
(e.g. a canopy's `canopy_height`) to anchor directly on `surface_temperature`
*at that height* for a plain two-point MOST log-law interpolation between
`(temperature_anchor_height, surface_temperature)` and
`(reference_height, reference_temperature)`.

`surface_relative_humidity` (only used when `temperature_anchor_height` is
given): the known relative humidity *at* `temperature_anchor_height`,
anchoring vapour pressure the same two-point way.

`stability_variable` ([`ObukhovLength`](@ref) or [`InverseObukhovLength`](@ref))
selects what the stability iteration solves for. Both are returned, as
`obukhov_length` and `inverse_obukhov_length`. Passing `stability` (as
`stability_variable`) skips the iteration and evaluates the profile at that
stability; `stability_update` in the result is then what the Obukhov
relation gives back for it, equal to `stability` when self-consistent.
"""
function atmospheric_surface_profile!(boundary_layer_model::MoninObukhov, buffers;
    site,
    environment_instant,
    surface_temperature,
    vapour_pressure_equation=GoffGratch(),
    reference_height=last(buffers.heights),
    displacement_height=0.0u"m",
    roughness_length=site.roughness_height,
    thermal_roughness_model=boundary_layer_model.thermal_roughness_model,
    temperature_anchor_height=nothing,
    surface_relative_humidity=nothing,
    stability_variable=ObukhovLength(),
    stability=nothing,
)
    (; elevation) = site
    (; karman_constant, stability_functions, min_friction_velocity) = boundary_layer_model
    (; atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity) = environment_instant

    (; heights, air_temperature, wind_speed, relative_humidity,
       obukhov_length_prev, warned_below_roughness, warned_above_reference) = buffers
    N_heights = length(heights)
    # `displacement_height`/`roughness_length` shift the log-profile origin to a
    # canopy-displaced height. Heights are measured from the true ground, so all
    # height-above-origin terms subtract `displacement_height`.
    reference_height = u"m"(reference_height)
    reference_height_above_origin = reference_height - displacement_height
    # Pin a reference-height slot to the measured values.
    # TODO remove with the min_friction_velocity floor (#137): the floor is all that makes
    # the profile differ from the forcing there.
    reference_index = findfirst(==(reference_height), heights)

    if !warned_below_roughness[] && minimum(heights) - displacement_height < roughness_length
        @warn """Some heights are within the roughness length ($roughness_length) of the displaced origin
        ($displacement_height above ground) -- below MOST's valid range (expected when canopy is present).
    wind → 0; air temperature interpolated linearly between T_surface and T_z0.
    Clamp downstream for a non-zero wind floor (e.g. convection terms)."""
        warned_below_roughness[] = true
    end
    if !warned_above_reference[] && maximum(heights) > reference_height
        @warn """Some requested heights exceed the reference measurement height ($reference_height).
    Profiles will be extrapolated above the measurement level using MOST."""
        warned_above_reference[] = true
    end

    reference_temperature = u"K"(reference_temperature)
    surface_temperature = u"K"(surface_temperature)
    κ = karman_constant

    # ρ·cₚ·T/(κ·g), the factor in the Obukhov relation L = factor·u*³/H. A constant, as in
    # the original Fortran. TODO calculate it from the air properties.
    obukhov_length_factor = u"J*s^2/m^4"(6.003e-8u"cal*minute^2/cm^4")

    temperature_difference = reference_temperature - surface_temperature
    mean_temperature = (surface_temperature + reference_temperature) / 2
    # TODO call calc_air_volumetric_heat_capacity method specific to elevation and RH in final version but do it this way for NicheMapR comparison
    air_volumetric_heat_capacity = calc_air_volumetric_heat_capacity(mean_temperature)#, elevation, reference_humidity)
    kinematic_viscosity = dry_air_properties(mean_temperature, atmospheric_pressure).kinematic_viscosity

    stability_result = surface_layer_stability(stability_variable, stability, obukhov_length_prev, reference_temperature,
        surface_temperature, reference_wind_speed, roughness_length, reference_height_above_origin, obukhov_length_factor, karman_constant, temperature_difference, air_volumetric_heat_capacity, kinematic_viscosity;
        stability_functions, min_friction_velocity, thermal_roughness_model)
    stability = stability_result.value
    roughness_height_temperature = stability_result.roughness_height_temperature
    convective_heat_flux = stability_result.convective_heat_flux
    friction_velocity = stability_result.friction_velocity
    ψ_h = stability_result.heat_correction
    # See docstring: default anchors on roughness_length/roughness_height_temperature;
    # temperature_anchor_height switches to a direct two-point interpolation.
    use_anchor = !isnothing(temperature_anchor_height)
    anchor_height_above_origin = use_anchor ? u"m"(temperature_anchor_height) - displacement_height : u"m"(roughness_length)
    anchor_temperature = use_anchor ? surface_temperature : roughness_height_temperature
    ψ_h_anchor = heat_stability_correction(stability_functions, stability_parameter_at_height(stability_variable, anchor_height_above_origin, stability))
    ψ_m_roughness = momentum_stability_correction(stability_functions, stability_parameter_at_height(stability_variable, roughness_length, stability))
    temperature_denominator = log(reference_height_above_origin / anchor_height_above_origin) + ψ_h_anchor - ψ_h

    reference_vapor_pressure = wet_air_properties(reference_temperature, reference_humidity, atmospheric_pressure; vapour_pressure_equation).vapour_pressure
    # No surface_relative_humidity: constant vapour pressure with height, as before.
    # Given: reuses temperature's own log/ψ_h ratio below (RHabove's approach).
    use_humidity_anchor = use_anchor && !isnothing(surface_relative_humidity)
    anchor_vapor_pressure = use_humidity_anchor ?
        vapour_pressure(vapour_pressure_equation, surface_temperature) * surface_relative_humidity : reference_vapor_pressure

    for i in 1:N_heights
        height_above_origin = heights[i] - displacement_height
        if height_above_origin < roughness_length
            # Goes negative below the displaced origin -- clamp so the
            # interpolation saturates at surface_temperature instead of extrapolating.
            height_fraction = clamp(ustrip(height_above_origin / roughness_length), 0.0, 1.0)
            wind_speed[i] = zero(friction_velocity)
            air_temperature[i] = surface_temperature + (roughness_height_temperature - surface_temperature) * height_fraction
            relative_humidity[i] = clamp(reference_vapor_pressure / vapour_pressure(vapour_pressure_equation, air_temperature[i]), 0.0, 1.0)
        elseif i == reference_index
            wind_speed[i] = reference_wind_speed
            air_temperature[i] = reference_temperature
            relative_humidity[i] = reference_humidity
        else
            height_stability_parameter = stability_parameter_at_height(stability_variable, height_above_origin, stability)
            wind_speed[i] = (friction_velocity / κ) * (log(height_above_origin / roughness_length) - momentum_stability_correction(stability_functions, height_stability_parameter) + ψ_m_roughness)
            # Temperature uses its own anchor (which may differ from roughness_length -- see above)
            temperature_log_term = log(height_above_origin / anchor_height_above_origin) + ψ_h_anchor - heat_stability_correction(stability_functions, height_stability_parameter)
            air_temperature[i] = anchor_temperature + (reference_temperature - anchor_temperature) * temperature_log_term / temperature_denominator
            vapor_pressure_h = use_humidity_anchor ?
                anchor_vapor_pressure + (reference_vapor_pressure - anchor_vapor_pressure) * temperature_log_term / temperature_denominator :
                reference_vapor_pressure
            relative_humidity[i] = clamp(vapor_pressure_h / vapour_pressure(vapour_pressure_equation, air_temperature[i]), 0.0, 1.0)
        end
    end

    return (;
        wind_speed,
        air_temperature,
        relative_humidity,
        convective_heat_flux=u"W/m^2"(convective_heat_flux),
        friction_velocity,
        obukhov_length=to_obukhov_length(stability_variable, stability),
        inverse_obukhov_length=to_inverse_obukhov_length(stability_variable, stability),
        stability_update=stability_result.stability_update,
    )
end

# No stability given: iterate it, warm-started from (and updating) the buffer's previous value
function surface_layer_stability(stability_variable, ::Nothing, obukhov_length_prev, args...; kw...)
    out = monin_obukhov_iteration(stability_variable, args...; max_iterations=30, tolerance=1e-2,
        initial_value=from_obukhov_length(stability_variable, obukhov_length_prev[]), kw...)
    obukhov_length_prev[] = to_obukhov_length(stability_variable, out.value)
    return out
end
# Stability given: one step at that value
surface_layer_stability(stability_variable, stability, obukhov_length_prev, args...; kw...) =
    (; value=stability, monin_obukhov_step(stability_variable, stability, args...; kw...)...)

function surface_fluxes(boundary_layer_model::MoninObukhov;
    surface_temperature, air_temperature, wind_speed, zenith_angle,
    roughness_height, reference_height, atmospheric_pressure, obukhov_length_prev=nothing,
)
    karman_constant = boundary_layer_model.karman_constant
    loglaw_ratio = reference_height / roughness_height
    loglaw_ratio > 1 || throw(DomainError(ustrip(loglaw_ratio),
        "surface_fluxes: reference_height ($reference_height) must exceed roughness_height " *
        "($roughness_height) for Monin-Obukhov log-law inversion"))
    loglaw_ratio < _MIN_LOGLAW_RATIO && @warn(
        "surface_fluxes: reference_height is close to roughness_height -- Monin-Obukhov " *
        "log-law inversion is unreliable in this roughness-sublayer regime",
        reference_height, roughness_height, loglaw_ratio, maxlog=10,
    )
    log_height_ratio = log(reference_height / roughness_height)
    temperature_difference = air_temperature - surface_temperature
    mean_temperature = (surface_temperature + air_temperature) / 2
    air_volumetric_heat_capacity = calc_air_volumetric_heat_capacity(mean_temperature)
    kinematic_viscosity = dry_air_properties(mean_temperature, atmospheric_pressure).kinematic_viscosity

    if air_temperature ≥ surface_temperature || zenith_angle ≥ 90°
        # Stable / nocturnal: neutral log-law.
        friction_velocity = max(
            calc_friction_velocity(; reference_wind_speed=wind_speed, log_height_ratio, karman_constant), boundary_layer_model.min_friction_velocity
        )
        convective_heat_flux = uconvert(u"W/m^2",
            calc_convection(; friction_velocity, log_height_ratio, temperature_difference, air_volumetric_heat_capacity, roughness_height, height=reference_height, karman_constant,
                kinematic_viscosity, thermal_roughness_model=boundary_layer_model.thermal_roughness_model))
    else
        # Unstable: iterate Obukhov length. 
        obukhov_length_factor = 6.003e-8u"cal*minute^2/cm^4"
        previous_obukhov_length = isnothing(obukhov_length_prev) ? -0.3u"m" : obukhov_length_prev[]
        initial_obukhov_length = previous_obukhov_length >= 0.0u"m" ? -0.3u"m" : previous_obukhov_length # Resetting to the default unstable guess?
        out = calc_obukhov_length(air_temperature, surface_temperature, wind_speed,
            roughness_height, reference_height, obukhov_length_factor, karman_constant, temperature_difference, air_volumetric_heat_capacity, kinematic_viscosity;
            initial_obukhov_length, boundary_layer_model.stability_functions, boundary_layer_model.min_friction_velocity, boundary_layer_model.thermal_roughness_model)
        friction_velocity = out.friction_velocity
        convective_heat_flux = out.convective_heat_flux
    end

    return (; convective_heat_flux, friction_velocity, temperature_difference)
end

function init_surface_fluxes!(boundary_layer_model::MoninObukhov, buffers, forcing, site, heights, soil_temperatures, i)
    t_next = ((i - 1) * 60)u"minute"
    (; air_temperature, wind_speed, zenith_angle) = interpolate_forcings(forcing, t_next)
    surface_temperature = soil_temperatures[1]
    if air_temperature < surface_temperature && zenith_angle < 90u"°"
        karman_constant = boundary_layer_model.karman_constant
        roughness_height = site.roughness_height
        reference_height = last(heights)
        temperature_difference = air_temperature - surface_temperature
        mean_temperature = (surface_temperature + air_temperature) / 2
        air_volumetric_heat_capacity = calc_air_volumetric_heat_capacity(mean_temperature)
        kinematic_viscosity = dry_air_properties(mean_temperature, site.atmospheric_pressure).kinematic_viscosity
        obukhov_length_factor = 6.003e-8u"cal*minute^2/cm^4"
        previous_obukhov_length = buffers.soil_energy_balance.obukhov_length_prev[]
        initial_obukhov_length = previous_obukhov_length >= 0.0u"m" ? -0.3u"m" : previous_obukhov_length
        out = calc_obukhov_length(air_temperature, surface_temperature, wind_speed,
            roughness_height, reference_height, obukhov_length_factor, karman_constant, temperature_difference, air_volumetric_heat_capacity, kinematic_viscosity;
            initial_obukhov_length, boundary_layer_model.stability_functions, boundary_layer_model.min_friction_velocity, boundary_layer_model.thermal_roughness_model)
        buffers.soil_energy_balance.obukhov_length_prev[] = out.obukhov_length end
end

"""
    calc_air_volumetric_heat_capacity(mean_temperature)

Compute the volumetric heat capacity of air (ρ·cₚ) as a function of mean temperature.

# Arguments
- `mean_temperature`: Mean air temperature (`Unitful.Temperature`), in Kelvin.

# Returns
- Volumetric heat capacity (`cal / (cm³·K)`).

This is a simplified empirical regression based only on temperature,
without accounting for moisture or elevation effects.
"""
function calc_air_volumetric_heat_capacity(mean_temperature)
    return u"J/m^3/K"(u"(cal*g)/(g*cm^3*K)" * (0.08472 / ustrip(u"K", mean_temperature)))
end

"""
    calc_air_volumetric_heat_capacity(mean_temperature, elevation, relative_humidity, atmospheric_pressure)

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
function calc_air_volumetric_heat_capacity(mean_temperature, elevation, relative_humidity, atmospheric_pressure; vapour_pressure_equation=GoffGratch())
    dry_air_out = dry_air_properties(u"K"(mean_temperature), atmospheric_pressure)
    wet_air_out = wet_air_properties(u"K"(mean_temperature), relative_humidity, atmospheric_pressure; vapour_pressure_equation)
    air_density = dry_air_out.density
    air_heat_capacity = wet_air_out.specific_heat
    return air_density * air_heat_capacity
end

"""
    calc_friction_velocity(; reference_wind_speed, log_height_ratio, karman_constant=0.4)

Compute the friction velocity (u*) from a reference wind speed using the
logarithmic wind profile.

# Arguments
- `reference_wind_speed::Quantity{<:Real,𝐋/𝐓}`: Wind speed at the reference height (e.g. `m/s`, `cm/min`).
- `log_height_ratio::Real`: Precomputed log height ratio, typically `log(height / roughness_height + 1.0)`.
- `karman_constant::Real`: von Kármán constant (default = 0.4).

# Returns
- Friction velocity `friction_velocity::Quantity{<:Real,𝐋/𝐓}`.

# See also
[`calc_convection`](@ref)
"""
function calc_friction_velocity(; reference_wind_speed, log_height_ratio, karman_constant=0.4)
    return karman_constant * reference_wind_speed / log_height_ratio
end

"""
    calc_convection(; friction_velocity, log_height_ratio, temperature_difference, air_volumetric_heat_capacity, roughness_height, height, karman_constant, kinematic_viscosity, thermal_roughness_model=SublayerStantonRoughness())

Calculate the convective heat flux (sensible heat exchange between surface and air).

# Arguments
- `friction_velocity::Quantity{<:Real,𝐋/𝐓}`: Friction velocity (e.g. `m/s`, `cm/min`).
- `log_height_ratio::Real`: Precomputed logarithmic height ratio, typically `log(height / roughness_height + 1.0)`.
- `temperature_difference::Quantity{<:Real,Θ}`: Temperature difference between reference air and surface (Kelvin).
- `air_volumetric_heat_capacity::Quantity{<:Real,(𝐌*𝐋^-1*𝐓^-2)}`: Volumetric heat capacity of air (e.g. `J/m³/K`, `cal/cm³/K`).
- `roughness_height::Quantity{<:Real,𝐋}`: Surface roughness length.

# Returns
- Convective heat flux as `Quantity{<:Real,(𝐌*𝐓^-3)}` (e.g. `W/m²`, `cal/min/cm²`).

Uses bulk and sublayer Stanton numbers to account for turbulence near the surface.

# See also
[`calc_friction_velocity`](@ref), [`sublayer_stanton`](@ref), [`bulk_stanton`](@ref), [`convective_flux`](@ref)
"""
function calc_convection(; friction_velocity, log_height_ratio, temperature_difference, air_volumetric_heat_capacity, roughness_height, height, karman_constant, kinematic_viscosity, thermal_roughness_model=SublayerStantonRoughness())
    return _stable_convective_heat_flux(thermal_roughness_model, air_volumetric_heat_capacity, temperature_difference, friction_velocity, karman_constant, roughness_height, height, log_height_ratio, kinematic_viscosity)
end

"""
    convective_flux(air_volumetric_heat_capacity, temperature_difference, friction_velocity, bulk_stanton_number, sublayer_stanton_number)

Compute convective heat flux given bulk and sublayer Stanton numbers.
"""
@inline function convective_flux(air_volumetric_heat_capacity, temperature_difference, friction_velocity, bulk_stanton_number, sublayer_stanton_number)
    return air_volumetric_heat_capacity * temperature_difference * friction_velocity * bulk_stanton_number / (1 + bulk_stanton_number / sublayer_stanton_number)
end

# Scalar-roughness-ratio resistance, ln(z/zh) - ψ_h(z/L) + ψ_h(zh/L)
function _scalar_roughness_resistance(ratio, karman_constant, friction_velocity, roughness_height, height, heat_correction_at_height, heat_correction_at_scalar_roughness)
    scalar_roughness_height = ratio * roughness_height
    return (log(height / scalar_roughness_height) - heat_correction_at_height + heat_correction_at_scalar_roughness) / (karman_constant * friction_velocity)
end

# calc_convection's stable/neutral direct branch: no Obukhov iteration, ψ_h implicitly 0.
_stable_convective_heat_flux(::SublayerStantonRoughness, air_volumetric_heat_capacity, temperature_difference, friction_velocity, karman_constant, roughness_height, height, log_height_ratio, kinematic_viscosity) =
    convective_flux(air_volumetric_heat_capacity, temperature_difference, friction_velocity, bulk_stanton(log_height_ratio), sublayer_stanton(roughness_height, friction_velocity, kinematic_viscosity))
function _stable_convective_heat_flux(thermal_roughness_model::ScalarRoughnessRatio, air_volumetric_heat_capacity, temperature_difference, friction_velocity, karman_constant, roughness_height, height, log_height_ratio, kinematic_viscosity)
    resistance = _scalar_roughness_resistance(thermal_roughness_model.ratio, karman_constant, friction_velocity, roughness_height, height, 0.0, 0.0)
    return air_volumetric_heat_capacity * temperature_difference / resistance
end

# calc_obukhov_length's iterative branch. Also returns roughness_height_temperature
# (the below-z0 profile's own boundary value) since SublayerStantonRoughness derives
# it from the same bulk/sublayer split as the flux itself.
function _unstable_convective_heat_flux(::SublayerStantonRoughness, stability_functions, air_volumetric_heat_capacity, temperature_difference, friction_velocity, karman_constant, roughness_height, height, log_ratio_corrected, ψ_h, stability_parameter, kinematic_viscosity, reference_temperature, surface_temperature)
    sublayer_stanton_number = sublayer_stanton(roughness_height, friction_velocity, kinematic_viscosity)
    bulk_stanton_number = bulk_stanton(log_ratio_corrected, stability_parameter)
    convective_heat_flux = convective_flux(air_volumetric_heat_capacity, temperature_difference, friction_velocity, bulk_stanton_number, sublayer_stanton_number)
    roughness_height_temperature = (reference_temperature * bulk_stanton_number + surface_temperature * sublayer_stanton_number) / (bulk_stanton_number + sublayer_stanton_number)
    return (; convective_heat_flux, roughness_height_temperature)
end
function _unstable_convective_heat_flux(thermal_roughness_model::ScalarRoughnessRatio, stability_functions, air_volumetric_heat_capacity, temperature_difference, friction_velocity, karman_constant, roughness_height, height, log_ratio_corrected, ψ_h, stability_parameter, kinematic_viscosity, reference_temperature, surface_temperature)
    (; ratio) = thermal_roughness_model
    heat_correction_at_scalar_roughness = heat_stability_correction(stability_functions, stability_parameter * ratio * roughness_height / height)
    resistance = _scalar_roughness_resistance(ratio, karman_constant, friction_velocity, roughness_height, height, ψ_h, heat_correction_at_scalar_roughness)
    return (; convective_heat_flux=air_volumetric_heat_capacity * temperature_difference / resistance, roughness_height_temperature=surface_temperature)
end


"""
    calc_heat_transfer_coefficient(convective_heat_flux, temperature_difference)

Heat transfer coefficient from convective heat flux and temperature difference.
Minimum 0.5 W/m²/K.
"""
calc_heat_transfer_coefficient(convective_heat_flux, temperature_difference) =
    iszero(temperature_difference) ? 0.5u"W/m^2/K" : max(abs(convective_heat_flux / temperature_difference), 0.5u"W/m^2/K")

"""
    LEWIS_HEAT_TO_MASS_RATIO

(Pr/Sc)^0.666 = (0.71/0.60)^0.666, the Lewis-relation factor converting a
heat transfer coefficient to a mass (vapor) transfer coefficient.
"""
const LEWIS_HEAT_TO_MASS_RATIO = (0.71 / 0.60)^0.666

"""
    calc_mass_transfer_coefficient(heat_transfer_coefficient, air_specific_heat, air_density)

Mass transfer coefficient from heat transfer coefficient via the Lewis relation.
"""
calc_mass_transfer_coefficient(heat_transfer_coefficient, air_specific_heat, air_density) =
    (heat_transfer_coefficient / (air_specific_heat * air_density)) * LEWIS_HEAT_TO_MASS_RATIO

"""
    sublayer_stanton(roughness_height, friction_velocity, kinematic_viscosity)

Stanton number for the roughness sublayer immediately above the surface,
`0.62·Re*^-0.45` with roughness Reynolds number `Re* = roughness_height·u*/ν` (form
resembles Owen & Thomson (1963)) Not suitable for a canopy.
"""
function sublayer_stanton(roughness_height, friction_velocity, kinematic_viscosity)
    roughness_reynolds_number = uconvert(NoUnits, roughness_height * friction_velocity / kinematic_viscosity)
    return 0.62 / roughness_reynolds_number^0.45
end

"""
    bulk_stanton(log_height_ratio)

Compute the bulk Stanton number for stable conditions.
"""
@inline function bulk_stanton(log_height_ratio)
    return 0.64 / log_height_ratio
end

"""
    bulk_stanton(log_height_ratio, stability_parameter)

Compute the bulk Stanton number for unstable conditions.
"""
@inline function bulk_stanton(log_height_ratio, stability_parameter)
    return (0.64 / log_height_ratio) * (1 - 0.1 * stability_parameter)
end


"""
    AbstractStabilityVariable

How atmospheric stability is represented: [`ObukhovLength`](@ref) or
[`InverseObukhovLength`](@ref).
"""
abstract type AbstractStabilityVariable end

"""
    ObukhovLength()

Stability as the Obukhov length (m), infinite at neutral.
"""
struct ObukhovLength <: AbstractStabilityVariable end

"""
    InverseObukhovLength()

Stability as the inverse Obukhov length (1/m). It is zero rather than
infinite at neutral, so it and its derivatives stay finite through neutral
conditions.
"""
struct InverseObukhovLength <: AbstractStabilityVariable end

"""
    stability_parameter_at_height(stability_variable, height, value)

The Monin–Obukhov stability parameter `height / obukhov_length`, from either
stability variable. Zero at neutral.
"""
stability_parameter_at_height(::ObukhovLength, height, obukhov_length) = uconvert(NoUnits, height / obukhov_length)
stability_parameter_at_height(::InverseObukhovLength, height, inverse_obukhov_length) = uconvert(NoUnits, height * inverse_obukhov_length)

# Obukhov relation, obukhov_length = obukhov_length_factor·friction_velocity³/convective_heat_flux
stability_from_fluxes(::ObukhovLength, obukhov_length_factor, friction_velocity, convective_heat_flux) =
    u"m"(obukhov_length_factor * friction_velocity^3 / convective_heat_flux)
stability_from_fluxes(::InverseObukhovLength, obukhov_length_factor, friction_velocity, convective_heat_flux) =
    u"m^-1"(convective_heat_flux / (obukhov_length_factor * friction_velocity^3))

# Convergence measure: relative change in the Obukhov length, or absolute
# change in the stability parameter for the inverse
stability_change(::ObukhovLength, height, new, old) = abs((new - old) / old)
stability_change(::InverseObukhovLength, height, new, old) = abs(uconvert(NoUnits, height * (new - old)))

# Conversions to and from the Obukhov length kept in warm-start buffers
from_obukhov_length(::ObukhovLength, obukhov_length) = obukhov_length
from_obukhov_length(::InverseObukhovLength, obukhov_length) = inv(obukhov_length)
to_obukhov_length(::ObukhovLength, obukhov_length) = obukhov_length
to_obukhov_length(::InverseObukhovLength, inverse_obukhov_length) = inv(inverse_obukhov_length)
to_inverse_obukhov_length(::ObukhovLength, obukhov_length) = inv(obukhov_length)
to_inverse_obukhov_length(::InverseObukhovLength, inverse_obukhov_length) = inverse_obukhov_length


"""
    paulson_x(stability_parameter, γ)

The Paulson (1970) variable `x = (1 - γ z / L)^(1/4)`, used as the argument to `momentum_stability_correction`
and `heat_stability_correction`. It is the reciprocal of the Businger–Dyer stability function
`φₘ = (1 - γ z / L)^(-1/4)`.

# Arguments
- `stability_parameter`: `z/L`, negative (unstable), see [`stability_parameter_at_height`](@ref).
- `γ`: Empirical constant (dimensionless, often ≈16).

# References
- Paulson, C. A. (1970). The mathematical representation of wind speed and
  temperature profiles in the unstable atmospheric surface layer.
  *Journal of Applied Meteorology*, 9(6), 857–861.
- Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971).
  Flux–profile relationships in the atmospheric surface layer.
  *Journal of the Atmospheric Sciences*, 28(2), 181–189.
- Dyer, A. J. (1974). A review of flux–profile relationships.
  *Boundary-Layer Meteorology*, 7(3), 363–372.
"""
@inline function paulson_x(stability_parameter, γ)
    # sqrt is faster than ^1/4
    return sqrt(sqrt(1.0 - γ * stability_parameter))
end


"""
    heat_gradient_stability_factor(stability_parameter, γ)
    heat_gradient_stability_factor(stability_functions::AbstractStabilityFunctions, stability_parameter)

Dimensionless heat gradient Φ_h, the eddy diffusivity's stability divisor.
2-arg form: the unstable Businger–Dyer expression, `x⁻²` of `paulson_x`'s
Paulson x-substitution (unlike `paulson_x` itself, which returns `x`, not
the momentum φₘ = x⁻¹). With `stability_functions`: either regime.

# References
- Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971).
  Flux–profile relationships in the atmospheric surface layer.
  *Journal of the Atmospheric Sciences*, 28(2), 181–189.
"""
function heat_gradient_stability_factor(stability_parameter, γ)
    return paulson_x(stability_parameter, γ)^(-2)
end


"""
    momentum_stability_correction(x)
    momentum_stability_correction(stability_functions::AbstractStabilityFunctions, stability_parameter)

Stability correction ψₘ for momentum. 1-arg form: the unstable Businger–Dyer
expression, `ψₘ(x) = 2 ln((1+x)/2) + ln((1+x²)/2) - 2 atan(x) + π/2`, `x`
typically `(1 - γ z / L)^(1/4)`. With `stability_functions`: either regime,
continuous at neutral.

# References
- Businger et al. (1971).
- Dyer (1974).
"""
@inline function momentum_stability_correction(x)
    return 2.0 * log((1.0 + x) / 2.0) + log((1.0 + x^2) / 2.0) - 2.0 * atan(x) + π / 2.0
end


"""
    heat_stability_correction(x)
    heat_stability_correction(stability_functions::AbstractStabilityFunctions, stability_parameter)

Stability correction ψ_h for heat and vapour. 1-arg form: the unstable
Businger–Dyer expression, `ψ_h(x) = 2 ln((1+x²)/2)`, `x` typically
`(1 - γ z / L)^(1/4)`. With `stability_functions`: either regime,
continuous at neutral.

# References
- Businger et al. (1971).
- Dyer (1974).
"""
@inline function heat_stability_correction(x)
    return 2.0 * log((1.0 + x^2.0) / 2.0)
end

# NicheMapR's bound on the stability parameter
bounded_stability_parameter(stability_parameter) = clamp(stability_parameter, -10.0, 10.0)

function momentum_stability_correction(stability_functions::BusingerDyer, stability_parameter)
    stability_parameter = bounded_stability_parameter(stability_parameter)
    stability_parameter < 0 && return momentum_stability_correction(paulson_x(stability_parameter, stability_functions.dyer_constant))
    return clamp(-stability_functions.stable_beta * stability_parameter, -4.0, 3.0)
end

function heat_stability_correction(stability_functions::BusingerDyer, stability_parameter)
    stability_parameter = bounded_stability_parameter(stability_parameter)
    stability_parameter < 0 && return heat_stability_correction(paulson_x(stability_parameter, stability_functions.dyer_constant))
    return clamp(-(stability_functions.stable_beta * stability_parameter) / stability_functions.turbulent_prandtl_number, -4.0, 3.0)
end

function heat_gradient_stability_factor(stability_functions::BusingerDyer, stability_parameter)
    stability_parameter = bounded_stability_parameter(stability_parameter)
    stability_parameter < 0 && return heat_gradient_stability_factor(stability_parameter, stability_functions.dyer_constant)
    (; stable_heat_gradient_coefficient, min_stable_heat_gradient, max_stable_heat_gradient) = stability_functions
    return clamp(1.0 + stable_heat_gradient_coefficient * stability_parameter / (1.0 + stability_parameter), min_stable_heat_gradient, max_stable_heat_gradient)
end

function momentum_stability_correction(stability_functions::ChengBrutsaert, stability_parameter)
    stability_parameter < 0 && return momentum_stability_correction(paulson_x(stability_parameter, stability_functions.momentum_dyer_constant))
    (; momentum_scale, momentum_exponent) = stability_functions
    return -momentum_scale * log(stability_parameter + (1 + stability_parameter^momentum_exponent)^(1 / momentum_exponent))
end

function heat_stability_correction(stability_functions::ChengBrutsaert, stability_parameter)
    stability_parameter < 0 && return heat_stability_correction(paulson_x(stability_parameter, stability_functions.heat_dyer_constant))
    (; heat_scale, heat_exponent) = stability_functions
    return -heat_scale * log(stability_parameter + (1 + stability_parameter^heat_exponent)^(1 / heat_exponent))
end

function heat_gradient_stability_factor(stability_functions::ChengBrutsaert, stability_parameter)
    stability_parameter < 0 && return heat_gradient_stability_factor(stability_parameter, stability_functions.heat_dyer_constant)
    # Φ_h = 1 - ζ·dψ_h/dζ, for ψ_h = -c·ln g, g = ζ + (1 + ζᵈ)^(1/d), with ζ the stability parameter
    (; heat_scale, heat_exponent) = stability_functions
    q = 1 + stability_parameter^heat_exponent
    g = stability_parameter + q^(1 / heat_exponent)
    return 1 + heat_scale * stability_parameter * (1 + stability_parameter^(heat_exponent - 1) * q^(1 / heat_exponent - 1)) / g
end


"""
    monin_obukhov_step(stability_variable, value, reference_temperature, surface_temperature, reference_wind_speed, roughness_length, reference_height_above_origin, obukhov_length_factor, karman_constant, temperature_difference,
                       air_volumetric_heat_capacity, kinematic_viscosity; stability_functions=ChengBrutsaert(),
                       min_friction_velocity=0.02u"m/s", thermal_roughness_model=SublayerStantonRoughness())

One evaluation of the Monin–Obukhov relations at stability `value` (as
`stability_variable`, [`ObukhovLength`](@ref) or [`InverseObukhovLength`](@ref)):
friction velocity and convective heat flux from the stability corrections,
then `stability_update`, the stability those fluxes imply through the
Obukhov relation. `value` is self-consistent when `stability_update == value`.
`kinematic_viscosity` feeds `sublayer_stanton`'s roughness Reynolds number.
Works for either sign of `temperature_difference` -- `convective_flux`'s sign already tracks
`temperature_difference`'s sign (negative temperature_difference, surface warmer than reference → negative heat flux
→ unstable; positive temperature_difference → stable).
"""
function monin_obukhov_step(stability_variable, value, reference_temperature, surface_temperature, reference_wind_speed, roughness_length, reference_height_above_origin, obukhov_length_factor, karman_constant, temperature_difference,
    air_volumetric_heat_capacity, kinematic_viscosity;
    stability_functions=ChengBrutsaert(),
    min_friction_velocity=0.02u"m/s", thermal_roughness_model=SublayerStantonRoughness(),
)
    reference_stability_parameter = stability_parameter_at_height(stability_variable, reference_height_above_origin, value)
    ψ_m = momentum_stability_correction(stability_functions, reference_stability_parameter)
    ψ_h = heat_stability_correction(stability_functions, reference_stability_parameter)
    log_ratio_corrected = log(reference_height_above_origin / roughness_length) - ψ_m +
        momentum_stability_correction(stability_functions, reference_stability_parameter * roughness_length / reference_height_above_origin)
    friction_velocity = max(karman_constant * reference_wind_speed / log_ratio_corrected, min_friction_velocity)
    convective_heat_flux, roughness_height_temperature = _unstable_convective_heat_flux(thermal_roughness_model,
        stability_functions, air_volumetric_heat_capacity, temperature_difference, friction_velocity, karman_constant, roughness_length, reference_height_above_origin, log_ratio_corrected, ψ_h, reference_stability_parameter, kinematic_viscosity,
        reference_temperature, surface_temperature)
    stability_update = stability_from_fluxes(stability_variable, obukhov_length_factor, friction_velocity, convective_heat_flux)
    return (;
        stability_update,
        friction_velocity,
        heat_correction=ψ_h,
        convective_heat_flux=u"W/m^2"(convective_heat_flux),
        roughness_height_temperature,
    )
end

"""
    monin_obukhov_iteration(stability_variable, args...; initial_value, max_iterations=50, tolerance=1e-2, relaxation=0.5, kw...)

Iterate [`monin_obukhov_step`](@ref) (same `args` and `kw`) from
`initial_value` until the stability stops changing, under-relaxing each
update by `relaxation`. Returns the last step's fluxes and `value`, the
stability that step was evaluated at.
"""
function monin_obukhov_iteration(stability_variable, reference_temperature, surface_temperature, reference_wind_speed, roughness_length, reference_height_above_origin, obukhov_length_factor, karman_constant, temperature_difference,
    air_volumetric_heat_capacity, kinematic_viscosity; initial_value, max_iterations=50, tolerance=1e-2, relaxation=0.5, kw...,
)
    value = initial_value
    iteration = 1
    while true
        step = monin_obukhov_step(stability_variable, value, reference_temperature, surface_temperature, reference_wind_speed, roughness_length, reference_height_above_origin, obukhov_length_factor, karman_constant,
            temperature_difference, air_volumetric_heat_capacity, kinematic_viscosity; kw...)
        change = stability_change(stability_variable, reference_height_above_origin, step.stability_update, value)
        (change <= tolerance || iteration >= max_iterations) && return (; value, step...)
        value = relaxation * step.stability_update + (1.0 - relaxation) * value
        iteration += 1
    end
end

"""
    calc_obukhov_length(reference_temperature, surface_temperature, reference_wind_speed, roughness_length, reference_height_above_origin, obukhov_length_factor, karman_constant, temperature_difference, air_volumetric_heat_capacity, kinematic_viscosity;
                         initial_obukhov_length=-0.3u"m", kw...)

[`monin_obukhov_iteration`](@ref) on the Obukhov length. Converges to `±Inf`
(neutral) as `temperature_difference → 0`.
"""
function calc_obukhov_length(args...; initial_obukhov_length=-0.3u"m", kw...)
    (; value, friction_velocity, heat_correction, convective_heat_flux, roughness_height_temperature) =
        monin_obukhov_iteration(ObukhovLength(), args...; initial_value=initial_obukhov_length, kw...)
    return (; obukhov_length=value, friction_velocity, heat_correction, convective_heat_flux, roughness_height_temperature)
end

"""
    inverse_obukhov_length_limits(stability_limit, stability_functions, height, roughness_length)

The range of inverse Obukhov length `stability_limit` allows, as `(; unstable, stable)`.

[`NoStabilityLimit`](@ref) is unbounded. [`MomentumCorrectionLimit`](@ref) gives the values
nearest to neutral where the momentum stability correction between `roughness_length` and
`height` (above the displacement height), `ψ_m(z₀/L) - ψ_m(z/L)`, reaches `fraction` of the log
term `ln(z/z₀)`, or an infinite value where it never does. Within these limits the friction
velocity changes with stability by at most a factor `1 ± fraction`.
"""
inverse_obukhov_length_limits(::NoStabilityLimit, stability_functions, height, roughness_length) =
    (; unstable=-Inf * u"m^-1", stable=Inf * u"m^-1")
function inverse_obukhov_length_limits(limit::MomentumCorrectionLimit, stability_functions, height, roughness_length)
    (; fraction) = limit
    target = fraction * log(height / roughness_length)
    correction(inverse_length) = momentum_stability_correction(stability_functions, roughness_length * inverse_length) -
        momentum_stability_correction(stability_functions, height * inverse_length)
    nearest_to_neutral(sign) = begin
        # The first doubling of |1/L| past the target, then bisection inside it
        previous = 0.0u"m^-1"
        for exponent in -30:40
            inverse_length = sign * 2.0^exponent * u"m^-1"
            if sign * correction(inverse_length) >= target
                below_target, above_target = previous, inverse_length
                for _ in 1:60
                    middle = (below_target + above_target) / 2
                    sign * correction(middle) >= target ? (above_target = middle) : (below_target = middle)
                end
                return above_target
            end
            previous = inverse_length
        end
        return sign * Inf * u"m^-1"
    end
    return (; unstable=nearest_to_neutral(-1.0), stable=nearest_to_neutral(1.0))
end

# The inverse Obukhov length within `limits`: unchanged near neutral, approaching the limit
# on its side smoothly, `limit·tanh(1/L / limit)`
function limit_inverse_obukhov_length(inverse_obukhov_length, limits)
    limit = inverse_obukhov_length < zero(inverse_obukhov_length) ? limits.unstable : limits.stable
    isinf(limit) && return inverse_obukhov_length
    return limit * tanh(inverse_obukhov_length / limit)
end
