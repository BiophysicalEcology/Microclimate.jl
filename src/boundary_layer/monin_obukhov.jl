"""
    AbstractUnstableStability

Monin–Obukhov stability corrections in unstable air, `ζ = z/L < 0`: variants implement
[`momentum_stability_correction`](@ref), [`heat_stability_correction`](@ref) and
[`heat_gradient_stability_factor`](@ref) of `ζ`.
"""
abstract type AbstractUnstableStability end

"""
    AbstractStableStability

Monin–Obukhov stability corrections in stable air, `ζ = z/L ≥ 0`: variants implement
[`momentum_stability_correction`](@ref), [`heat_stability_correction`](@ref) and
[`heat_gradient_stability_factor`](@ref) of `ζ`.
"""
abstract type AbstractStableStability end

"""
    MoninObukhov(; karman_constant=0.4, unstable_stability=PaulsonDyer(), stable_stability=ChengBrutsaert(),
                   min_friction_velocity=0.1u"m/s", thermal_roughness_model=SublayerStantonRoughness(),
                   stability_limit=NoStabilityLimit())

Monin–Obukhov similarity theory boundary-layer formulation. Stability corrections come from
`unstable_stability` ([`AbstractUnstableStability`](@ref)) in unstable air and
`stable_stability` ([`AbstractStableStability`](@ref)) in stable air. The default pairs
Paulson's integrated Businger–Dyer forms with Cheng & Brutsaert (2005);
`PaulsonDyer()` with [`BeljaarsHoltslag`](@ref) is ECMWF's pair, and
[`NicheMapRUnstable`](@ref) with [`NicheMapRStable`](@ref) is NicheMapR's.
`stability_limit` bounds the canopy-top stability ([`NoStabilityLimit`](@ref) or
[`MomentumCorrectionLimit`](@ref)).

`thermal_roughness_model` sets heat roughness from momentum roughness.
"""
@kwdef struct MoninObukhov{KC,US<:AbstractUnstableStability,SS<:AbstractStableStability,MINUF,TRM,SL} <: AbstractBoundaryLayerModel
    karman_constant::KC = 0.4
    unstable_stability::US = PaulsonDyer()
    stable_stability::SS = ChengBrutsaert()
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
velocity by at most a factor `1/(1 + fraction)` in stable air and `1/(1 - fraction)` in
unstable air.

On calm nights a canopy can cool faster than stability reduces its heat flux, and then no
stability agrees with that flux. This limit gives one. Stability approaches the limit smoothly.
"""
@kwdef struct MomentumCorrectionLimit{F} <: AbstractStabilityLimit
    fraction::F = 0.9
end

"""
    PaulsonDyer(; momentum_dyer_constant=16.0, heat_dyer_constant=16.0)

Paulson's (1970) integrals of the Businger–Dyer flux-gradient relations,
`φ_m = (1 - γ_m·ζ)^(-1/4)` and `Φ_h = (1 - γ_h·ζ)^(-1/2)`, with `γ_m` the
`momentum_dyer_constant` and `γ_h` the `heat_dyer_constant` (Dyer 1974).

# References
- Paulson, C. A. (1970). The mathematical representation of wind speed and
  temperature profiles in the unstable atmospheric surface layer.
  *Journal of Applied Meteorology*, 9(6), 857–861.
- Dyer, A. J. (1974). A review of flux–profile relationships.
  *Boundary-Layer Meteorology*, 7(3), 363–372.
"""
@kwdef struct PaulsonDyer{DM,DH} <: AbstractUnstableStability
    momentum_dyer_constant::DM = 16.0
    heat_dyer_constant::DH = 16.0
end

"""
    ChengBrutsaert(; momentum_scale=6.1, momentum_exponent=2.5, heat_scale=5.3, heat_exponent=1.1)

Cheng & Brutsaert's (2005) stable corrections, fitted to CASES-99,
`ψ_m = -a·ln(ζ + (1 + ζᵇ)^(1/b))` with `a = momentum_scale` and `b = momentum_exponent`,
`ψ_h = -c·ln(ζ + (1 + ζᵈ)^(1/d))` with `c = heat_scale` and `d = heat_exponent`. They
level off in very stable air, so every stable state has a solution. Paired with an
unstable form they are continuous at neutral, but their slopes there differ from it.

# References
Cheng, Y., & Brutsaert, W. (2005). Flux-profile relationships for wind speed
and temperature in the stable atmospheric boundary layer.
*Boundary-Layer Meteorology*, 114, 519–538.
"""
@kwdef struct ChengBrutsaert{A,B,C,D} <: AbstractStableStability
    momentum_scale::A = 6.1
    momentum_exponent::B = 2.5
    heat_scale::C = 5.3
    heat_exponent::D = 1.1
end

"""
    BusingerLinear(; coefficient=5.0)

The log-linear stable corrections, `ψ_m = ψ_h = -coefficient·ζ` and
`Φ_h = 1 + coefficient·ζ` (Dyer 1974). They overcorrect past `ζ ≈ 1`, where turbulence
under them collapses at a critical bulk Richardson number of `1/coefficient`.

# References
- Dyer, A. J. (1974). A review of flux–profile relationships.
  *Boundary-Layer Meteorology*, 7(3), 363–372.
"""
@kwdef struct BusingerLinear{C} <: AbstractStableStability
    coefficient::C = 5.0
end

"""
    BeljaarsHoltslag(; a=1.0, b=2/3, c=5.0, d=0.35)

Beljaars & Holtslag's (1991) stable corrections, as ECMWF's IFS uses them,
`ψ_m = -(a·ζ + b·(ζ - c/d)·exp(-d·ζ) + b·c/d)` and
`ψ_h = -((1 + 2a·ζ/3)^(3/2) + b·(ζ - c/d)·exp(-d·ζ) + b·c/d - 1)`. Like
[`ChengBrutsaert`](@ref) they weaken the correction at strong stability.

# References
Beljaars, A. C. M., & Holtslag, A. A. M. (1991). Flux parameterization over land
surfaces for atmospheric models. *Journal of Applied Meteorology*, 30(3), 327–341.
"""
@kwdef struct BeljaarsHoltslag{A,B,C,D} <: AbstractStableStability
    a::A = 1.0
    b::B = 2 / 3
    c::C = 5.0
    d::D = 0.35
end

"""
    NicheMapRUnstable(; dyer_constant=16.0)

NicheMapR's unstable corrections: [`PaulsonDyer`](@ref) with one `dyer_constant` for
momentum and heat, and `ζ` bounded at -10.
"""
@kwdef struct NicheMapRUnstable{DC} <: AbstractUnstableStability
    dyer_constant::DC = 16.0
end

"""
    NicheMapRStable(; stable_beta=4.7, turbulent_prandtl_number=0.74,
        stable_heat_gradient_coefficient=6.0, min_stable_heat_gradient=0.5, max_stable_heat_gradient=1.5)

NicheMapR's stable corrections: Businger et al.'s (1971) linear form,
`ψ_m = -stable_beta·ζ` and `ψ_h = -(stable_beta/turbulent_prandtl_number)·ζ`, clamped to
`[-4, 3]`, and a separate saturating diffusivity multiplier
`Φ_h = clamp(1 + stable_heat_gradient_coefficient·ζ/(1+ζ), min_stable_heat_gradient, max_stable_heat_gradient)`, with
`ζ` bounded at 10. The clamps put kinks in every function of stability.

# References
Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971).
Flux–profile relationships in the atmospheric surface layer.
*Journal of the Atmospheric Sciences*, 28(2), 181–189.
"""
@kwdef struct NicheMapRStable{SB,PR,PHC,MINPH,MAXPH} <: AbstractStableStability
    stable_beta::SB = 4.7
    turbulent_prandtl_number::PR = 0.74
    stable_heat_gradient_coefficient::PHC = 6.0
    min_stable_heat_gradient::MINPH = 0.5
    max_stable_heat_gradient::MAXPH = 1.5
end

"""
    AbstractThermalRoughnessModel

How heat roughness relates to momentum roughness.
"""
abstract type AbstractThermalRoughnessModel end

"""
    SublayerStantonRoughness()

Default, for bare ground: [`sublayer_stanton`](@ref) in series with the log-law Stanton
number ([`convective_flux`](@ref)).
"""
struct SublayerStantonRoughness <: AbstractThermalRoughnessModel end

"""
    ScalarRoughnessRatio(; ratio=0.2)

Heat roughness is a fixed `ratio` of momentum roughness.
"""
@kwdef struct ScalarRoughnessRatio{R} <: AbstractThermalRoughnessModel
    ratio::R = 0.2
end

const _MIN_LOGLAW_RATIO = 5.0 # avoid the log-law inversion becoming ill-conditioned

# The unstable first guess, before any previous value exists to start from
const INITIAL_OBUKHOV_LENGTH = -0.3u"m"

# ρ·cₚ·T/(κ·g), the factor in the Obukhov relation L = factor·u*³/H. A constant, as in
# the original Fortran. TODO calculate it from the air properties.
const OBUKHOV_LENGTH_FACTOR = u"J*s^2/m^4"(6.003e-8u"cal*minute^2/cm^4")

function allocate_profile(::MoninObukhov, heights)
    wind_speed = similar(heights, typeof(0.0u"m/s")) # output wind speeds
    air_temperature = similar(heights, typeof(0.0u"K")) # output temperatures, need to do this otherwise get InexactError
    relative_humidity = similar(heights, Float64) # output relative humidities
    obukhov_length_prev = Ref(INITIAL_OBUKHOV_LENGTH)  # warm-start across timesteps
    warned_below_roughness = Ref(false)  # emit warning at most once per buffer
    warned_above_reference = Ref(false)
    return (; heights, air_temperature, wind_speed, relative_humidity,
              obukhov_length_prev, warned_below_roughness, warned_above_reference)
end

"""
    atmospheric_surface_profile(boundary_layer_model::MoninObukhov; heights=DEFAULT_HEIGHTS, reference_height=last(heights), kwargs...)

Compute vertical profiles of wind speed, air temperature, and relative humidity
using Monin–Obukhov similarity theory (MOST). Profiles are returned at all requested
`heights`; heights below the roughness length return zero wind and linearly interpolated
temperature; heights above `reference_height` are extrapolated via the log-law.

# Keyword Arguments

- `heights`: Heights at which profiles are returned (default: `DEFAULT_HEIGHTS`).
- `reference_height`: Measurement height for met inputs (default: `last(heights)`).
- `surface_temperature`: Soil/surface temperature.
- `roughness_length`: roughness length of the surface, e.g. the ground roughness length.
- `environment_instant`: Named tuple with `reference_temperature`, `reference_wind_speed`,
  `reference_humidity` and `atmospheric_pressure`.

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
    roughness_length = 0.004u"m",
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
    atmospheric_surface_profile!(boundary_layer_model::MoninObukhov, buffers; roughness_length, environment_instant,
                                  surface_temperature, vapour_pressure_equation=GoffGratch(),
                                  reference_height=last(buffers.heights),
                                  displacement_height=0.0u"m",
                                  thermal_roughness_model=boundary_layer_model.thermal_roughness_model,
                                  temperature_anchor_height=nothing, surface_relative_humidity=nothing,
                                  stability_variable=ObukhovLength(), stability=nothing)

`displacement_height` and `roughness_length` move the log profile's origin, e.g. for a
canopy. Heights are still measured from the ground.

`temperature_anchor_height` fixes `surface_temperature`, and `surface_relative_humidity`
if given, at that height instead of at the roughness length.

`stability_variable` ([`ObukhovLength`](@ref) or [`InverseObukhovLength`](@ref))
selects what the stability iteration solves for. Both are returned, as
`obukhov_length` and `inverse_obukhov_length`. Passing `stability` (as
`stability_variable`) skips the iteration and evaluates the profile at that
stability.
"""
function atmospheric_surface_profile!(boundary_layer_model::MoninObukhov, buffers;
    roughness_length,
    environment_instant,
    surface_temperature,
    vapour_pressure_equation=GoffGratch(),
    reference_height=last(buffers.heights),
    displacement_height=0.0u"m",
    thermal_roughness_model=boundary_layer_model.thermal_roughness_model,
    temperature_anchor_height=nothing,
    surface_relative_humidity=nothing,
    stability_variable=ObukhovLength(),
    stability=nothing,
)
    (; karman_constant) = boundary_layer_model
    (; atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity) = environment_instant

    (; heights, air_temperature, wind_speed, relative_humidity,
       obukhov_length_prev, warned_below_roughness, warned_above_reference) = buffers
    N_heights = length(heights)
    # `displacement_height`/`roughness_length` shift the log-profile origin to a
    # canopy-displaced height. Heights are measured from the true ground, so all
    # height-above-origin terms subtract `displacement_height`.
    reference_height = u"m"(reference_height)
    reference_height_above_origin = reference_height - displacement_height

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

    temperature_difference = reference_temperature - surface_temperature
    mean_temperature = (surface_temperature + reference_temperature) / 2
    # TODO call calc_air_volumetric_heat_capacity method specific to elevation and RH in final version but do it this way for NicheMapR comparison
    air_volumetric_heat_capacity = calc_air_volumetric_heat_capacity(mean_temperature)
    kinematic_viscosity = dry_air_properties(mean_temperature, atmospheric_pressure).kinematic_viscosity

    stability_result = surface_layer_stability(boundary_layer_model, stability_variable, stability, obukhov_length_prev, reference_temperature,
        surface_temperature, reference_wind_speed, roughness_length, reference_height_above_origin, temperature_difference, air_volumetric_heat_capacity, kinematic_viscosity;
        thermal_roughness_model)
    stability = stability_result.value
    roughness_length_temperature = stability_result.roughness_length_temperature
    convective_heat_flux = stability_result.convective_heat_flux
    friction_velocity = stability_result.friction_velocity
    reference_heat_correction = stability_result.heat_correction
    # See docstring: default anchors on roughness_length/roughness_length_temperature;
    # temperature_anchor_height switches to a direct two-point interpolation.
    use_anchor = !isnothing(temperature_anchor_height)
    anchor_height_above_origin = use_anchor ? u"m"(temperature_anchor_height) - displacement_height : u"m"(roughness_length)
    anchor_temperature = use_anchor ? surface_temperature : roughness_length_temperature
    anchor_heat_correction = heat_stability_correction(boundary_layer_model, stability_parameter_at_height(stability_variable, anchor_height_above_origin, stability))
    roughness_momentum_correction = momentum_stability_correction(boundary_layer_model, stability_parameter_at_height(stability_variable, roughness_length, stability))
    temperature_denominator = log(reference_height_above_origin / anchor_height_above_origin) + anchor_heat_correction - reference_heat_correction
    # Wind follows the reference wind, not the friction velocity, which has a floor in calm air
    reference_momentum_term = log(reference_height_above_origin / roughness_length) -
        momentum_stability_correction(boundary_layer_model, stability_parameter_at_height(stability_variable, reference_height_above_origin, stability)) + roughness_momentum_correction

    reference_vapor_pressure = wet_air_properties(reference_temperature, reference_humidity, atmospheric_pressure; vapour_pressure_equation).vapour_pressure
    # No surface_relative_humidity: constant vapour pressure with height, as before.
    # Given: reuses temperature's own log/reference_heat_correction ratio below (RHabove's approach).
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
            air_temperature[i] = surface_temperature + (roughness_length_temperature - surface_temperature) * height_fraction
            relative_humidity[i] = clamp(reference_vapor_pressure / vapour_pressure(vapour_pressure_equation, air_temperature[i]), 0.0, 1.0)
        else
            height_stability_parameter = stability_parameter_at_height(stability_variable, height_above_origin, stability)
            momentum_term = log(height_above_origin / roughness_length) - momentum_stability_correction(boundary_layer_model, height_stability_parameter) + roughness_momentum_correction
            wind_speed[i] = reference_wind_speed * momentum_term / reference_momentum_term
            # Temperature uses its own anchor (which may differ from roughness_length -- see above)
            temperature_log_term = log(height_above_origin / anchor_height_above_origin) + anchor_heat_correction - heat_stability_correction(boundary_layer_model, height_stability_parameter)
            air_temperature[i] = anchor_temperature + (reference_temperature - anchor_temperature) * temperature_log_term / temperature_denominator
            vapor_pressure_h = if use_humidity_anchor
                anchor_vapor_pressure + (reference_vapor_pressure - anchor_vapor_pressure) * temperature_log_term / temperature_denominator
            else
                reference_vapor_pressure
            end
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
    )
end

# No stability given: iterate it, warm-started from (and updating) the buffer's previous value
function surface_layer_stability(boundary_layer_model, stability_variable, ::Nothing, obukhov_length_prev, args...; kw...)
    out = monin_obukhov_iteration(boundary_layer_model, stability_variable, args...; max_iterations=30, tolerance=1e-2,
        initial_value=from_obukhov_length(stability_variable, obukhov_length_prev[]), kw...)
    obukhov_length_prev[] = to_obukhov_length(stability_variable, out.value)
    return out
end
# Stability given: one step at that value
surface_layer_stability(boundary_layer_model, stability_variable, stability, obukhov_length_prev, args...; kw...) =
    (; value=stability, monin_obukhov_step(boundary_layer_model, stability_variable, stability, args...; kw...)...)

function surface_fluxes(boundary_layer_model::MoninObukhov;
    surface_temperature, air_temperature, wind_speed, zenith_angle,
    roughness_length, reference_height, atmospheric_pressure, obukhov_length_prev=nothing,
)
    karman_constant = boundary_layer_model.karman_constant
    loglaw_ratio = reference_height / roughness_length
    loglaw_ratio > 1 || throw(DomainError(ustrip(loglaw_ratio),
        "surface_fluxes: reference_height ($reference_height) must exceed roughness_length " *
        "($roughness_length) for Monin-Obukhov log-law inversion"))
    loglaw_ratio < _MIN_LOGLAW_RATIO && @warn(
        "surface_fluxes: reference_height is close to roughness_length -- Monin-Obukhov " *
        "log-law inversion is unreliable in this roughness-sublayer regime",
        reference_height, roughness_length, loglaw_ratio, maxlog=10,
    )
    log_height_ratio = log(reference_height / roughness_length)
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
            calc_convection(; friction_velocity, log_height_ratio, temperature_difference, air_volumetric_heat_capacity, roughness_length, height=reference_height, karman_constant,
                kinematic_viscosity, thermal_roughness_model=boundary_layer_model.thermal_roughness_model))
    else
        # Unstable: iterate Obukhov length.
        previous_obukhov_length = isnothing(obukhov_length_prev) ? INITIAL_OBUKHOV_LENGTH : obukhov_length_prev[]
        initial_obukhov_length = previous_obukhov_length >= 0.0u"m" ? INITIAL_OBUKHOV_LENGTH : previous_obukhov_length
        out = calc_obukhov_length(boundary_layer_model, air_temperature, surface_temperature, wind_speed,
            roughness_length, reference_height, temperature_difference, air_volumetric_heat_capacity, kinematic_viscosity;
            initial_obukhov_length)
        friction_velocity = out.friction_velocity
        convective_heat_flux = out.convective_heat_flux
    end

    return (; convective_heat_flux, friction_velocity, temperature_difference)
end

function init_surface_fluxes!(boundary_layer_model::MoninObukhov, buffers, forcing, roughness_length, heights, soil_temperatures, i)
    t_next = ((i - 1) * 60)u"minute"
    (; air_temperature, wind_speed, zenith_angle, atmospheric_pressure) = interpolate_forcings(forcing, t_next)
    surface_temperature = soil_temperatures[1]
    if air_temperature < surface_temperature && zenith_angle < 90u"°"
        karman_constant = boundary_layer_model.karman_constant
        reference_height = last(heights)
        temperature_difference = air_temperature - surface_temperature
        mean_temperature = (surface_temperature + air_temperature) / 2
        air_volumetric_heat_capacity = calc_air_volumetric_heat_capacity(mean_temperature)
        kinematic_viscosity = dry_air_properties(mean_temperature, atmospheric_pressure).kinematic_viscosity
        previous_obukhov_length = buffers.soil_energy_balance.obukhov_length_prev[]
        initial_obukhov_length = previous_obukhov_length >= 0.0u"m" ? INITIAL_OBUKHOV_LENGTH : previous_obukhov_length
        out = calc_obukhov_length(boundary_layer_model, air_temperature, surface_temperature, wind_speed,
            roughness_length, reference_height, temperature_difference, air_volumetric_heat_capacity, kinematic_viscosity;
            initial_obukhov_length)
        buffers.soil_energy_balance.obukhov_length_prev[] = out.obukhov_length end
end

"""
    calc_air_volumetric_heat_capacity(mean_temperature)

Volumetric heat capacity of air (`J/m^3/K`) from temperature alone, ignoring moisture
and elevation.
"""
function calc_air_volumetric_heat_capacity(mean_temperature)
    return u"J/m^3/K"(u"(cal*g)/(g*cm^3*K)" * (0.08472 / ustrip(u"K", mean_temperature)))
end

"""
    calc_air_volumetric_heat_capacity(mean_temperature, relative_humidity, atmospheric_pressure;
                                      vapour_pressure_equation=GoffGratch())

Volumetric heat capacity of moist air: dry air density (`dry_air_properties`) times moist
air specific heat (`wet_air_properties`).
"""
function calc_air_volumetric_heat_capacity(mean_temperature, relative_humidity, atmospheric_pressure; vapour_pressure_equation=GoffGratch())
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
- `log_height_ratio::Real`: `log(height / roughness_length)`.
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
    calc_convection(; friction_velocity, log_height_ratio, temperature_difference, air_volumetric_heat_capacity, roughness_length, height, karman_constant, kinematic_viscosity, thermal_roughness_model=SublayerStantonRoughness())

Neutral convective heat flux between the surface and the air at `height`, with no
stability correction. `log_height_ratio` is `log(height / roughness_length)` and
`temperature_difference` is air minus surface. `thermal_roughness_model` sets how heat
roughness is found (see [`SublayerStantonRoughness`](@ref) and [`ScalarRoughnessRatio`](@ref)).

# See also
[`calc_friction_velocity`](@ref), [`sublayer_stanton`](@ref), [`bulk_stanton`](@ref), [`convective_flux`](@ref)
"""
function calc_convection(; friction_velocity, log_height_ratio, temperature_difference, air_volumetric_heat_capacity, roughness_length, height, karman_constant, kinematic_viscosity, thermal_roughness_model=SublayerStantonRoughness())
    return _stable_convective_heat_flux(thermal_roughness_model, air_volumetric_heat_capacity, temperature_difference, friction_velocity, karman_constant, roughness_length, height, log_height_ratio, kinematic_viscosity)
end

"""
    convective_flux(air_volumetric_heat_capacity, temperature_difference, friction_velocity, bulk_stanton_number, sublayer_stanton_number)

Compute convective heat flux given bulk and sublayer Stanton numbers.
"""
@inline function convective_flux(air_volumetric_heat_capacity, temperature_difference, friction_velocity, bulk_stanton_number, sublayer_stanton_number)
    return air_volumetric_heat_capacity * temperature_difference * friction_velocity * bulk_stanton_number /
        (1 + bulk_stanton_number / sublayer_stanton_number)
end

# Scalar-roughness-ratio resistance, ln(z/zh) - ψ_h(z/L) + ψ_h(zh/L)
function _scalar_roughness_resistance(ratio, karman_constant, friction_velocity, roughness_length, height, heat_correction_at_height, heat_correction_at_scalar_roughness)
    scalar_roughness_length = ratio * roughness_length
    return (log(height / scalar_roughness_length) - heat_correction_at_height + heat_correction_at_scalar_roughness) /
        (karman_constant * friction_velocity)
end

# calc_convection's stable/neutral direct branch: no Obukhov iteration, reference_heat_correction implicitly 0.
_stable_convective_heat_flux(::SublayerStantonRoughness, air_volumetric_heat_capacity, temperature_difference, friction_velocity, karman_constant, roughness_length, height, log_height_ratio, kinematic_viscosity) =
    convective_flux(air_volumetric_heat_capacity, temperature_difference, friction_velocity, bulk_stanton(log_height_ratio), sublayer_stanton(roughness_length, friction_velocity, kinematic_viscosity))
function _stable_convective_heat_flux(thermal_roughness_model::ScalarRoughnessRatio, air_volumetric_heat_capacity, temperature_difference, friction_velocity, karman_constant, roughness_length, height, log_height_ratio, kinematic_viscosity)
    resistance = _scalar_roughness_resistance(thermal_roughness_model.ratio, karman_constant, friction_velocity, roughness_length, height, 0.0, 0.0)
    return air_volumetric_heat_capacity * temperature_difference / resistance
end

# calc_obukhov_length's iterative branch. Also returns roughness_length_temperature
# (the below-z0 profile's own boundary value) since SublayerStantonRoughness derives
# it from the same bulk/sublayer split as the flux itself.
function _unstable_convective_heat_flux(::SublayerStantonRoughness, boundary_layer_model, air_volumetric_heat_capacity, temperature_difference, friction_velocity, karman_constant, roughness_length, height, log_ratio_corrected, heat_correction, stability_parameter, kinematic_viscosity, reference_temperature, surface_temperature)
    sublayer_stanton_number = sublayer_stanton(roughness_length, friction_velocity, kinematic_viscosity)
    bulk_stanton_number = _bulk_stanton(boundary_layer_model, log_ratio_corrected, heat_correction, stability_parameter, roughness_length, height)
    convective_heat_flux = convective_flux(air_volumetric_heat_capacity, temperature_difference, friction_velocity, bulk_stanton_number, sublayer_stanton_number)
    roughness_length_temperature = (reference_temperature * bulk_stanton_number + surface_temperature * sublayer_stanton_number) /
        (bulk_stanton_number + sublayer_stanton_number)
    return (; convective_heat_flux, roughness_length_temperature)
end
function _unstable_convective_heat_flux(thermal_roughness_model::ScalarRoughnessRatio, boundary_layer_model, air_volumetric_heat_capacity, temperature_difference, friction_velocity, karman_constant, roughness_length, height, log_ratio_corrected, heat_correction, stability_parameter, kinematic_viscosity, reference_temperature, surface_temperature)
    (; ratio) = thermal_roughness_model
    heat_correction_at_scalar_roughness = heat_stability_correction(boundary_layer_model, stability_parameter * ratio * roughness_length / height)
    resistance = _scalar_roughness_resistance(ratio, karman_constant, friction_velocity, roughness_length, height, heat_correction, heat_correction_at_scalar_roughness)
    return (; convective_heat_flux=air_volumetric_heat_capacity * temperature_difference / resistance, roughness_length_temperature=surface_temperature)
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
    sublayer_stanton(roughness_length, friction_velocity, kinematic_viscosity)

Stanton number for the roughness sublayer immediately above the surface,
`0.62·Re^-0.45`, with `Re = roughness_length · friction_velocity / kinematic_viscosity`.
Not suitable for a canopy.
"""
function sublayer_stanton(roughness_length, friction_velocity, kinematic_viscosity)
    roughness_reynolds_number = uconvert(NoUnits, roughness_length * friction_velocity / kinematic_viscosity)
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

# The bulk Stanton number at either stability. Unstable air uses the NicheMapR form above.
# In stable air (1 - 0.1·z/L) goes negative past z/L = 10, so stable air uses the heat log
# law with its own stability correction instead, 0.64/(ln(z/z0) - ψ_h(z/L) + ψ_h(z0/L)).
function _bulk_stanton(boundary_layer_model, log_ratio_corrected, heat_correction, stability_parameter, roughness_length, height)
    if stability_parameter > 0
        return 0.64 / (log(height / roughness_length) - heat_correction +
            heat_stability_correction(boundary_layer_model, stability_parameter * roughness_length / height))
    else
        return bulk_stanton(log_ratio_corrected, stability_parameter)
    end
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

# Convergence measure: the change in the stability parameter. A relative change in the
# Obukhov length never settles at neutral, where it heads to infinity.
stability_change(::ObukhovLength, height, new, old) = abs(uconvert(NoUnits, height * (inv(new) - inv(old))))
stability_change(::InverseObukhovLength, height, new, old) = abs(uconvert(NoUnits, height * (new - old)))

# Conversions to and from the Obukhov length kept in warm-start buffers
from_obukhov_length(::ObukhovLength, obukhov_length) = obukhov_length
from_obukhov_length(::InverseObukhovLength, obukhov_length) = inv(obukhov_length)
to_obukhov_length(::ObukhovLength, obukhov_length) = obukhov_length
to_obukhov_length(::InverseObukhovLength, inverse_obukhov_length) = inv(inverse_obukhov_length)
to_inverse_obukhov_length(::ObukhovLength, obukhov_length) = inv(obukhov_length)
to_inverse_obukhov_length(::InverseObukhovLength, inverse_obukhov_length) = inverse_obukhov_length
from_inverse_obukhov_length(::ObukhovLength, inverse_obukhov_length) = inv(inverse_obukhov_length)
from_inverse_obukhov_length(::InverseObukhovLength, inverse_obukhov_length) = inverse_obukhov_length


"""
    paulson_x(stability_parameter, dyer_constant)

`(1 - dyer_constant * stability_parameter)^(1/4)`, the input to the unstable stability corrections.

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
@inline function paulson_x(stability_parameter, dyer_constant)
    # sqrt is faster than ^1/4
    return sqrt(sqrt(1.0 - dyer_constant * stability_parameter))
end


"""
    heat_gradient_stability_factor(stability_parameter, dyer_constant)
    heat_gradient_stability_factor(model, stability_parameter)

Factor that divides the heat mixing coefficient: below 1 in unstable air, above 1 in stable air.
`model` is one side's corrections, an [`AbstractUnstableStability`](@ref) or
[`AbstractStableStability`](@ref), or a [`MoninObukhov`](@ref), which takes the side by the
sign of `stability_parameter`.

# References
- Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971).
  Flux–profile relationships in the atmospheric surface layer.
  *Journal of the Atmospheric Sciences*, 28(2), 181–189.
"""
function heat_gradient_stability_factor(stability_parameter, dyer_constant)
    return paulson_x(stability_parameter, dyer_constant)^(-2)
end


"""
    momentum_stability_correction(paulson_variable)
    momentum_stability_correction(model, stability_parameter)

Stability correction to the log wind profile. Zero in neutral air. `model` is as for
[`heat_gradient_stability_factor`](@ref).

# References
- Businger et al. (1971).
- Dyer (1974).
"""
@inline function momentum_stability_correction(paulson_variable)
    return 2.0 * log((1.0 + paulson_variable) / 2.0) + log((1.0 + paulson_variable^2) / 2.0) -
        2.0 * atan(paulson_variable) + π / 2.0
end


"""
    heat_stability_correction(paulson_variable)
    heat_stability_correction(model, stability_parameter)

Stability correction to the log temperature and vapour profiles. Zero in neutral air. `model`
is as for [`heat_gradient_stability_factor`](@ref).

# References
- Businger et al. (1971).
- Dyer (1974).
"""
@inline function heat_stability_correction(paulson_variable)
    return 2.0 * log((1.0 + paulson_variable^2.0) / 2.0)
end

# The boundary layer model's unstable or stable corrections, by the sign of the stability parameter
momentum_stability_correction(boundary_layer_model::MoninObukhov, stability_parameter) = stability_parameter < 0 ?
    momentum_stability_correction(boundary_layer_model.unstable_stability, stability_parameter) :
    momentum_stability_correction(boundary_layer_model.stable_stability, stability_parameter)
heat_stability_correction(boundary_layer_model::MoninObukhov, stability_parameter) = stability_parameter < 0 ?
    heat_stability_correction(boundary_layer_model.unstable_stability, stability_parameter) :
    heat_stability_correction(boundary_layer_model.stable_stability, stability_parameter)
heat_gradient_stability_factor(boundary_layer_model::MoninObukhov, stability_parameter) = stability_parameter < 0 ?
    heat_gradient_stability_factor(boundary_layer_model.unstable_stability, stability_parameter) :
    heat_gradient_stability_factor(boundary_layer_model.stable_stability, stability_parameter)

momentum_stability_correction(model::PaulsonDyer, stability_parameter) =
    momentum_stability_correction(paulson_x(stability_parameter, model.momentum_dyer_constant))
heat_stability_correction(model::PaulsonDyer, stability_parameter) =
    heat_stability_correction(paulson_x(stability_parameter, model.heat_dyer_constant))
heat_gradient_stability_factor(model::PaulsonDyer, stability_parameter) =
    heat_gradient_stability_factor(stability_parameter, model.heat_dyer_constant)

# NicheMapR bounds the stability parameter to ±10
momentum_stability_correction(model::NicheMapRUnstable, stability_parameter) =
    momentum_stability_correction(paulson_x(max(stability_parameter, -10.0), model.dyer_constant))
heat_stability_correction(model::NicheMapRUnstable, stability_parameter) =
    heat_stability_correction(paulson_x(max(stability_parameter, -10.0), model.dyer_constant))
heat_gradient_stability_factor(model::NicheMapRUnstable, stability_parameter) =
    heat_gradient_stability_factor(max(stability_parameter, -10.0), model.dyer_constant)

momentum_stability_correction(model::NicheMapRStable, stability_parameter) =
    clamp(-model.stable_beta * min(stability_parameter, 10.0), -4.0, 3.0)
heat_stability_correction(model::NicheMapRStable, stability_parameter) =
    clamp(-(model.stable_beta * min(stability_parameter, 10.0)) / model.turbulent_prandtl_number, -4.0, 3.0)
function heat_gradient_stability_factor(model::NicheMapRStable, stability_parameter)
    bounded = min(stability_parameter, 10.0)
    (; stable_heat_gradient_coefficient, min_stable_heat_gradient, max_stable_heat_gradient) = model
    return clamp(1.0 + stable_heat_gradient_coefficient * bounded / (1.0 + bounded), min_stable_heat_gradient, max_stable_heat_gradient)
end

momentum_stability_correction(model::BusingerLinear, stability_parameter) = -model.coefficient * stability_parameter
heat_stability_correction(model::BusingerLinear, stability_parameter) = -model.coefficient * stability_parameter
heat_gradient_stability_factor(model::BusingerLinear, stability_parameter) = 1 + model.coefficient * stability_parameter

function momentum_stability_correction(model::BeljaarsHoltslag, stability_parameter)
    (; a, b, c, d) = model
    return -(a * stability_parameter + b * (stability_parameter - c / d) * exp(-d * stability_parameter) + b * c / d)
end

function heat_stability_correction(model::BeljaarsHoltslag, stability_parameter)
    (; a, b, c, d) = model
    return -((1 + 2a * stability_parameter / 3)^1.5 + b * (stability_parameter - c / d) * exp(-d * stability_parameter) + b * c / d - 1)
end

function heat_gradient_stability_factor(model::BeljaarsHoltslag, stability_parameter)
    # 1 - stability_parameter times the derivative of heat_stability_correction
    (; a, b, c, d) = model
    return 1 + stability_parameter *
        (a * sqrt(1 + 2a * stability_parameter / 3) + b * (1 + c - d * stability_parameter) * exp(-d * stability_parameter))
end

function momentum_stability_correction(model::ChengBrutsaert, stability_parameter)
    (; momentum_scale, momentum_exponent) = model
    return -momentum_scale * log(stability_parameter + (1 + stability_parameter^momentum_exponent)^(1 / momentum_exponent))
end

function heat_stability_correction(model::ChengBrutsaert, stability_parameter)
    (; heat_scale, heat_exponent) = model
    return -heat_scale * log(stability_parameter + (1 + stability_parameter^heat_exponent)^(1 / heat_exponent))
end

function heat_gradient_stability_factor(model::ChengBrutsaert, stability_parameter)
    # 1 - stability_parameter times the derivative of heat_stability_correction
    (; heat_scale, heat_exponent) = model
    power_sum = 1 + stability_parameter^heat_exponent
    log_argument = stability_parameter + power_sum^(1 / heat_exponent)
    # stability_parameter^heat_exponent rather than stability_parameter times its (heat_exponent - 1)
    # power, whose derivative is infinite at neutral, so automatic differentiation stays finite there
    return 1 + heat_scale * (stability_parameter + stability_parameter^heat_exponent * power_sum^(1 / heat_exponent - 1)) /
        log_argument
end


"""
    monin_obukhov_step(boundary_layer_model, stability_variable, value, reference_temperature, surface_temperature, reference_wind_speed,
                       roughness_length, reference_height_above_origin, temperature_difference, air_volumetric_heat_capacity, kinematic_viscosity;
                       thermal_roughness_model=boundary_layer_model.thermal_roughness_model)

One evaluation of the Monin–Obukhov relations at stability `value` (as
`stability_variable`, [`ObukhovLength`](@ref) or [`InverseObukhovLength`](@ref)):
friction velocity and convective heat flux from the stability corrections,
then `stability_update`, the stability those fluxes imply through the
Obukhov relation. `value` is self-consistent when `stability_update == value`.
Works for either sign of `temperature_difference`.
"""
function monin_obukhov_step(boundary_layer_model, stability_variable, value, reference_temperature, surface_temperature, reference_wind_speed,
    roughness_length, reference_height_above_origin, temperature_difference, air_volumetric_heat_capacity, kinematic_viscosity;
    thermal_roughness_model=boundary_layer_model.thermal_roughness_model,
)
    (; karman_constant, min_friction_velocity) = boundary_layer_model
    reference_stability_parameter = stability_parameter_at_height(stability_variable, reference_height_above_origin, value)
    reference_momentum_correction = momentum_stability_correction(boundary_layer_model, reference_stability_parameter)
    reference_heat_correction = heat_stability_correction(boundary_layer_model, reference_stability_parameter)
    log_ratio_corrected = log(reference_height_above_origin / roughness_length) - reference_momentum_correction +
        momentum_stability_correction(boundary_layer_model, reference_stability_parameter * roughness_length / reference_height_above_origin)
    friction_velocity = max(karman_constant * reference_wind_speed / log_ratio_corrected, min_friction_velocity)
    convective_heat_flux, roughness_length_temperature = _unstable_convective_heat_flux(thermal_roughness_model,
        boundary_layer_model, air_volumetric_heat_capacity, temperature_difference, friction_velocity, karman_constant, roughness_length, reference_height_above_origin, log_ratio_corrected, reference_heat_correction, reference_stability_parameter, kinematic_viscosity,
        reference_temperature, surface_temperature)
    stability_update = stability_from_fluxes(stability_variable, OBUKHOV_LENGTH_FACTOR, friction_velocity, convective_heat_flux)
    return (;
        stability_update,
        friction_velocity,
        heat_correction=reference_heat_correction,
        convective_heat_flux=u"W/m^2"(convective_heat_flux),
        roughness_length_temperature,
    )
end

"""
    monin_obukhov_iteration(boundary_layer_model, stability_variable, args...; initial_value, max_iterations=50, tolerance=1e-2, relaxation=0.5, kw...)

Iterate [`monin_obukhov_step`](@ref) (same `args` and `kw`) from
`initial_value` until the stability stops changing, under-relaxing each
update by `relaxation` in inverse Obukhov length. Returns the last step's fluxes and `value`, the
stability that step was evaluated at.
"""
function monin_obukhov_iteration(boundary_layer_model, stability_variable, reference_temperature, surface_temperature, reference_wind_speed,
    roughness_length, reference_height_above_origin, temperature_difference, air_volumetric_heat_capacity, kinematic_viscosity;
    initial_value, max_iterations=50, tolerance=1e-2, relaxation=0.5, kw...,
)
    value = initial_value
    iteration = 1
    while true
        step = monin_obukhov_step(boundary_layer_model, stability_variable, value, reference_temperature, surface_temperature, reference_wind_speed,
            roughness_length, reference_height_above_origin, temperature_difference, air_volumetric_heat_capacity, kinematic_viscosity; kw...)
        change = stability_change(stability_variable, reference_height_above_origin, step.stability_update, value)
        (change <= tolerance || iteration >= max_iterations) && return (; value, step...)
        # Relax in 1/L, which passes through zero at neutral. Averaging two Obukhov lengths
        # of opposite sign can land near L = 0, an extreme stability.
        inverse_value = relaxation * to_inverse_obukhov_length(stability_variable, step.stability_update) +
            (1.0 - relaxation) * to_inverse_obukhov_length(stability_variable, value)
        value = from_inverse_obukhov_length(stability_variable, inverse_value)
        iteration += 1
    end
end

"""
    calc_obukhov_length(boundary_layer_model, reference_temperature, surface_temperature, reference_wind_speed, roughness_length,
                        reference_height_above_origin, temperature_difference, air_volumetric_heat_capacity, kinematic_viscosity;
                        initial_obukhov_length, kw...)

[`monin_obukhov_iteration`](@ref) on the Obukhov length. Converges to `±Inf`
(neutral) as `temperature_difference → 0`.
"""
function calc_obukhov_length(boundary_layer_model, args...; initial_obukhov_length, kw...)
    (; value, friction_velocity, heat_correction, convective_heat_flux, roughness_length_temperature) =
        monin_obukhov_iteration(boundary_layer_model, ObukhovLength(), args...; initial_value=initial_obukhov_length, kw...)
    return (; obukhov_length=value, friction_velocity, heat_correction, convective_heat_flux, roughness_length_temperature)
end

"""
    inverse_obukhov_length_limits(stability_limit, boundary_layer_model, height, roughness_length)

The range of inverse Obukhov length `stability_limit` allows, as `(; unstable, stable)`.

[`NoStabilityLimit`](@ref) is unbounded. [`MomentumCorrectionLimit`](@ref) gives the values
nearest to neutral where the momentum stability correction between `roughness_length` and
`height` (above the displacement height), `ψ_m(z₀/L) - ψ_m(z/L)`, reaches `fraction` of the log
term `ln(z/z₀)`, or an infinite value where it never does. Within these limits the friction
velocity changes with stability by at most a factor `1/(1 + fraction)` in stable air and
`1/(1 - fraction)` in unstable air.
"""
# |1/L| from 2⁻³⁰ to 2⁴⁰ m⁻¹ spans every stability a surface layer can have
const LIMIT_SEARCH_EXPONENTS = -30:40
# The limits' precision, relative to their size
const LIMIT_RELATIVE_TOLERANCE = 1e-12

inverse_obukhov_length_limits(::NoStabilityLimit, boundary_layer_model, height, roughness_length) =
    (; unstable=-Inf * u"m^-1", stable=Inf * u"m^-1")
function inverse_obukhov_length_limits(limit::MomentumCorrectionLimit, boundary_layer_model, height, roughness_length)
    (; fraction) = limit
    target = fraction * log(height / roughness_length)
    # 1/L in m⁻¹
    correction(inverse_length) = momentum_stability_correction(boundary_layer_model, roughness_length * inverse_length * u"m^-1") -
        momentum_stability_correction(boundary_layer_model, height * inverse_length * u"m^-1")
    nearest_to_neutral(sign) = begin
        excess(inverse_length, _) = sign * correction(inverse_length) - target
        # The first doubling of |1/L| past the target brackets it
        previous = 0.0
        for exponent in LIMIT_SEARCH_EXPONENTS
            inverse_length = sign * 2.0^exponent
            if excess(inverse_length, nothing) >= 0
                problem = IntervalNonlinearProblem{false}(excess, minmax(previous, inverse_length))
                return solve(problem, ITP(); abstol=LIMIT_RELATIVE_TOLERANCE * abs(inverse_length)).u * u"m^-1"
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
