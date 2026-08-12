"""
    MoninObukhov(; karman_constant=0.4, dyer_constant=16.0, stable_beta=4.7,
                   turbulent_prandtl_number=0.74, stable_Φ_h_coefficient=6.0,
                   min_stable_Φ_h=0.5, max_stable_Φ_h=1.5, min_friction_velocity=0.02u"m/s",
                   thermal_roughness_model=SublayerStantonRoughness())

Monin–Obukhov similarity theory boundary-layer formulation. Holds the
empirical constants of the unstable (Businger-Dyer/Paulson, `dyer_constant`)
and stable Φ relations. Previously these constants were buried as fields on
`MicroTerrain`, where they masqueraded as terrain properties.

Two distinct stable-branch parameterizations, matching R's `microclimlearn`
having two distinct functions for them (`dpsih` vs `dphih`) rather than one:
- `stable_beta`/`turbulent_prandtl_number`: linear form (`ψ_m = -stable_beta·z/L`,
  `ψ_h = -(stable_beta/turbulent_prandtl_number)·z/L`, `ζ=z/L≥0`) for the
  log-profile offsets [`calc_ψ_m`](@ref)/[`calc_ψ_h`](@ref) use, matching R's
  `dpsim`/`dpsih` (Businger et al. 1971).
- `stable_Φ_h_coefficient`/`min_stable_Φ_h`/`max_stable_Φ_h`: saturating form
  (`Φ_h = clamp(1 + stable_Φ_h_coefficient·ζ/(1+ζ), min_stable_Φ_h,
  max_stable_Φ_h)`) for the bulk diffusivity multiplier [`calc_Φ_h`](@ref)
  uses, matching R's `dphih`. Saturates by construction
  (→ `1+stable_Φ_h_coefficient` as `ζ→∞`) even before the explicit clamp,
  keeping `RaupachLTheoryAirProfile`'s `T_L` bounded under strongly stable,
  low-wind conditions.

# References
Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971).
Flux–profile relationships in the atmospheric surface layer.
*Journal of the Atmospheric Sciences*, 28(2), 181–189.

Dyer, A. J. (1974). A review of flux–profile relationships.
*Boundary-Layer Meteorology*, 7(3), 363–372.
"""
@kwdef struct MoninObukhov{KC,DC,SB,PR,PHC,MINPH,MAXPH,MINUF,TRM} <: AbstractBoundaryLayerModel
    karman_constant::KC = 0.4
    dyer_constant::DC = 16.0
    stable_beta::SB = 4.7
    turbulent_prandtl_number::PR = 0.74
    stable_Φ_h_coefficient::PHC = 6.0
    min_stable_Φ_h::MINPH = 0.5
    max_stable_Φ_h::MAXPH = 1.5
    # Floors 1/friction_velocity terms (T_L, aerodynamic resistances) at a
    # "never fully calm" value.
    min_friction_velocity::MINUF = 0.02u"m/s"
    # z0-to-z0h correction for the sensible heat transfer coefficient. Default
    # matches prior behaviour; canopy-top z0 scales should use ScalarRoughnessRatio
    # instead (see its docstring).
    thermal_roughness_model::TRM = SublayerStantonRoughness()
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
Default, matching this package's prior behaviour. Confirmed (see
`ScalarRoughnessRatio`) to produce a far larger excess resistance than
microclimlearn's own approach once `z0` reaches canopy-top scale (~1 m);
appropriate at the bare-ground/short-roughness scale it's more plausibly
suited to.
"""
struct SublayerStantonRoughness <: AbstractThermalRoughnessModel end

"""
    ScalarRoughnessRatio(; ratio=0.2)

Fixed scalar-to-momentum roughness ratio `z0h = ratio·z0` (`kB⁻¹ =
-log(ratio)`, constant, not `Re*`-dependent) -- matches microclimlearn's
`canopyresistance` (`zh = 0.2·zm`) exactly. Confirmed against R's own
`solve_wholecanopy` output for a real canopy case (`z0≈1.9 m`): recovers
`L≈-64 m` against R's `-65 m`, versus `SublayerStantonRoughness`'s `≈-1500 m`
for the same inputs.
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
`convective_heat_flux` (`W/m²`), `friction_velocity` (`m/s`).

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
atmospheric_surface_profile(bl::MoninObukhov; heights=DEFAULT_HEIGHTS, reference_height=last(heights), kw...) =
    atmospheric_surface_profile!(bl, allocate_profile(bl, heights); reference_height, kw...)
"""
    atmospheric_surface_profile!(bl::MoninObukhov, buffers; site, environment_instant,
                                  surface_temperature, vapour_pressure_equation=GoffGratch(),
                                  reference_height=last(buffers.heights),
                                  displacement_height=0.0u"m", roughness_length=site.roughness_height)

`displacement_height`/`roughness_length` shift the effective log-profile
origin to a canopy-displaced height with a canopy-derived roughness length,
instead of bare ground. Defaults (`0.0u"m"` / `site.roughness_height`)
give the bare-ground profile. Heights are measured from the true ground (as
in `MicroModel.heights`); `displacement_height` is subtracted internally
wherever the log-law needs height above the displaced origin.

`temperature_anchor_height`: by default (`nothing`) the temperature profile
is anchored at `z0` using the aerodynamic (thermal-roughness-corrected)
`roughness_height_temperature`, as before. Pass an absolute height (e.g. a
canopy's `canopy_height`) to instead anchor directly on `surface_temperature`
*at that height* -- for a case like canopy top, `surface_temperature` is
already a genuine known air temperature at a real height, not a skin
temperature needing a `z0`/`z0h` sublayer correction, so this bypasses that
correction and does a plain two-point MOST log-law interpolation between
`(temperature_anchor_height, surface_temperature)` and
`(reference_height, reference_temperature)` instead -- matching
microclimlearn's `Tabove` (`micropoint_new2.cpp`), which anchors the
above-canopy profile the same way, including the `ψ_h` evaluated at the
anchor height itself (not assumed ≈0, since unlike `z0` this height need
not be roughness-sublayer-scale).

`surface_relative_humidity` (only used when `temperature_anchor_height` is
given): the known relative humidity *at* `temperature_anchor_height`,
anchoring vapour pressure the same two-point way (reusing the same log/ψ_h
ratio computed for temperature) instead of the default assumption of
constant vapour pressure with height -- matches microclimlearn's `RHabove`,
which is `Tabove`'s vapour-pressure counterpart.
"""
function atmospheric_surface_profile!(bl::MoninObukhov, buffers;
    site,
    environment_instant,
    surface_temperature,
    vapour_pressure_equation=GoffGratch(),
    reference_height=last(buffers.heights),
    displacement_height=0.0u"m",
    roughness_length=site.roughness_height,
    thermal_roughness_model=bl.thermal_roughness_model,
    temperature_anchor_height=nothing,
    surface_relative_humidity=nothing,
)
    (; elevation) = site
    (; karman_constant, dyer_constant, stable_beta, turbulent_prandtl_number, min_friction_velocity) = bl
    (; atmospheric_pressure, reference_temperature, reference_wind_speed, reference_humidity) = environment_instant

    (; heights, air_temperature, wind_speed, relative_humidity,
       obukhov_length_prev, warned_below_roughness, warned_above_reference) = buffers
    N_heights = length(heights)
    # `displacement_height`/`roughness_length` shift the log-profile origin to a
    # canopy-displaced height. Heights are measured from the true ground, so all
    # height-above-origin terms use `h - displacement_height`.
    ref_height = u"m"(reference_height)
    z0 = roughness_length
    z = ref_height - displacement_height
    # If the reference height appears in the heights array, pin its slot to measured values
    # to avoid the small inconsistency in ψ_h between the penultimate and final Obukhov iterations.
    ref_idx = findfirst(==(ref_height), heights)

    if !warned_below_roughness[] && minimum(heights) - displacement_height < z0
        @warn """Some requested heights are below the roughness length ($z0) above the displacement height ($displacement_height).
    Monin-Obukhov similarity theory is not valid in this sublayer region -- expected
    for below-canopy heights canopy model is running.
    Assumptions applied:
      • wind speed → 0 (log-law gives u = 0 at z = z0 by definition)
      • air temperature linearly interpolated between T_surface (z = 0) and T_z0 (z = z0)
    For a non-zero wind floor (e.g. for convection terms), clamp the result downstream."""
        warned_below_roughness[] = true
    end
    if !warned_above_reference[] && maximum(heights) > ref_height
        @warn """Some requested heights exceed the reference measurement height ($ref_height).
    Profiles will be extrapolated above the measurement level using MOST.
    This is physically reasonable within the surface layer but note that
    the log-law becomes less accurate as height approaches the boundary-layer depth."""
        warned_above_reference[] = true
    end

    reference_temp = u"K"(reference_temperature)
    surface_temp = u"K"(surface_temperature)
    κ = karman_constant
    γ = dyer_constant
    v_ref_height = reference_wind_speed

    # compute rcptkg (was a constant in original Fortran version)
    # dry_air_out = dry_air_properties(u"K"(reference_temperature), P_atmos)
    # wet_air_out = wet_air_properties(u"K"(reference_temperature), rh=reference_humidity, P_atmos)
    # ρ = dry_air_out.ρ_air
    # c_p = wet_air_out.specific_heat
    # TODO make this work with SI units
    #ρcpTκg = u"cal*minute^2/cm^4"(ρ * c_p * T_ref_height / (κ * g_n))
    ρcpTκg = u"J*s^2/m^4"(6.003e-8u"cal*minute^2/cm^4")

    ΔT = reference_temp - surface_temp
    mean_temp = (surface_temp + reference_temp) / 2
    # TODO call calc_ρ_cp method specific to elevation and RH in final version but do it this way for NicheMapR comparison
    ρ_cp = calc_ρ_cp(mean_temp)#, elevation, reference_humidity)
    kinematic_viscosity = dry_air_properties(mean_temp, atmospheric_pressure).kinematic_viscosity

    Obukhov_out = calc_Obukhov_length(reference_temp, surface_temp, v_ref_height, z0, z, ρcpTκg, κ, ΔT, ρ_cp, kinematic_viscosity;
        γ, stable_beta, turbulent_prandtl_number, max_iter=30, tol=1e-2, initial_obukhov_length=obukhov_length_prev[],
        min_friction_velocity, thermal_roughness_model)
    obukhov_length = Obukhov_out.obukhov_length
    obukhov_length_prev[] = obukhov_length
    roughness_height_temp = Obukhov_out.roughness_height_temperature
    convective_heat_flux = Obukhov_out.convective_heat_flux
    friction_velocity = Obukhov_out.friction_velocity
    ψ_h = Obukhov_out.ψ_h
    # See docstring: default anchors on z0/roughness_height_temp; temperature_anchor_height
    # switches to a direct two-point interpolation, ψ_h included at the anchor.
    use_anchor = !isnothing(temperature_anchor_height)
    z_temp_anchor = use_anchor ? u"m"(temperature_anchor_height) - displacement_height : u"m"(z0)
    anchor_temp = use_anchor ? surface_temp : roughness_height_temp
    ψ_h_anchor = use_anchor ? calc_ψ_h(z_temp_anchor, γ, obukhov_length, stable_beta, turbulent_prandtl_number) : 0.0
    # Floor: under strong instability ψ_h can approach/exceed log(z/z_temp_anchor),
    # blowing up every non-reference air_temperature.
    ref_temp_log_ratio = log(z / z_temp_anchor)
    temp_denom = max(ref_temp_log_ratio + ψ_h_anchor - ψ_h, 0.1 * ref_temp_log_ratio, 1e-6)

    reference_vapor_pressure = wet_air_properties(reference_temp, reference_humidity, atmospheric_pressure; vapour_pressure_equation).vapour_pressure
    # No surface_relative_humidity: constant vapour pressure with height, as before.
    # Given: reuses temperature's own log/ψ_h ratio below (RHabove's approach).
    use_humidity_anchor = use_anchor && !isnothing(surface_relative_humidity)
    anchor_vapor_pressure = use_humidity_anchor ?
        vapour_pressure(vapour_pressure_equation, surface_temp) * surface_relative_humidity : reference_vapor_pressure

    for i in 1:N_heights
        h = heights[i] - displacement_height
        if h < z0
            # h goes negative below the displaced origin -- clamp so the
            # interpolation saturates at surface_temp instead of extrapolating.
            frac = clamp(ustrip(h / z0), 0.0, 1.0)
            wind_speed[i] = zero(friction_velocity)
            air_temperature[i] = surface_temp + (roughness_height_temp - surface_temp) * frac
            relative_humidity[i] = clamp(reference_vapor_pressure / vapour_pressure(vapour_pressure_equation, air_temperature[i]), 0.0, 1.0)
        elseif i == ref_idx
            wind_speed[i] = v_ref_height
            air_temperature[i] = reference_temp
            relative_humidity[i] = reference_humidity
        else
            ψ_m1 = calc_ψ_m(h, γ, obukhov_length, stable_beta)
            ψ_h2 = calc_ψ_h(h, γ, obukhov_length, stable_beta, turbulent_prandtl_number)
            h_ratio = h / z0  # dimensionless h/z0
            # Floor at a fraction of the neutral log term (see calc_Obukhov_length) --
            # an absolute constant is meaningless when log(h_ratio) is itself O(1).
            log_h_ratio = log(h_ratio)
            wind_log_arg = max(log_h_ratio - ψ_m1, 0.1 * log_h_ratio, 1e-6)
            wind_speed[i] = (friction_velocity / κ) * wind_log_arg
            # Temperature uses its own anchor/ratio (z_temp_anchor may differ
            # from z0 -- see above); only floored when h is above the anchor,
            # the guaranteed regime for actual anchor-mode callers (below-anchor
            # heights, e.g. below canopy top in solve_air!'s full-grid profile,
            # are unused downstream and left unfloored but finite).
            log_temp_h_ratio = log(h / z_temp_anchor)
            temp_num = log_temp_h_ratio + ψ_h_anchor - ψ_h2
            temp_log_arg = log_temp_h_ratio >= 0 ?
                max(temp_num, 0.1 * log_temp_h_ratio, 1e-6) : temp_num
            air_temperature[i] = anchor_temp + (reference_temp - anchor_temp) * temp_log_arg / temp_denom
            vapor_pressure_h = use_humidity_anchor ?
                anchor_vapor_pressure + (reference_vapor_pressure - anchor_vapor_pressure) * temp_log_arg / temp_denom :
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
        obukhov_length,
    )
end

function surface_fluxes(bl::MoninObukhov;
    surface_temperature, air_temperature, wind_speed, zenith_angle,
    roughness_height, reference_height, atmospheric_pressure, obukhov_length_prev=nothing,
)
    κ = bl.karman_constant
    loglaw_ratio = reference_height / roughness_height
    loglaw_ratio > 1 || throw(DomainError(ustrip(loglaw_ratio),
        "surface_fluxes: reference_height ($reference_height) must exceed roughness_height " *
        "($roughness_height) for Monin-Obukhov log-law inversion"))
    loglaw_ratio < _MIN_LOGLAW_RATIO && @warn(
        "surface_fluxes: reference_height is close to roughness_height -- Monin-Obukhov " *
        "log-law inversion is unreliable in this roughness-sublayer regime",
        reference_height, roughness_height, loglaw_ratio, maxlog=10,
    )
    log_z_ratio = log(reference_height / roughness_height)
    ΔT = air_temperature - surface_temperature
    mean_temp = (surface_temperature + air_temperature) / 2
    ρ_cp = calc_ρ_cp(mean_temp)
    kinematic_viscosity = dry_air_properties(mean_temp, atmospheric_pressure).kinematic_viscosity

    if air_temperature ≥ surface_temperature || zenith_angle ≥ 90°
        # Stable / nocturnal: neutral log-law. Real reanalysis wind speed
        # occasionally rounds to exactly 0 (e.g. BARRA reports in 0.25 m/s
        # steps and does report dead calm). Without a floor, friction_velocity=0
        # sends sublayer_stanton (∝ u*^-0.45) to Inf, blowing up
        # convective_heat_flux and destabilising the ODE. Mirrors the
        # existing guard in calc_Obukhov_length's unstable branch below.
        friction_velocity = max(
            calc_friction_velocity(; reference_wind_speed=wind_speed, log_z_ratio, κ), bl.min_friction_velocity
        )
        convective_heat_flux = uconvert(u"W/m^2",
            calc_convection(; friction_velocity, log_z_ratio, ΔT, ρ_cp, z0=roughness_height, z=reference_height, κ,
                kinematic_viscosity, thermal_roughness_model=bl.thermal_roughness_model))
    else
        # Unstable: iterate Obukhov length. A warm-started value from a
        # previous hour may have converged to positive or near-zero, which
        # would cause NaN in calc_φ_m via sqrt of a negative number — guard
        # against that by resetting to the default unstable guess.
        ρcpTκg = 6.003e-8u"cal*minute^2/cm^4"
        L0_raw = isnothing(obukhov_length_prev) ? -0.3u"m" : obukhov_length_prev[]
        L0 = L0_raw >= 0.0u"m" ? -0.3u"m" : L0_raw
        out = calc_Obukhov_length(air_temperature, surface_temperature, wind_speed,
            roughness_height, reference_height, ρcpTκg, κ, ΔT, ρ_cp, kinematic_viscosity;
            initial_obukhov_length=L0, min_friction_velocity=bl.min_friction_velocity,
            thermal_roughness_model=bl.thermal_roughness_model)
        friction_velocity = out.friction_velocity
        convective_heat_flux = out.convective_heat_flux
    end

    return (; convective_heat_flux, friction_velocity, ΔT, ρ_cp)
end

function init_surface_fluxes!(bl::MoninObukhov, buffers, forcing, site, heights, T0, i)
    t_next = ((i - 1) * 60)u"minute"
    (; air_temperature, wind_speed, zenith_angle) = interpolate_forcings(forcing, t_next)
    surface_temperature = T0[1]
    if air_temperature < surface_temperature && zenith_angle < 90u"°"
        κ = bl.karman_constant
        roughness_height = site.roughness_height
        reference_height = last(heights)
        ΔT = air_temperature - surface_temperature
        mean_temp = (surface_temperature + air_temperature) / 2
        ρ_cp = calc_ρ_cp(mean_temp)
        kinematic_viscosity = dry_air_properties(mean_temp, site.atmospheric_pressure).kinematic_viscosity
        ρcpTκg = 6.003e-8u"cal*minute^2/cm^4"
        L0_raw = buffers.soil_energy_balance.obukhov_length_prev[]
        L0 = L0_raw >= 0.0u"m" ? -0.3u"m" : L0_raw
        out = calc_Obukhov_length(air_temperature, surface_temperature, wind_speed,
            roughness_height, reference_height, ρcpTκg, κ, ΔT, ρ_cp, kinematic_viscosity;
            initial_obukhov_length=L0, min_friction_velocity=bl.min_friction_velocity,
            thermal_roughness_model=bl.thermal_roughness_model)
        buffers.soil_energy_balance.obukhov_length_prev[] = out.obukhov_length end
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
    calc_friction_velocity(; reference_wind_speed, log_z_ratio, κ=0.4)

Compute the friction velocity (u*) from a reference wind speed using the
logarithmic wind profile.

# Arguments
- `reference_wind_speed::Quantity{<:Real,𝐋/𝐓}`: Wind speed at the reference height (e.g. `m/s`, `cm/min`).
- `log_z_ratio::Real`: Precomputed log height ratio, typically `log(z/z0 + 1.0)`.
- `κ::Real`: von Kármán constant (default = 0.4).

# Returns
- Friction velocity `friction_velocity::Quantity{<:Real,𝐋/𝐓}`.

# See also
[`calc_convection`](@ref), [`calc_wind`](@ref)
"""
function calc_friction_velocity(; reference_wind_speed, log_z_ratio, κ=0.4)
    v_ref_height = reference_wind_speed
    return κ * v_ref_height / log_z_ratio
end

"""
    calc_wind(z, z0, κ, friction_velocity, b)

Calculate wind speed at height `z` using the logarithmic wind profile.

# Arguments
- `z::Quantity{<:Real,𝐋}`: Height above the surface (e.g. `m`, `cm`).
- `z0::Quantity{<:Real,𝐋}`: Roughness length (e.g. `m`, `cm`).
- `κ::Real`: von Kármán constant.
- `friction_velocity::Quantity{<:Real,𝐋/𝐓}`: Friction velocity.
- `b::Real`: Offset term (e.g. `1.0` for neutral stability, or stability correction).

# Returns
- Wind speed at height `z::Quantity{<:Real,𝐋/𝐓}`.

# See also
[`calc_friction_velocity`](@ref), [`calc_convection`](@ref)
"""
function calc_wind(z, z0, κ, friction_velocity, b)
    return (friction_velocity / κ) * log(z / z0 + b)
end


"""
    calc_convection(; friction_velocity, log_z_ratio, ΔT, ρ_cp, z0)

Calculate the convective heat flux (sensible heat exchange between surface and air).

# Arguments
- `friction_velocity::Quantity{<:Real,𝐋/𝐓}`: Friction velocity (e.g. `m/s`, `cm/min`).
- `log_z_ratio::Real`: Precomputed logarithmic height ratio, typically `log(z/z0 + 1.0)`.
- `ΔT::Quantity{<:Real,Θ}`: Temperature difference between reference air and surface (Kelvin).
- `ρ_cp::Quantity{<:Real,(𝐌*𝐋^-1*𝐓^-2)}`: Volumetric heat capacity of air (e.g. `J/m³/K`, `cal/cm³/K`).
- `z0::Quantity{<:Real,𝐋}`: Surface roughness length (length).

# Returns
- Convective heat flux as `Quantity{<:Real,(𝐌*𝐓^-3)}` (e.g. `W/m²`, `cal/min/cm²`).

Uses bulk and sublayer Stanton numbers to account for turbulence near the surface.

# See also
[`calc_friction_velocity`](@ref), [`calc_wind`](@ref), [`sublayer_stanton`](@ref), [`bulk_stanton`](@ref), [`convective_flux`](@ref)
"""
function calc_convection(; friction_velocity, log_z_ratio, ΔT, ρ_cp, z0, z, κ, kinematic_viscosity, thermal_roughness_model=SublayerStantonRoughness())
    return _stable_convective_heat_flux(thermal_roughness_model, ρ_cp, ΔT, friction_velocity, κ, z0, z, log_z_ratio, kinematic_viscosity)
end

"""
    convective_flux(ρ_cp, ΔT, friction_velocity, St_bulk, St_sublayer)

Compute convective heat flux given bulk and sublayer Stanton numbers.
"""
@inline function convective_flux(ρ_cp, ΔT, friction_velocity, bulk_stanton_number, sublayer_stanton_number)
    return ρ_cp * ΔT * friction_velocity * bulk_stanton_number / (1 + bulk_stanton_number / sublayer_stanton_number)
end

# Scalar-roughness-ratio resistance, shared by the stable/neutral direct
# branch (ψ_h=0) and the unstable iterative branch below.
@inline function _scalar_roughness_resistance(ratio, κ, friction_velocity, z0, z, ψ_h)
    zh = ratio * z0
    log_z_ratio_h = log(z / zh)
    log_ratio_h_corrected = max(log_z_ratio_h - ψ_h, 0.1 * log_z_ratio_h, 1.0e-6)
    return log_ratio_h_corrected / (κ * friction_velocity)
end

# calc_convection's stable/neutral direct branch: no Obukhov iteration, ψ_h implicitly 0.
_stable_convective_heat_flux(::SublayerStantonRoughness, ρ_cp, ΔT, friction_velocity, κ, z0, z, log_z_ratio, kinematic_viscosity) =
    convective_flux(ρ_cp, ΔT, friction_velocity, bulk_stanton(log_z_ratio), sublayer_stanton(z0, friction_velocity, kinematic_viscosity))
function _stable_convective_heat_flux(m::ScalarRoughnessRatio, ρ_cp, ΔT, friction_velocity, κ, z0, z, log_z_ratio, kinematic_viscosity)
    resistance = _scalar_roughness_resistance(m.ratio, κ, friction_velocity, z0, z, 0.0)
    return ρ_cp * ΔT / resistance
end

# calc_Obukhov_length's iterative branch. Also returns roughness_height_temperature
# (the below-z0 profile's own boundary value) since SublayerStantonRoughness derives
# it from the same bulk/sublayer split as the flux itself.
function _unstable_convective_heat_flux(::SublayerStantonRoughness, ρ_cp, ΔT, friction_velocity, κ, z0, z, log_ratio_corrected, ψ_h, obukhov_length, kinematic_viscosity, reference_temp, surface_temp)
    sublayer_stanton_number = sublayer_stanton(z0, friction_velocity, kinematic_viscosity)
    bulk_stanton_number = bulk_stanton(log_ratio_corrected, z, _floor_obukhov_length(z, obukhov_length))
    convective_heat_flux = convective_flux(ρ_cp, ΔT, friction_velocity, bulk_stanton_number, sublayer_stanton_number)
    roughness_height_temperature = (reference_temp * bulk_stanton_number + surface_temp * sublayer_stanton_number) / (bulk_stanton_number + sublayer_stanton_number)
    return (; convective_heat_flux, roughness_height_temperature)
end
function _unstable_convective_heat_flux(m::ScalarRoughnessRatio, ρ_cp, ΔT, friction_velocity, κ, z0, z, log_ratio_corrected, ψ_h, obukhov_length, kinematic_viscosity, reference_temp, surface_temp)
    resistance = _scalar_roughness_resistance(m.ratio, κ, friction_velocity, z0, z, ψ_h)
    return (; convective_heat_flux=ρ_cp * ΔT / resistance, roughness_height_temperature=surface_temp)
end


"""
    calc_heat_transfer_coefficient(convective_heat_flux, ΔT)

Heat transfer coefficient from convective heat flux and temperature difference.
Minimum 0.5 W/m²/K.
"""
calc_heat_transfer_coefficient(convective_heat_flux, ΔT) =
    iszero(ΔT) ? 0.5u"W/m^2/K" : max(abs(convective_heat_flux / ΔT), 0.5u"W/m^2/K")

"""
    calc_mass_transfer_coefficient(heat_transfer_coefficient, air_specific_heat, air_density)

Mass transfer coefficient from heat transfer coefficient via the Lewis relation.
The (Pr/Sc)^0.666 = (0.71/0.60)^0.666 factor converts from heat to mass transfer.
"""
calc_mass_transfer_coefficient(heat_transfer_coefficient, air_specific_heat, air_density) =
    (heat_transfer_coefficient / (air_specific_heat * air_density)) * (0.71 / 0.60)^0.666

"""
    sublayer_stanton(z0, friction_velocity, kinematic_viscosity)

Stanton number for the roughness sublayer immediately above the surface,
`0.62·Re*^-0.45` with roughness Reynolds number `Re* = z0·u*/ν` (form
resembles Owen & Thomson (1963); citation not independently verified). This
is one way of expressing the momentum-to-scalar roughness (`z0`-to-`z0h`)
correction; see [`ScalarRoughnessRatio`](@ref) for microclimlearn's own
(fixed-ratio, `Re*`-independent) alternative, used instead of this at
canopy-top roughness scales.
"""
@inline function sublayer_stanton(z0, friction_velocity, kinematic_viscosity)
    roughness_reynolds_number = uconvert(NoUnits, z0 * friction_velocity / kinematic_viscosity)
    return 0.62 / roughness_reynolds_number^0.45
end

"""
    bulk_stanton(log_z_ratio)

Compute the bulk Stanton number for stable conditions.
"""
@inline function bulk_stanton(log_z_ratio)
    return 0.64 / log_z_ratio
end

"""
    bulk_stanton(log_z_ratio, z, obukhov_length)

Compute the bulk Stanton number for unstable conditions.
"""
@inline function bulk_stanton(log_z_ratio, z, obukhov_length)
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

Returns the Paulson (1970) x-substitution variable x = (1 - γ z / L)^(1/4),
used as the argument to `calc_ψ_m` and `calc_ψ_h`. Note this is the reciprocal
of the true Businger–Dyer φₘ stability function, φₘ = (1 - γ z / L)^(-1/4).

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
@inline function calc_φ_m(z, γ, obukhov_length)
    # Businger-Dyer/Paulson relations are only valid for |z/L| up to ~10;
    # unbounded, ζ→neutral can diverge and blow up friction_velocity/wind_speed.
    # Floor |L| at z/10, always negative (unstable) since only called from there.
    L_floor = z / 10.0
    L = obukhov_length < -L_floor ? obukhov_length : -L_floor
    # sqrt is faster than ^1/4
    return sqrt(sqrt(1.0 - γ * z / L))
end


"""
    calc_Φ_h(z, γ, obukhov_length)

True Businger–Dyer stability correction function for heat, Φ_h = x⁻², the
reciprocal square of `calc_φ_m`'s Paulson x-substitution (unlike `calc_φ_m`
itself, which returns `x`, not the momentum φₘ = x⁻¹).

# References
- Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971).
  Flux–profile relationships in the atmospheric surface layer.
  *Journal of the Atmospheric Sciences*, 28(2), 181–189.
"""
@inline function calc_Φ_h(z, γ, obukhov_length)
    return calc_φ_m(z, γ, obukhov_length)^(-2)
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
@inline function calc_ψ_m(x)
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
@inline function calc_ψ_h(x)
    return 2.0 * log((1.0 + x^2.0) / 2.0)
end

# Numerical floor on |obukhov_length|, sign-preserving -- avoids ζ=z/L blowing
# up ψ/Φ/bulk_stanton near L=0 (extreme instability *or* extreme stability,
# both of which are real, reachable states once the stable branch below is
# allowed to iterate down toward neutral). Same z/10 magnitude as calc_φ_m's
# own (unstable-only) floor, generalized to both signs. Not a physical
# parameter.
@inline function _floor_obukhov_length(z, obukhov_length)
    L_floor = z / 10.0
    return obukhov_length >= zero(obukhov_length) ? max(obukhov_length, L_floor) : min(obukhov_length, -L_floor)
end

"""
    calc_ψ_m(z, γ, obukhov_length, stable_beta)

General momentum stability correction ψ_m, valid for both unstable
(`obukhov_length < 0`, Businger–Dyer/Paulson form via [`calc_φ_m`](@ref)/
`calc_ψ_m(x)`) and stable (`obukhov_length ≥ 0`, Dyer (1974) linear form
`ψ_m = -stable_beta·z/L`) conditions -- the two meet continuously at
`obukhov_length → ±Inf` (neutral), so this is safe to call unconditionally
without pre-classifying the regime. Clamped to `[-4, 3]` in the stable
branch (same numerical bound Businger–Dyer's own unstable form is only
valid within, |ζ| ≲ 10).

# References
- Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971).
  Flux–profile relationships in the atmospheric surface layer.
  *Journal of the Atmospheric Sciences*, 28(2), 181–189.
- Dyer, A. J. (1974). A review of flux–profile relationships.
  *Boundary-Layer Meteorology*, 7(3), 363–372.
"""
@inline function calc_ψ_m(z, γ, obukhov_length, stable_beta)
    if obukhov_length < zero(obukhov_length)
        return calc_ψ_m(calc_φ_m(z, γ, obukhov_length))
    else
        ζ = z / _floor_obukhov_length(z, obukhov_length)
        return clamp(-stable_beta * ζ, -4.0, 3.0)
    end
end

"""
    calc_ψ_h(z, γ, obukhov_length, stable_beta, turbulent_prandtl_number)

General heat/scalar stability correction ψ_h, stable+unstable -- see
[`calc_ψ_m`](@ref). Stable form
`ψ_h = -(stable_beta/turbulent_prandtl_number)·z/L` (Businger et al. 1971's
turbulent Prandtl number relates the heat and momentum coefficients).
"""
@inline function calc_ψ_h(z, γ, obukhov_length, stable_beta, turbulent_prandtl_number)
    if obukhov_length < zero(obukhov_length)
        return calc_ψ_h(calc_φ_m(z, γ, obukhov_length))
    else
        ζ = z / _floor_obukhov_length(z, obukhov_length)
        return clamp(-(stable_beta * ζ) / turbulent_prandtl_number, -4.0, 3.0)
    end
end

"""
    calc_Φ_h(z, γ, obukhov_length, stable_Φ_h_coefficient, min_stable_Φ_h, max_stable_Φ_h)

General bulk heat stability multiplier Φ_h (eddy-diffusivity scaling in
[`KTheoryAirProfile`](@ref)/[`RaupachLTheoryAirProfile`](@ref)), stable+unstable
-- see [`calc_ψ_m`](@ref) for the (unrelated) log-profile analog. Stable form
`Φ_h = clamp(1 + stable_Φ_h_coefficient·ζ/(1+ζ), min_stable_Φ_h, max_stable_Φ_h)`,
matching R's `dphih` exactly -- a saturating rational form, not `calc_ψ_h`'s
linear one, since this feeds a multiplicative diffusivity scale rather than
an additive log-profile offset: unclamped linear growth here let very stable,
low-wind hours inflate `RaupachLTheoryAirProfile`'s `T_L` by orders of
magnitude (`_floor_obukhov_length` bounds `ζ` but not `Φ_h` itself).
"""
@inline function calc_Φ_h(z, γ, obukhov_length, stable_Φ_h_coefficient, min_stable_Φ_h, max_stable_Φ_h)
    if obukhov_length < zero(obukhov_length)
        return calc_Φ_h(z, γ, obukhov_length)
    else
        ζ = z / _floor_obukhov_length(z, obukhov_length)
        return clamp(1.0 + stable_Φ_h_coefficient * ζ / (1.0 + ζ), min_stable_Φ_h, max_stable_Φ_h)
    end
end


"""
    calc_Obukhov_length(reference_temp, surface_temp, v_ref_height, z0, z, ρcpTκg, κ, ΔT, ρ_cp, kinematic_viscosity;
                         γ=16.0, stable_beta=4.7, turbulent_prandtl_number=0.74,
                         max_iter=30, tol=1e-2, initial_obukhov_length=-0.3u"m")

`kinematic_viscosity` feeds `sublayer_stanton`'s roughness Reynolds number;
see that function's docstring for why it must be a real, temperature-dependent
value rather than a fixed constant.

Iteratively solve for the Monin–Obukhov length and convective heat flux, for
either sign of `ΔT` -- `convective_flux`'s sign already tracks `ΔT`'s sign
correctly either way (negative ΔT, surface warmer than reference → negative
heat flux → negative L, unstable; positive ΔT → positive L, stable), so no
pre-classification of the regime is needed before calling this: `calc_ψ_m`/
`calc_ψ_h` (used internally) handle whichever sign `obukhov_length` iterates
to. Converges to `±Inf` (neutral) as `ΔT → 0`.
"""
@inline function calc_Obukhov_length(
    reference_temp, surface_temp, v_ref_height, z0, z, ρcpTκg, κ, ΔT, ρ_cp, kinematic_viscosity;
    γ=16.0, stable_beta=4.7, turbulent_prandtl_number=0.74, max_iter=30, tol=1e-2, initial_obukhov_length=-0.3u"m",
    min_friction_velocity=0.02u"m/s", thermal_roughness_model=SublayerStantonRoughness(),
)
    obukhov_length = initial_obukhov_length

    # initialise with zeros
    convective_heat_flux = 0.0u"W/m^2"
    roughness_height_temperature = surface_temp
    ψ_h = 0.0
    friction_velocity = 0.0u"m/s"
    obukhov_length_new = 0.0u"m"

    relative_error = 1.0
    count = 0
    just_above_zero = 1.0e-6
    while relative_error > tol && count < max_iter
        count += 1
        ψ_m = calc_ψ_m(z, γ, obukhov_length, stable_beta)
        ψ_h = calc_ψ_h(z, γ, obukhov_length, stable_beta, turbulent_prandtl_number)
        # Floor at a fraction of the neutral log term, not an absolute
        # constant -- log(z/z0) can be O(1) for canopy geometry, where an
        # absolute 1e-6 floor lets friction_velocity blow up ~1e6-fold.
        log_z_ratio_i = log(z / z0)
        log_ratio_corrected = max(log_z_ratio_i - ψ_m, 0.1 * log_z_ratio_i, just_above_zero)
        friction_velocity = max(κ * v_ref_height / log_ratio_corrected, min_friction_velocity)
        convective_heat_flux, roughness_height_temperature = _unstable_convective_heat_flux(thermal_roughness_model,
            ρ_cp, ΔT, friction_velocity, κ, z0, z, log_ratio_corrected, ψ_h, obukhov_length, kinematic_viscosity,
            reference_temp, surface_temp)
        obukhov_length_new = ρcpTκg * friction_velocity^3 / convective_heat_flux
        relative_error = abs((obukhov_length_new - obukhov_length) / obukhov_length)
        obukhov_length = obukhov_length_new
    end

    return (;
        obukhov_length=u"m"(obukhov_length),
        friction_velocity,
        ψ_h,
        convective_heat_flux=u"W/m^2"(convective_heat_flux),
        roughness_height_temperature,
    )
end
