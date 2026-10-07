"""
    GarrattSegalCondensation()

Ground dewfall/frost from Garratt, J. R., & Segal, M. (1988). On the
contribution of atmospheric moisture to dew formation. Boundary-Layer
Meteorology, 45(3), 209-236. https://doi.org/10.1007/BF01066671, eq. 6b --
a Penman-style combination of an energy-availability term (net radiation
minus the surface's own actual current heat flux, weighted by the
saturation-specific-humidity-curve slope) and an aerodynamic term (the
*air's own* saturation deficit, scaled by a small neutral-log-law bulk
transfer coefficient and wind speed).

Unlike [`BulkTransferCondensation`](@ref), the energy term suppresses
condensation whenever the surface is net gaining energy (the usual daytime
case); the aerodynamic coefficient is typically an order of magnitude (or
more) smaller than the model's own Monin-Obukhov mass-transfer coefficient,
so magnitudes stay physically modest even when it does trigger.

Ground-only: the aerodynamic term's bulk transfer coefficient is a
bare-surface neutral log-law value, not meaningful for canopy leaves. Takes
`karman_constant` from `boundary_layer_model`.
"""
struct GarrattSegalCondensation <: AbstractCondensationModel end

function condensation_energy_flux(::GarrattSegalCondensation;
    evaporation_model, surface_temperature, air_temperature, relative_humidity, atmospheric_pressure,
    convective_heat_flux, mass_transfer_coefficient, actual_soil_wetness,
    wind_speed, reference_height, roughness_height, karman_constant,
    sky_temperature, surface_emissivity, absorbed_solar_radiation,
    vapour_pressure_equation=GoffGratch(), kw...,
)
    T_g = surface_temperature
    T_1 = air_temperature
    z = reference_height
    z0 = roughness_height
    ū = wind_speed

    surface_saturated_air = wet_air_properties(u"K"(T_g), 1.0, atmospheric_pressure; vapour_pressure_equation)
    ambient_air = wet_air_properties(u"K"(T_1), relative_humidity, atmospheric_pressure; vapour_pressure_equation)
    ambient_saturated_air = wet_air_properties(u"K"(T_1), 1.0, atmospheric_pressure; vapour_pressure_equation)

    # Specific humidities (kg water / kg moist air), eq. 5's q.
    q_g★ = surface_saturated_air.vapour_density / surface_saturated_air.density
    q_1  = ambient_air.vapour_density / ambient_air.density
    q_1★ = ambient_saturated_air.vapour_density / ambient_saturated_air.density
    δq = q_1★ - q_1  # air's own saturation deficit, eq. 6b

    ΔT = T_g - T_1
    ΔT = iszero(ΔT) ? 0.1u"K" : ΔT  # divide-by-zero guard, cf. DSUB.f:561
    s = (q_g★ - q_1★) / ΔT  # eq. 5

    λ = enthalpy_of_vaporisation(T_g)
    γ = surface_saturated_air.specific_heat / λ  # psychrometric constant

    net_longwave = σ * u"K"(sky_temperature)^4 - σ * surface_emissivity * u"K"(T_g)^4
    R_N = absorbed_solar_radiation + net_longwave

    # Actual (moisture-limited) evaporative loss -- distinct from the
    # always-saturated flux the aerodynamic term uses below.
    actual_evaporative_flux, _ = evaporation(evaporation_model;
        surface_temperature=T_g, air_temperature=T_1, relative_humidity,
        surface_relative_humidity=1.0, mass_transfer_coefficient, atmospheric_pressure,
        soil_wetness=actual_soil_wetness, saturated=false, vapour_pressure_equation,
    )
    G = R_N + convective_heat_flux - actual_evaporative_flux

    z_p = z0 / exp(2.0)  # eq. A10's "analogous scaling length"
    C_E = karman_constant^2 / (log(z / z0) * log(z / z_p))

    dew_frost_energy_flux = (s / (s + γ)) * (R_N - G) + (γ / (s + γ)) * (surface_saturated_air.density * λ * δq) * C_E * ū
    return dew_frost_energy_flux  # W/m², Monteith convention: +ve loss, -ve condensation
end
