"""
    BulkTransferEvaporation()

Bulk transfer evaporation surface evaporation (as implemented in NicheMapR): Penman-style mass-transfer
flux from surface and air vapour densities times wetness, with a cold-surface
mass-loss clamp. The combined `surface_convection_evaporation` picks the
stable / unstable atmospheric branch via Monin–Obukhov.
"""
struct BulkTransferEvaporation <: AbstractEvaporationModel end

function evaporation(::BulkTransferEvaporation;
    surface_temperature,
    air_temperature,
    relative_humidity,
    surface_relative_humidity,
    mass_transfer_coefficient,
    atmospheric_pressure,
    soil_wetness,
    saturated,
    vapour_pressure_equation=GoffGratch(),
)
    # Assumes all units are SI (Kelvin, Pascal, meters, seconds, kg, etc.)

    clamped_surface_temperature = surface_temperature < u"K"(-81.0u"°C") ? u"K"(-81.0u"°C") : surface_temperature

    # surface and air vapor densities
    surface_vapour_density = wet_air_properties(u"K"(clamped_surface_temperature), surface_relative_humidity, atmospheric_pressure; vapour_pressure_equation).vapour_density
    air_vapour_density = wet_air_properties(u"K"(air_temperature), relative_humidity, atmospheric_pressure; vapour_pressure_equation).vapour_density

    # Effective surface wetness fraction
    surface_wetness = saturated ? 1.0 : soil_wetness

    # Water evaporated from surface (kg/s/m^2)
    water_flux = surface_wetness * mass_transfer_coefficient * (surface_vapour_density - air_vapour_density)

    # Latent heat of vaporization (J/kg)
    latent_heat_vaporisation = enthalpy_of_vaporisation(clamped_surface_temperature)

    # Energy flux due to evaporation (W/m²). uconvert keeps the return type
    # stable as W/m² so callers can mix it with other W/m² fluxes without
    # triggering Unitful's runtime unit-promotion heap allocation.
    Q_evaporation = uconvert(u"W/m^2", water_flux * latent_heat_vaporisation)

    # Mass flux (g/s/m²)
    evaporation_mass_flux = u"g/s/m^2"(water_flux)

    # No water loss if surface temperature ≤ 0°C (e.g., melting snow only)
    if surface_temperature <= u"K"(0.0u"°C")
        evaporation_mass_flux = 0.0u"g/s/m^2"
    end

    return Q_evaporation, evaporation_mass_flux
end

function surface_convection_evaporation(model::BulkTransferEvaporation;
    boundary_layer_model,
    surface_temperature, air_temperature, wind_speed, relative_humidity,
    atmospheric_pressure, roughness_height, reference_height,
    zenith_angle=90.0u"°", soil_wetness, vapour_pressure_equation=GoffGratch(),
    obukhov_length_prev=nothing,
)
    (; convective_heat_flux, ΔT) = surface_fluxes(boundary_layer_model;
        surface_temperature, air_temperature, wind_speed, zenith_angle,
        roughness_height, reference_height, atmospheric_pressure, obukhov_length_prev,
    )
    heat_transfer_coefficient = calc_heat_transfer_coefficient(convective_heat_flux, ΔT)
    wet_air_out = wet_air_properties(u"K"(air_temperature), relative_humidity, atmospheric_pressure; vapour_pressure_equation)
    mass_transfer_coefficient = calc_mass_transfer_coefficient(heat_transfer_coefficient, wet_air_out.specific_heat, wet_air_out.density)
    Q_evaporation, _ = evaporation(model;
        surface_temperature=u"K"(surface_temperature),
        air_temperature=u"K"(air_temperature),
        relative_humidity,
        surface_relative_humidity=1.0,
        mass_transfer_coefficient,
        atmospheric_pressure,
        soil_wetness,
        saturated=false,
        vapour_pressure_equation,
    )
    return Q_evaporation, convective_heat_flux
end
