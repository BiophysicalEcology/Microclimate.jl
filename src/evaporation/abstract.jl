abstract type AbstractEvaporationModel end

"""
    evaporation(model; surface_temperature, air_temperature, relative_humidity, surface_relative_humidity, mass_transfer_coefficient, atmospheric_pressure, soil_wetness, saturated, vapour_pressure_equation)

Surface latent flux: returns `(Q_evaporation, evaporation_mass_flux)` in W/m²
and g/s/m². Variants supply the parameterisation.
"""
function evaporation end

"""
    surface_convection_evaporation(model; boundary_layer_model, surface_temperature, air_temperature, wind_speed, relative_humidity, atmospheric_pressure, roughness_height, reference_height, zenith_angle, soil_wetness, vapour_pressure_equation, obukhov_length_prev)

Combined entry point: picks stable vs unstable atmosphere, computes friction
velocity + convective heat flux via Monin–Obukhov, then calls `evaporation`
for the latent flux. Returns `(Q_evaporation, convective_heat_flux)`.
"""
function surface_convection_evaporation end
