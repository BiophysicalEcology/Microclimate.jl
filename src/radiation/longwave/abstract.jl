abstract type AbstractLongwaveScheme end

"""
    longwave_radiation(scheme, atmospheric_radiation_model; site, environment_instant, surface_temperature, vapour_pressure_equation)

Full surface longwave budget: net longwave, sky components, and ground term at
the given surface temperature.
"""
function longwave_radiation end
