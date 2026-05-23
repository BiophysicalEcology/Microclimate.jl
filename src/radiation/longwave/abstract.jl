abstract type AbstractLongwaveScheme end

"""
    longwave_radiation(scheme; site, environment_instant, surface_temperature, vapour_pressure_equation)

Full surface longwave budget: net longwave, sky components, and ground term at
the given surface temperature. Sub-models (e.g. the atmospheric clear-sky
emissivity model) live on the longwave `scheme` itself.
"""
function longwave_radiation end
