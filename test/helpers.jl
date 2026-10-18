# Shared test helpers

# Bytes allocated by a warmed-up call, measured in function scope so
# non-constant test globals don't add their own boxing
function allocations(f::F, args::Vararg{Any,N}; kw...) where {F,N}
    f(args...; kw...)
    return @allocated f(args...; kw...)
end

environment_instant(; zenith_angle = 30.0u"°", reference_temperature = 293.0u"K") = (;
    atmospheric_pressure = 101325.0u"Pa", reference_temperature, reference_wind_speed = 2.0u"m/s",
    reference_humidity = 0.5, zenith_angle, cloud_emissivity = 1.0, cloud_cover = 0.3,
    surface_emissivity = 0.95, shade = 0.0,
)

# Hourly canopy inputs for a clear midday hour over a 290 K ground
canopy_inputs(model; zenith_angle = 30.0u"°", reference_temperature = 293.0u"K",
    direct_horizontal_irradiance = 500.0u"W/m^2", diffuse_horizontal_irradiance = 100.0u"W/m^2",
    ground_temperature = 290.0u"K", canopy_source_temperature = reference_temperature, kw...) =
    Microclimate.CanopyEnergyBalanceInputs(model;
        environment_instant = environment_instant(; zenith_angle, reference_temperature), zenith_angle,
        direct_horizontal_irradiance, diffuse_horizontal_irradiance,
        ground_reflectance = 0.15, ground_temperature, ground_emissivity = 0.95,
        ground_relative_humidity = 0.5, canopy_source_temperature, kw...,
    )
