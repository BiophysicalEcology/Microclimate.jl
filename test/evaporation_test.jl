using Microclimate
using Microclimate: evaporation
using Unitful
using Test

model = BulkTransferEvaporation()
common = (;
    atmospheric_pressure = 101325.0u"Pa",
    mass_transfer_coefficient = 0.01u"m/s",
    soil_wetness = 1.0,
    saturated = true,
)

@testset "condensation (dew/frost) is not zeroed at or below 0°C" begin
    # Cold, saturated surface; warm, humid air -- condenses.
    Q, mass_flux = evaporation(model; common...,
        surface_temperature = u"K"(-5.0u"°C"),
        surface_relative_humidity = 1.0,
        air_temperature = u"K"(5.0u"°C"),
        relative_humidity = 0.9,
    )
    @test mass_flux < 0.0u"g/s/m^2"  # condensing, i.e. dew/frost forming
    @test Q < 0.0u"W/m^2"            # latent heat released to the surface
end

@testset "evaporative loss is still zeroed at or below 0°C" begin
    # Cold, saturated surface; equally cold, dry air -- would evaporate, but
    # the cold-surface clamp still suppresses mass loss (e.g. no phantom
    # sublimation of nonexistent frost from bare, dry ground).
    Q, mass_flux = evaporation(model; common...,
        surface_temperature = u"K"(-5.0u"°C"),
        surface_relative_humidity = 1.0,
        air_temperature = u"K"(-5.0u"°C"),
        relative_humidity = 0.1,
    )
    @test mass_flux == 0.0u"g/s/m^2"
end

@testset "evaporative loss above 0°C is unaffected" begin
    Q, mass_flux = evaporation(model; common...,
        surface_temperature = u"K"(20.0u"°C"),
        surface_relative_humidity = 1.0,
        air_temperature = u"K"(20.0u"°C"),
        relative_humidity = 0.1,
    )
    @test mass_flux > 0.0u"g/s/m^2"
end
