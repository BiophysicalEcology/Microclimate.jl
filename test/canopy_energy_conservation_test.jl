using Microclimate
using Microclimate: allocate_interception, allocate_longwave, allocate_shortwave,
    canopy_interception!, canopy_longwave!, canopy_shortwave!
using Unitful
using Test

# Checks canopy_shortwave!/canopy_longwave! against global energy/mass
# conservation, using the buffer arrays directly. Calls the sub-models
# standalone (not through canopy_energy_balance!'s Picard loop) with a fixed
# leaf_temperature array.

const SIGMA_SB = 5.670374419e-8u"W/m^2/K^4"

site = Microclimate.example_site()
canopy_height = 5.0u"m"
plant_area_index = 3.0
n_layers = 10
canopy_projection_ratio = 1.0
leaf_emissivity = 0.97
ground_emissivity = 0.95

function make_environment_instant(; zenith_angle, reference_temperature = 293.0u"K")
    (;
        atmospheric_pressure = 101325.0u"Pa", reference_temperature, reference_wind_speed = 2.0u"m/s",
        reference_humidity = 0.5, zenith_angle, cloud_emissivity = 1.0, cloud_cover = 0.3,
        surface_emissivity = 0.95, shade = 0.0,
    )
end

conditions = [
    ("midday, clear", 20.0u"°", 600.0u"W/m^2", 100.0u"W/m^2", 295.0u"K"),
    ("overcast", 60.0u"°", 0.0u"W/m^2", 150.0u"W/m^2", 290.0u"K"),
    ("night", 100.0u"°", 0.0u"W/m^2", 0.0u"W/m^2", 283.0u"K"),
]

@testset "canopy_shortwave! conservation" begin
    shortwave_model = TwoStreamRadiation()
    buffers = allocate_shortwave(shortwave_model, canopy_height, plant_area_index, n_layers, canopy_projection_ratio)

    for (label, zenith_angle, direct, diffuse, _) in conditions
        @testset "$label" begin
            result = canopy_shortwave!(buffers, shortwave_model, plant_area_index;
                zenith_angle, direct_horizontal_irradiance = direct, diffuse_horizontal_irradiance = diffuse,
                ground_reflectance = 0.15)

            # net_absorbed_shortwave already includes ground absorption
            # (everything not reflected back to sky) -- not leaf-only.
            reflected_to_sky = buffers.boundary_upward_shortwave[1]
            # rtol: direct-beam's particular solution doesn't exactly satisfy
            # the ground boundary condition (up to ~3e-3 at high direct beam).
            @test ustrip(u"W/m^2", direct + diffuse) ≈
                ustrip(u"W/m^2", result.net_absorbed_shortwave + reflected_to_sky) rtol=1e-3
            @test ustrip(u"W/m^2", sum(buffers.absorbed_shortwave) + result.ground_absorbed_shortwave) ≈
                ustrip(u"W/m^2", result.net_absorbed_shortwave) rtol=3e-3
        end
    end
end

@testset "canopy_longwave! conservation" begin
    longwave_model = LayeredLongwaveExchange()
    buffers = allocate_longwave(longwave_model, plant_area_index, n_layers, canopy_projection_ratio)

    for (label, zenith_angle, _, _, ground_T) in conditions
        @testset "$label" begin
            environment_instant = make_environment_instant(; zenith_angle)
            leaf_temperature = range(ground_T - 5.0u"K", environment_instant.reference_temperature; length = n_layers)

            result = canopy_longwave!(buffers, longwave_model, leaf_emissivity;
                leaf_temperature, ground_temperature = ground_T, ground_emissivity, site, environment_instant)

            ground_emission = ground_emissivity * SIGMA_SB * u"K"(ground_T)^4
            # Each layer radiates emission*(1-τ) both up and down -- 2x.
            leaf_emission_total = 2.0 * sum(
                leaf_emissivity * SIGMA_SB * u"K"(T)^4 * (1.0 - τ)
                for (T, τ) in zip(leaf_temperature, buffers.layer_transmission)
            )
            sky_incoming = buffers.boundary_downward_longwave[1]
            escaped_to_sky = buffers.boundary_upward_longwave[1]
            total_absorbed = result.ground_absorbed_longwave + sum(buffers.absorbed_longwave)

            @test ustrip(u"W/m^2", sky_incoming + ground_emission + leaf_emission_total) ≈
                ustrip(u"W/m^2", total_absorbed + escaped_to_sky) atol=1e-6
        end
    end
end

@testset "_cap_wet_surface_evaporation conserves mass and energy" begin
    air_T = 293.0u"K"
    flux = 1.0e-4u"kg/s"  # demands 0.36 kg/m^2 over the hour at wet_fraction=1
    sensible = 50.0u"W/m^2"

    @testset "ample storage: passes through unchanged" begin
        flux2, sensible2, actual = Microclimate._cap_wet_surface_evaporation(flux, sensible, 1.0, 1.0u"kg/m^2", air_T)
        @test flux2 == flux
        @test sensible2 == sensible
        @test ustrip(u"kg/m^2", actual) ≈
            ustrip(u"kg/m^2", uconvert(u"kg", flux * Microclimate.CANOPY_TIMESTEP) / 1.0u"m^2")
    end

    @testset "scarce storage: withdrawal capped, shortfall moved to sensible" begin
        available = 0.05u"kg/m^2"
        flux2, sensible2, actual = Microclimate._cap_wet_surface_evaporation(flux, sensible, 1.0, available, air_T)
        @test actual == available
        @test flux2 < flux
        latent_removed = uconvert(u"W", (flux - flux2) * Microclimate.enthalpy_of_vaporisation(air_T))
        sensible_added = (sensible2 - sensible) * 1.0u"m^2"
        @test ustrip(u"W", latent_removed) ≈ ustrip(u"W", sensible_added) rtol=1e-10
        @test ustrip(u"kg/m^2", actual) ≈
            ustrip(u"kg/m^2", uconvert(u"kg", flux2 * Microclimate.CANOPY_TIMESTEP) / 1.0u"m^2") atol=1e-12
    end

    @testset "no storage: wet-surface share fully rerouted, dry (transpiration) share untouched" begin
        wet_fraction = 0.5
        flux2, sensible2, actual = Microclimate._cap_wet_surface_evaporation(flux, sensible, wet_fraction, 0.0u"kg/m^2", air_T)
        @test actual == 0.0u"kg/m^2"
        @test ustrip(u"kg/s", flux2) ≈ ustrip(u"kg/s", flux * (1.0 - wet_fraction)) rtol=1e-10
    end
end

@testset "rain interception mass conservation" begin
    model = LayeredRainInterception(; leaf_water_storage_capacity = 0.2u"kg/m^2")
    boundary_layer_model = MoninObukhov()
    heights = vcat(0.5:0.5:5.0, [7.0]) .* u"m"
    buffers = allocate_interception(model, canopy_height, plant_area_index, heights, n_layers, boundary_layer_model)
    wind_speed = fill(1.0u"m/s", n_layers)

    for rainfall in (0.0u"kg/m^2", 0.5u"kg/m^2", 5.0u"kg/m^2", 50.0u"kg/m^2")
        stored_before = sum(buffers.leaf_surface_water)
        result = canopy_interception!(buffers, model, 1.0; rainfall, wind_speed)
        stored_after = sum(buffers.leaf_surface_water)
        @test ustrip(u"kg/m^2", rainfall) ≈
            ustrip(u"kg/m^2", result.ground_throughfall + (stored_after - stored_before)) atol=1e-9
    end
end

@testset "canopy_energy_balance!: wet-surface evaporation can't exceed stored water" begin
    heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"
    boundary_layer_model = MoninObukhov()
    environment_instant = make_environment_instant(; zenith_angle = 20.0u"°")

    run_with_capacity(capacity) = begin
        model = Microclimate.example_multilayer_canopy(;
            canopy_height = 1.0u"m", plant_area_index = 3.0,
            interception_model = LayeredRainInterception(; leaf_water_storage_capacity = capacity),
        )
        buffers = Microclimate.allocate_canopy(model, heights, boundary_layer_model)
        inputs = Microclimate.CanopyEnergyBalanceInputs(model;
            site, environment_instant, zenith_angle = 20.0u"°",
            direct_horizontal_irradiance = 600.0u"W/m^2", diffuse_horizontal_irradiance = 100.0u"W/m^2",
            ground_reflectance = 0.15, ground_temperature = 295.0u"K", ground_emissivity = 0.95,
            ground_relative_humidity = 0.5, canopy_source_temperature = 295.0u"K", rainfall = 5.0u"kg/m^2",
        )
        Microclimate.canopy_energy_balance!(buffers, model, boundary_layer_model, inputs)
        buffers
    end

    scarce = run_with_capacity(1.0e-4u"kg/m^2")
    abundant = run_with_capacity(10.0u"kg/m^2")

    @test all(>=(0.0u"kg/m^2"), scarce.interception.leaf_surface_water)
    @test all(isfinite, ustrip.(u"kg/s", scarce.leaf.evaporation_mass_flow))
    @test all(isfinite, ustrip.(u"W/m^2", scarce.leaf.sensible_heat_source))
    # scarce storage forces some latent demand into sensible heat instead
    @test sum(scarce.leaf.evaporation_mass_flow) < sum(abundant.leaf.evaporation_mass_flow)
    @test sum(scarce.leaf.sensible_heat_source) > sum(abundant.leaf.sensible_heat_source)
end
