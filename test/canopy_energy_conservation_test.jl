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
            # the ground boundary condition (documented, ~6e-4 max error).
            @test ustrip(u"W/m^2", direct + diffuse) ≈
                ustrip(u"W/m^2", result.net_absorbed_shortwave + reflected_to_sky) rtol=1e-3
            @test ustrip(u"W/m^2", sum(buffers.absorbed_shortwave) + result.ground_absorbed_shortwave) ≈
                ustrip(u"W/m^2", result.net_absorbed_shortwave) rtol=1e-3
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
