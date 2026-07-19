using Microclimate
using Unitful
using Test

site = Microclimate.example_site()
environment_instant = (;
    atmospheric_pressure = 101325.0u"Pa",
    reference_humidity = 0.5,
    reference_temperature = 293.0u"K",
    cloud_emissivity = 1.0,
    cloud_cover = 0.3,
    surface_emissivity = 0.95,
    shade = 0.0,
)

model = LayeredLongwaveExchange()
plant_area_index = 3.0
n_layers = 10
buffers = Microclimate.allocate_longwave(model, plant_area_index, n_layers)
leaf_emissivity = 0.97

@testset "energy conservation: net flux divergence closes to floating-point precision" begin
    leaf_temperature = fill(288.0u"K", n_layers)
    ground_temperature = 285.0u"K"
    ground_emissivity = 0.95

    result = Microclimate.canopy_longwave!(buffers, model, leaf_emissivity;
        leaf_temperature, ground_temperature, ground_emissivity, site, environment_instant)

    net_layers = sum(1:n_layers) do i
        buffers.absorbed_longwave[i] - 2.0 * leaf_emissivity * Microclimate.σ * u"K"(leaf_temperature[i])^4 * (1.0 - buffers.layer_transmission[i])
    end
    net_ground = result.ground_absorbed_longwave - ground_emissivity * Microclimate.σ * u"K"(ground_temperature)^4

    @test ustrip(u"W/m^2", net_layers + net_ground) ≈ ustrip(u"W/m^2", result.canopy_absorbed_longwave) atol=1e-8
end

@testset "warm sky, cold canopy: downward flux monotonically non-increasing, net gain" begin
    leaf_temperature = fill(250.0u"K", n_layers)  # much colder than the ~330+ K sky implied below
    ground_temperature = 250.0u"K"
    ground_emissivity = 0.95
    hot_environment = merge(environment_instant, (; reference_temperature = 320.0u"K", cloud_cover = 0.0))

    result = Microclimate.canopy_longwave!(buffers, model, leaf_emissivity;
        leaf_temperature, ground_temperature, ground_emissivity, site, environment_instant=hot_environment)

    @test issorted(ustrip.(u"W/m^2", buffers.boundary_downward_longwave), rev=true)
    @test result.canopy_absorbed_longwave > 0.0u"W/m^2"
    @test all(>(0.0u"W/m^2"), buffers.absorbed_longwave)
end

@testset "isothermal canopy at sky-equivalent temperature: near-zero net exchange" begin
    # If leaf/ground temperature equals the effective sky temperature AND both are
    # perfect emitters (emissivity 1), the whole profile is in radiative
    # equilibrium (`sky_temperature` is defined as `(incoming_longwave/σ)^(1/4)`,
    # i.e. exactly the blackbody temperature matching the incoming flux). A leaf
    # emissivity below 1 (checked elsewhere) breaks this deliberately: an
    # imperfect emitter at the sky-equivalent temperature radiates less than a
    # blackbody would, so a small systematic net gain is real, not a bug.
    sky = Microclimate.precompute_longwave_sky(model.atmospheric_radiation_model;
        site, environment_instant, shade=0.0)
    T_eq = sky.sky_temperature
    leaf_temperature = fill(T_eq, n_layers)

    result = Microclimate.canopy_longwave!(buffers, model, 1.0;
        leaf_temperature, ground_temperature=T_eq, ground_emissivity=1.0, site, environment_instant)

    @test abs(ustrip(u"W/m^2", result.canopy_absorbed_longwave)) < 1e-6
end
