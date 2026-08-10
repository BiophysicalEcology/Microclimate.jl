using Microclimate
using Microclimate: allocate_longwave, canopy_longwave!, precompute_longwave_sky
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
canopy_projection_ratio = 1.0
buffers = Microclimate.allocate_longwave(model, plant_area_index, n_layers, canopy_projection_ratio)
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

    @test ustrip(u"W/m^2", net_layers + net_ground) ≈ ustrip(u"W/m^2", result.net_absorbed_longwave) atol=1e-8
end

@testset "warm sky, cold canopy: downward flux monotonically non-increasing, net gain" begin
    leaf_temperature = fill(250.0u"K", n_layers)  # much colder than the ~330+ K sky implied below
    ground_temperature = 250.0u"K"
    ground_emissivity = 0.95
    hot_environment = merge(environment_instant, (; reference_temperature = 320.0u"K", cloud_cover = 0.0))

    result = Microclimate.canopy_longwave!(buffers, model, leaf_emissivity;
        leaf_temperature, ground_temperature, ground_emissivity, site, environment_instant=hot_environment)

    @test issorted(ustrip.(u"W/m^2", buffers.boundary_downward_longwave), rev=true)
    @test result.net_absorbed_longwave > 0.0u"W/m^2"
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

    @test abs(ustrip(u"W/m^2", result.net_absorbed_longwave)) < 1e-6
end

@testset "AllPairsLongwaveExchange: runs and gives physically sane fluxes" begin
    # Ported from micropoint/microclimlearn for direct comparison against that
    # reference (see module_comparison.jl); not held to the same conservation
    # bar as the two models below since it inherits the source's own
    # simplifications (no ground reflection, no self-view term).
    ap_model = AllPairsLongwaveExchange()
    ap_buffers = Microclimate.allocate_longwave(ap_model, plant_area_index, n_layers, canopy_projection_ratio)
    leaf_temperature = fill(288.0u"K", n_layers)
    ground_temperature = 285.0u"K"

    result = Microclimate.canopy_longwave!(ap_buffers, ap_model, leaf_emissivity;
        leaf_temperature, ground_temperature, ground_emissivity=0.95, site, environment_instant)

    @test all(>(0.0u"W/m^2"), ap_buffers.downward_longwave)
    @test all(>(0.0u"W/m^2"), ap_buffers.upward_longwave)
    @test result.ground_absorbed_longwave > 0.0u"W/m^2"
end

@testset "LayeredRadiosityExchange" begin
    lr_model = LayeredRadiosityExchange()
    lr_buffers = Microclimate.allocate_longwave(lr_model, plant_area_index, n_layers, canopy_projection_ratio)

    @testset "isothermal blackbody equilibrium: zero net exchange everywhere" begin
        sky = Microclimate.precompute_longwave_sky(lr_model.atmospheric_radiation_model;
            site, environment_instant, shade=0.0)
        T_eq = sky.sky_temperature
        leaf_temperature = fill(T_eq, n_layers)

        result = Microclimate.canopy_longwave!(lr_buffers, lr_model, 1.0;
            leaf_temperature, ground_temperature=T_eq, ground_emissivity=1.0, site, environment_instant)

        @test abs(ustrip(u"W/m^2", result.net_absorbed_longwave)) < 1e-6
        @test all(f -> abs(ustrip(u"W/m^2", f)) < 1e-6, lr_buffers.absorbed_longwave .- 2.0 .* Microclimate.σ .* T_eq^4 .* (1.0 .- lr_buffers.layer_transmission))
    end

    @testset "isothermal grey-body equilibrium: reflection doesn't break Kirchhoff's law" begin
        # Different emissivities at a UNIFORM temperature must still give zero net
        # exchange everywhere -- this is what catches a reflection split that
        # isn't energy-conserving (e.g. an asymmetric backscatter fraction without
        # compensating forward-scatter terms).
        sky = Microclimate.precompute_longwave_sky(lr_model.atmospheric_radiation_model;
            site, environment_instant, shade=0.0)
        T_eq = sky.sky_temperature
        leaf_temperature = fill(T_eq, n_layers)

        result = Microclimate.canopy_longwave!(lr_buffers, lr_model, 0.92;
            leaf_temperature, ground_temperature=T_eq, ground_emissivity=0.9, site, environment_instant)

        @test abs(ustrip(u"W/m^2", result.net_absorbed_longwave)) < 1e-6
        @test all(f -> abs(ustrip(u"W/m^2", f)) < 1e-6, lr_buffers.boundary_upward_longwave .- Microclimate.σ .* T_eq^4)
        @test all(f -> abs(ustrip(u"W/m^2", f)) < 1e-6, lr_buffers.boundary_downward_longwave .- Microclimate.σ .* T_eq^4)
    end

    @testset "global energy conservation, non-isothermal, emissivity<1" begin
        leaf_temperature = range(280.0u"K", 293.0u"K"; length=n_layers)
        ground_temperature = 281.0u"K"
        ground_emissivity = 0.95

        sky = Microclimate.precompute_longwave_sky(lr_model.atmospheric_radiation_model;
            site, environment_instant, shade=0.0)

        result = Microclimate.canopy_longwave!(lr_buffers, lr_model, leaf_emissivity;
            leaf_temperature, ground_temperature, ground_emissivity, site, environment_instant)

        total_leaf_emission = 2.0 * sum(1:n_layers) do i
            leaf_emissivity * (1.0 - lr_buffers.layer_transmission[i]) * Microclimate.σ * u"K"(leaf_temperature[i])^4
        end
        ground_emission = ground_emissivity * Microclimate.σ * u"K"(ground_temperature)^4
        escaping_to_sky = sky.incoming_longwave - result.net_absorbed_longwave

        total_in = sky.incoming_longwave + total_leaf_emission + ground_emission
        total_out = sum(lr_buffers.absorbed_longwave) + result.ground_absorbed_longwave + escaping_to_sky

        @test ustrip(u"W/m^2", total_in) ≈ ustrip(u"W/m^2", total_out) rtol=1e-6
    end

    @testset "black-leaf limit (emissivity=1) matches LayeredLongwaveExchange exactly" begin
        leaf_temperature = range(280.0u"K", 293.0u"K"; length=n_layers)
        ground_temperature = 281.0u"K"
        ground_emissivity = 0.95

        result_radiosity = Microclimate.canopy_longwave!(lr_buffers, lr_model, 1.0;
            leaf_temperature, ground_temperature, ground_emissivity, site, environment_instant)
        result_sequential = Microclimate.canopy_longwave!(buffers, model, 1.0;
            leaf_temperature, ground_temperature, ground_emissivity, site, environment_instant)

        @test ustrip.(u"W/m^2", lr_buffers.boundary_downward_longwave) ≈
              ustrip.(u"W/m^2", buffers.boundary_downward_longwave) atol=1e-9
        @test ustrip.(u"W/m^2", lr_buffers.boundary_upward_longwave) ≈
              ustrip.(u"W/m^2", buffers.boundary_upward_longwave) atol=1e-9
        @test ustrip(u"W/m^2", result_radiosity.ground_absorbed_longwave) ≈
              ustrip(u"W/m^2", result_sequential.ground_absorbed_longwave) atol=1e-9
    end
end
