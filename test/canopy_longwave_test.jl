using Microclimate
using Microclimate: allocate_longwave, canopy_longwave!, canopy_incoming_longwave, precompute_longwave_sky
using Unitful
using Test

include("helpers.jl")

site = Microclimate.example_site()
environment = environment_instant()

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
        leaf_temperature, ground=bare_ground_surface(; temperature=ground_temperature, emissivity=ground_emissivity),
        incoming_longwave=canopy_incoming_longwave(model; site, environment_instant=environment))

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
    hot_environment = merge(environment, (; reference_temperature = 320.0u"K", cloud_cover = 0.0))

    result = Microclimate.canopy_longwave!(buffers, model, leaf_emissivity;
        leaf_temperature, ground=bare_ground_surface(; temperature=ground_temperature, emissivity=ground_emissivity),
        incoming_longwave=canopy_incoming_longwave(model; site, environment_instant=hot_environment))

    @test issorted(ustrip.(u"W/m^2", buffers.boundary_downward_longwave), rev=true)
    @test result.net_absorbed_longwave > 0.0u"W/m^2"
    @test all(>(0.0u"W/m^2"), buffers.absorbed_longwave)
end

@testset "isothermal canopy at sky-equivalent temperature: near-zero net exchange" begin
    sky = Microclimate.precompute_longwave_sky(model.atmospheric_radiation_model;
        site, environment_instant=environment, shade=0.0)
    T_eq = sky.sky_temperature
    leaf_temperature = fill(T_eq, n_layers)

    for emissivity in (1.0, leaf_emissivity, 0.9)
        result = Microclimate.canopy_longwave!(buffers, model, emissivity;
            leaf_temperature, ground=bare_ground_surface(; temperature=T_eq, emissivity=1.0),
            incoming_longwave=canopy_incoming_longwave(model; site, environment_instant=environment))

        @test abs(ustrip(u"W/m^2", result.net_absorbed_longwave)) < 1e-6
        emitted = 2.0 .* emissivity .* Microclimate.σ .* T_eq^4 .* (1.0 .- buffers.layer_transmission)
        @test all(f -> abs(ustrip(u"W/m^2", f)) < 1e-6, buffers.absorbed_longwave .- emitted)
    end
end

@testset "LayeredRadiosityExchange" begin
    lr_model = LayeredRadiosityExchange()
    lr_buffers = Microclimate.allocate_longwave(lr_model, plant_area_index, n_layers, canopy_projection_ratio)

    @testset "isothermal blackbody equilibrium: zero net exchange everywhere" begin
        sky = Microclimate.precompute_longwave_sky(lr_model.atmospheric_radiation_model;
            site, environment_instant=environment, shade=0.0)
        T_eq = sky.sky_temperature
        leaf_temperature = fill(T_eq, n_layers)

        result = Microclimate.canopy_longwave!(lr_buffers, lr_model, 1.0;
            leaf_temperature, ground=bare_ground_surface(; temperature=T_eq, emissivity=1.0),
            incoming_longwave=canopy_incoming_longwave(lr_model; site, environment_instant=environment))

        @test abs(ustrip(u"W/m^2", result.net_absorbed_longwave)) < 1e-6
        @test all(f -> abs(ustrip(u"W/m^2", f)) < 1e-6, lr_buffers.absorbed_longwave .- 2.0 .* Microclimate.σ .* T_eq^4 .* (1.0 .- lr_buffers.layer_transmission))
    end

    @testset "isothermal grey-body equilibrium: reflection doesn't break Kirchhoff's law" begin
        # Different emissivities at a UNIFORM temperature must still give zero net
        # exchange everywhere -- this is what catches a reflection split that
        # isn't energy-conserving (e.g. an asymmetric backscatter fraction without
        # compensating forward-scatter terms).
        sky = Microclimate.precompute_longwave_sky(lr_model.atmospheric_radiation_model;
            site, environment_instant=environment, shade=0.0)
        T_eq = sky.sky_temperature
        leaf_temperature = fill(T_eq, n_layers)

        result = Microclimate.canopy_longwave!(lr_buffers, lr_model, 0.92;
            leaf_temperature, ground=bare_ground_surface(; temperature=T_eq, emissivity=0.9, longwave_reflectance=0.1),
            incoming_longwave=canopy_incoming_longwave(lr_model; site, environment_instant=environment))

        @test abs(ustrip(u"W/m^2", result.net_absorbed_longwave)) < 1e-6
        @test all(f -> abs(ustrip(u"W/m^2", f)) < 1e-6, lr_buffers.boundary_upward_longwave .- Microclimate.σ .* T_eq^4)
        @test all(f -> abs(ustrip(u"W/m^2", f)) < 1e-6, lr_buffers.boundary_downward_longwave .- Microclimate.σ .* T_eq^4)
    end

    @testset "global energy conservation, non-isothermal, emissivity<1" begin
        leaf_temperature = range(280.0u"K", 293.0u"K"; length=n_layers)
        ground_temperature = 281.0u"K"
        ground_emissivity = 0.95

        sky = Microclimate.precompute_longwave_sky(lr_model.atmospheric_radiation_model;
            site, environment_instant=environment, shade=0.0)

        result = Microclimate.canopy_longwave!(lr_buffers, lr_model, leaf_emissivity;
            leaf_temperature, ground=bare_ground_surface(; temperature=ground_temperature, emissivity=ground_emissivity),
            incoming_longwave=canopy_incoming_longwave(lr_model; site, environment_instant=environment))

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

        ground = bare_ground_surface(; temperature=ground_temperature, emissivity=ground_emissivity)
        result_radiosity = Microclimate.canopy_longwave!(lr_buffers, lr_model, 1.0;
            leaf_temperature, ground, incoming_longwave=canopy_incoming_longwave(lr_model; site, environment_instant=environment))
        result_sequential = Microclimate.canopy_longwave!(buffers, model, 1.0;
            leaf_temperature, ground, incoming_longwave=canopy_incoming_longwave(model; site, environment_instant=environment))

        @test ustrip.(u"W/m^2", lr_buffers.boundary_downward_longwave) ≈
              ustrip.(u"W/m^2", buffers.boundary_downward_longwave) atol=1e-9
        @test ustrip.(u"W/m^2", lr_buffers.boundary_upward_longwave) ≈
              ustrip.(u"W/m^2", buffers.boundary_upward_longwave) atol=1e-9
        @test ustrip(u"W/m^2", result_radiosity.ground_absorbed_longwave) ≈
              ustrip(u"W/m^2", result_sequential.ground_absorbed_longwave) atol=1e-9
    end
end

@testset "the longwave at the ground is linear in what the ground sends up" begin
    leaf_temperature = collect(range(280.0u"K", 293.0u"K"; length=n_layers))
    for longwave_model in (LayeredLongwaveExchange(), LayeredRadiosityExchange(), BlackBodyLongwaveExchange()), emissivity in (0.9, 0.97)
        longwave_buffers = Microclimate.allocate_longwave(longwave_model, plant_area_index, n_layers, canopy_projection_ratio)
        incoming_longwave = canopy_incoming_longwave(longwave_model; site, environment_instant=environment)
        (; downward_longwave, longwave_return_fraction) = Microclimate.ground_longwave_conditions(longwave_buffers,
            longwave_model, emissivity; leaf_temperature, incoming_longwave)
        result = Microclimate.canopy_longwave!(longwave_buffers, longwave_model, emissivity;
            leaf_temperature, ground=bare_ground_surface(; temperature=285.0u"K", emissivity=0.95), incoming_longwave)
        ground_upward = last(longwave_buffers.boundary_upward_longwave)
        @test result.ground_downward_longwave ≈ downward_longwave + longwave_return_fraction * ground_upward
    end
end

@testset "Bare soil under longwave from its own temperature" begin
    emission = Microclimate.σ * (300.0u"K")^4
    emissivity = 0.9
    for (absorption, net_longwave) in ((EmissivityAbsorption(), 0.0u"W/m^2"), (FullAbsorption(), (1 - emissivity) * emission))
        reflected = Microclimate.soil_longwave_reflectance(ground(; emissivity, absorption)) * emission
        @test emission - reflected - emissivity * emission ≈ net_longwave atol=1e-9u"W/m^2"
    end
    # Only the unshaded ground exchanges longwave without a canopy; under one, all of it does
    shade = 0.3
    @test Microclimate.ground_longwave_terms(emission, shade, nothing) == (emission, 1 - shade)
    @test Microclimate.ground_longwave_terms(emission, shade, 2 * emission) == (2 * emission, 1.0)
end

@testset "BlackBodyLongwaveExchange" begin
    black_body_model = BlackBodyLongwaveExchange()
    black_body_buffers = Microclimate.allocate_longwave(black_body_model, plant_area_index, n_layers, canopy_projection_ratio)
    leaf_temperature = collect(range(280.0u"K", 293.0u"K"; length=n_layers))
    incoming_longwave = canopy_incoming_longwave(black_body_model; site, environment_instant=environment)
    ground_temperature = 285.0u"K"
    grey_ground = bare_ground_surface(; temperature=ground_temperature, emissivity=0.95)
    result = Microclimate.canopy_longwave!(black_body_buffers, black_body_model, leaf_emissivity;
        leaf_temperature, ground=grey_ground, incoming_longwave)

    # micropoint's longwavemodelCpp, written out: each layer passes on what it transmits and sends black-body leaf longwave, over the ground's own surface longwave
    transmission = exp.(-fill(plant_area_index / n_layers, n_layers))
    leaf_flux = Microclimate.σ .* leaf_temperature .^ 4
    downward = accumulate((flux, i) -> transmission[i] * flux + (1 - transmission[i]) * leaf_flux[i], 1:n_layers; init=incoming_longwave)
    upward = reverse(accumulate((flux, i) -> transmission[i] * flux + (1 - transmission[i]) * leaf_flux[i], n_layers:-1:1;
        init=grey_ground.emitted_longwave + grey_ground.longwave_reflectance * last(downward)))
    @test black_body_buffers.boundary_downward_longwave[2:end] ≈ downward
    @test black_body_buffers.boundary_upward_longwave[1:n_layers] ≈ upward
    @test result.ground_absorbed_longwave ≈ (1 - grey_ground.longwave_reflectance) * last(downward)

    # With black leaves and a black ground it is LayeredLongwaveExchange
    layered_model = LayeredLongwaveExchange()
    layered_buffers = Microclimate.allocate_longwave(layered_model, plant_area_index, n_layers, canopy_projection_ratio)
    black_ground = bare_ground_surface(; temperature=ground_temperature, emissivity=1.0)
    Microclimate.canopy_longwave!(black_body_buffers, black_body_model, 1.0; leaf_temperature, ground=black_ground, incoming_longwave)
    Microclimate.canopy_longwave!(layered_buffers, layered_model, 1.0; leaf_temperature, ground=black_ground, incoming_longwave)
    @test black_body_buffers.boundary_upward_longwave ≈ layered_buffers.boundary_upward_longwave
    @test black_body_buffers.boundary_downward_longwave ≈ layered_buffers.boundary_downward_longwave
    @test black_body_buffers.absorbed_longwave ≈ layered_buffers.absorbed_longwave
end
