using Microclimate
using Microclimate: allocate_canopy, canopy_energy_balance!
using Unitful
using Test

site = Microclimate.example_site()
canopy_height = 1.0u"m"
plant_area_index = 3.0
heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"
boundary_layer_model = MoninObukhov()

model = Microclimate.example_multilayer_canopy(;
    canopy_height, plant_area_index,
    convergence_model = PicardCanopyConvergence(;
        convergence = IterationToleranceConvergence(; tolerance = 0.01u"K", max_iterations_per_day = 15)),
)
buffers = Microclimate.allocate_canopy(model, heights, boundary_layer_model)

function make_environment_instant(; zenith_angle, reference_temperature = 293.0u"K")
    (;
        atmospheric_pressure = 101325.0u"Pa",
        reference_temperature,
        reference_wind_speed = 2.0u"m/s",
        reference_humidity = 0.5,
        zenith_angle,
        cloud_emissivity = 1.0,
        cloud_cover = 0.3,
        surface_emissivity = 0.95,
        shade = 0.0,
    )
end

@testset "midday: runs, converges, physically reasonable" begin
    environment_instant = make_environment_instant(; zenith_angle = 30.0u"°")
    inputs = Microclimate.CanopyEnergyBalanceInputs(model;
        site, environment_instant, zenith_angle = 30.0u"°",
        direct_horizontal_irradiance = 500.0u"W/m^2", diffuse_horizontal_irradiance = 100.0u"W/m^2",
        ground_reflectance = 0.15, ground_temperature = 290.0u"K", ground_emissivity = 0.95,
        ground_relative_humidity = 0.5, canopy_source_temperature = 293.0u"K",
    )
    result = Microclimate.canopy_energy_balance!(buffers, model, boundary_layer_model, inputs)

    @test result.iterations <= Microclimate.max_iterations(model.convergence_model.convergence)
    @test all(isfinite, ustrip.(u"K", buffers.leaf.leaf_temperature))
    @test all(t -> 250.0u"K" < t < 340.0u"K", buffers.leaf.leaf_temperature)
    @test all(isfinite, ustrip.(u"K", buffers.air_profile.air_temperature))
    @test result.landscape_absorbed_shortwave > 0.0u"W/m^2"
    @test result.canopy_potential_transpiration >= sum(buffers.leaf.evaporation_mass_flow) / 1.0u"m^2"

    # Visual check — run manually (not in CI)
    #= using Plots
    layer_heights = sort(heights[heights .<= canopy_height]; rev=true)
    plot(ustrip.(u"°C", buffers.air_profile.air_temperature), ustrip.(u"m", layer_heights);
        xlabel="Air temperature (°C)", ylabel="Height (m)", label="canopy air", marker=:circle) |> display
    plot(ustrip.(u"m/s", buffers.wind.wind_speed), ustrip.(u"m", layer_heights);
        xlabel="Wind speed (m/s)", ylabel="Height (m)", label="canopy wind", marker=:circle) |> display
    =#
end

@testset "update_canopy_energy_balance_inputs!: in-place refresh drives day -> night" begin
    inputs = Microclimate.CanopyEnergyBalanceInputs(model;
        site, environment_instant = make_environment_instant(; zenith_angle = 30.0u"°"), zenith_angle = 30.0u"°",
        direct_horizontal_irradiance = 500.0u"W/m^2", diffuse_horizontal_irradiance = 100.0u"W/m^2",
        ground_reflectance = 0.15, ground_temperature = 290.0u"K", ground_emissivity = 0.95,
        ground_relative_humidity = 0.5, canopy_source_temperature = 293.0u"K",
    )
    day_result = Microclimate.canopy_energy_balance!(buffers, model, boundary_layer_model, inputs)
    @test day_result.landscape_absorbed_shortwave > 0.0u"W/m^2"

    Microclimate.update_canopy_energy_balance_inputs!(inputs;
        environment_instant = make_environment_instant(; zenith_angle = 100.0u"°", reference_temperature = 283.0u"K"),
        zenith_angle = 100.0u"°",
        direct_horizontal_irradiance = 0.0u"W/m^2", diffuse_horizontal_irradiance = 0.0u"W/m^2",
        ground_reflectance = 0.15, ground_temperature = 280.0u"K", ground_emissivity = 0.95,
        ground_relative_humidity = 0.5, canopy_source_temperature = 283.0u"K",
    )
    night_result = Microclimate.canopy_energy_balance!(buffers, model, boundary_layer_model, inputs)
    @test night_result.landscape_absorbed_shortwave == 0.0u"W/m^2"
    @test all(isfinite, ustrip.(u"K", buffers.leaf.leaf_temperature))
    @test all(t -> 250.0u"K" < t < 310.0u"K", buffers.leaf.leaf_temperature)
end

@testset "leaf_water_potential (Campbell coupling input) warms canopy leaves" begin
    mean_leaf_temperature(lwp) = begin
        m = Microclimate.example_multilayer_canopy(;
            canopy_height, plant_area_index, convergence_model = PicardCanopyConvergence(; convergence = FixedIterationConvergence(20)),
        )
        b = Microclimate.allocate_canopy(m, heights, boundary_layer_model)
        inputs = Microclimate.CanopyEnergyBalanceInputs(m;
            site, environment_instant = make_environment_instant(; zenith_angle = 30.0u"°", reference_temperature = 298.0u"K"),
            zenith_angle = 30.0u"°",
            direct_horizontal_irradiance = 500.0u"W/m^2", diffuse_horizontal_irradiance = 100.0u"W/m^2",
            ground_reflectance = 0.15, ground_temperature = 295.0u"K", ground_emissivity = 0.95,
            ground_relative_humidity = 0.5, canopy_source_temperature = 298.0u"K", leaf_water_potential = lwp,
        )
        Microclimate.canopy_energy_balance!(b, m, boundary_layer_model, inputs)
        sum(ustrip.(u"K", b.leaf.leaf_temperature)) / length(b.leaf.leaf_temperature)
    end
    hydrated = mean_leaf_temperature(0.0u"J/kg")
    stressed = mean_leaf_temperature(-1.0e6u"J/kg")
    @test stressed > hydrated  # less evaporative cooling under water stress
end

@testset "canopy_energy_balance! is allocation-light" begin
    # Fixed iteration count for a reproducible byte budget (IterationToleranceConvergence's
    # own convergence-dependent iteration count would make this threshold nondeterministic).
    fixed_model = Microclimate.example_multilayer_canopy(;
        canopy_height, plant_area_index, convergence_model = PicardCanopyConvergence(; convergence = FixedIterationConvergence(3)),
    )
    fixed_buffers = Microclimate.allocate_canopy(fixed_model, heights, boundary_layer_model)
    inputs = Microclimate.CanopyEnergyBalanceInputs(fixed_model;
        site, environment_instant = make_environment_instant(; zenith_angle = 30.0u"°"), zenith_angle = 30.0u"°",
        direct_horizontal_irradiance = 500.0u"W/m^2", diffuse_horizontal_irradiance = 100.0u"W/m^2",
        ground_reflectance = 0.15, ground_temperature = 290.0u"K", ground_emissivity = 0.95,
        ground_relative_humidity = 0.5, canopy_source_temperature = 293.0u"K",
    )
    f() = Microclimate.canopy_energy_balance!(fixed_buffers, fixed_model, boundary_layer_model, inputs)
    f() # warm up
    # ~3000-3500 bytes measured for 3 iterations x 10 layers (confirmed unchanged
    # before/after adding the potential-transpiration calc); generous headroom
    # above that, not a tight bound.
    @test (@allocated f()) < 8_000
end

@testset "vector plant_area_index matches equivalent scalar and stays allocation-light" begin
    run_once(pai) = begin
        m = Microclimate.example_multilayer_canopy(; canopy_height, plant_area_index = pai,
            convergence_model = PicardCanopyConvergence(; convergence = FixedIterationConvergence(3)))
        b = Microclimate.allocate_canopy(m, heights, boundary_layer_model)
        inputs = Microclimate.CanopyEnergyBalanceInputs(m;
            site, environment_instant = make_environment_instant(; zenith_angle = 30.0u"°"), zenith_angle = 30.0u"°",
            direct_horizontal_irradiance = 500.0u"W/m^2", diffuse_horizontal_irradiance = 100.0u"W/m^2",
            ground_reflectance = 0.15, ground_temperature = 290.0u"K", ground_emissivity = 0.95,
            ground_relative_humidity = 0.5, canopy_source_temperature = 293.0u"K",
        )
        f() = Microclimate.canopy_energy_balance!(b, m, boundary_layer_model, inputs)
        f() # warm up
        (; leaf_temperature = b.leaf.leaf_temperature, bytes = @allocated f())
    end

    scalar = run_once(plant_area_index)
    uniform_vector = run_once(fill(plant_area_index / 10, 10))
    @test all(isapprox.(ustrip.(u"K", scalar.leaf_temperature), ustrip.(u"K", uniform_vector.leaf_temperature); atol=1e-9))
    @test uniform_vector.bytes < 8_000

    # Non-uniform (top-heavy) profile: runs finite and stays allocation-light too.
    top_heavy = [0.6, 0.5, 0.45, 0.4, 0.35, 0.3, 0.2, 0.15, 0.1, 0.05]
    result = run_once(top_heavy)
    @test all(isfinite, ustrip.(u"K", result.leaf_temperature))
    @test result.bytes < 8_000
end

@testset "fixed-iteration convergence strategy still runs to completion" begin
    fixed_model = Microclimate.example_multilayer_canopy(;
        canopy_height, plant_area_index, convergence_model = PicardCanopyConvergence(; convergence = FixedIterationConvergence(3)),
    )
    fixed_buffers = Microclimate.allocate_canopy(fixed_model, heights, boundary_layer_model)
    inputs = Microclimate.CanopyEnergyBalanceInputs(fixed_model;
        site, environment_instant = make_environment_instant(; zenith_angle = 30.0u"°"), zenith_angle = 30.0u"°",
        direct_horizontal_irradiance = 500.0u"W/m^2", diffuse_horizontal_irradiance = 100.0u"W/m^2",
        ground_reflectance = 0.15, ground_temperature = 290.0u"K", ground_emissivity = 0.95,
        ground_relative_humidity = 0.5, canopy_source_temperature = 293.0u"K",
    )
    result = Microclimate.canopy_energy_balance!(fixed_buffers, fixed_model, boundary_layer_model, inputs)
    @test result.iterations == 3
end
