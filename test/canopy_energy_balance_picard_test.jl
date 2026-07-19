using Microclimate
using Unitful
using Test

site = Microclimate.example_site()
canopy_height = 1.0u"m"
plant_area_index = 3.0
heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"
boundary_layer_model = MoninObukhov()

model = Microclimate.example_multilayer_canopy(;
    canopy_height, plant_area_index,
    convergence = SoilTemperatureConvergenceTolerance(; tolerance = 0.01u"K", max_iterations_per_day = 15),
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
        canopy_source_temperature = 293.0u"K",
    )
    result = Microclimate.canopy_energy_balance!(buffers, model, boundary_layer_model, inputs)

    @test result.iterations <= Microclimate.max_iterations(model.convergence)
    @test all(isfinite, ustrip.(u"K", buffers.leaf.leaf_temperature))
    @test all(t -> 250.0u"K" < t < 340.0u"K", buffers.leaf.leaf_temperature)
    @test all(isfinite, ustrip.(u"K", buffers.air_profile.air_temperature))
    @test result.canopy_absorbed_shortwave > 0.0u"W/m^2"

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
        canopy_source_temperature = 293.0u"K",
    )
    day_result = Microclimate.canopy_energy_balance!(buffers, model, boundary_layer_model, inputs)
    @test day_result.canopy_absorbed_shortwave > 0.0u"W/m^2"

    Microclimate.update_canopy_energy_balance_inputs!(inputs;
        environment_instant = make_environment_instant(; zenith_angle = 100.0u"°", reference_temperature = 283.0u"K"),
        zenith_angle = 100.0u"°",
        direct_horizontal_irradiance = 0.0u"W/m^2", diffuse_horizontal_irradiance = 0.0u"W/m^2",
        ground_reflectance = 0.15, ground_temperature = 280.0u"K", ground_emissivity = 0.95,
        canopy_source_temperature = 283.0u"K",
    )
    night_result = Microclimate.canopy_energy_balance!(buffers, model, boundary_layer_model, inputs)
    @test night_result.canopy_absorbed_shortwave == 0.0u"W/m^2"
    @test all(isfinite, ustrip.(u"K", buffers.leaf.leaf_temperature))
    @test all(t -> 250.0u"K" < t < 310.0u"K", buffers.leaf.leaf_temperature)
end

@testset "fixed-iteration convergence strategy still runs to completion" begin
    fixed_model = Microclimate.example_multilayer_canopy(;
        canopy_height, plant_area_index, convergence = FixedSoilTemperatureIterations(3),
    )
    fixed_buffers = Microclimate.allocate_canopy(fixed_model, heights, boundary_layer_model)
    inputs = Microclimate.CanopyEnergyBalanceInputs(fixed_model;
        site, environment_instant = make_environment_instant(; zenith_angle = 30.0u"°"), zenith_angle = 30.0u"°",
        direct_horizontal_irradiance = 500.0u"W/m^2", diffuse_horizontal_irradiance = 100.0u"W/m^2",
        ground_reflectance = 0.15, ground_temperature = 290.0u"K", ground_emissivity = 0.95,
        canopy_source_temperature = 293.0u"K",
    )
    result = Microclimate.canopy_energy_balance!(fixed_buffers, fixed_model, boundary_layer_model, inputs)
    @test result.iterations == 3
end
