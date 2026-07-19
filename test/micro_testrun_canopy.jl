using Microclimate
using Unitful
using Test

heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"

@testset "MultilayerCanopy wired into the main solve loop: finite, physically sane" begin
    canopy_model = example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 2.0)
    problem = example_microclimate_problem(; heights, canopy_model)
    output = solve(problem)

    @test all(isfinite, ustrip.(u"K", output.soil_temperature))
    @test all(t -> 200.0u"K" < t < 340.0u"K", output.soil_temperature)
    @test all(isfinite, output.soil_moisture)

    @test size(output.canopy.leaf_temperature, 2) == 10
    @test all(isfinite, ustrip.(u"K", output.canopy.leaf_temperature))
    @test all(t -> 200.0u"K" < t < 340.0u"K", output.canopy.leaf_temperature)
    @test all(isfinite, ustrip.(u"K", output.canopy.air_temperature))
    @test all(isfinite, ustrip.(u"m/s", output.canopy.wind_speed))
    @test all(>=(0.0u"W/m^2"), output.canopy.ground_absorbed_shortwave)
    @test all(>(0), output.canopy.iterations)
end

@testset "NoCanopy leaves canopy output zero-columned" begin
    output = solve(example_microclimate_problem(; heights))
    @test size(output.canopy.leaf_temperature) == (size(output.soil_temperature, 1), 0)
    @test all(iszero, output.canopy.ground_absorbed_shortwave)
end

@testset "MultilayerCanopy changes soil temperature relative to NoCanopy" begin
    canopy_model = example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 3.0)
    canopy_output = solve(example_microclimate_problem(; heights, canopy_model))
    bare_output = solve(example_microclimate_problem(; heights))

    @test canopy_output.soil_temperature != bare_output.soil_temperature
end

@testset "MultilayerCanopy leaf temperature responds to soil moisture via Campbell coupling" begin
    canopy_model = example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 3.0)
    config = MicroConfig(; soil_moisture_strategy = DynamicSoilMoisture())
    n = length(Microclimate.DEFAULT_DEPTHS)

    wet_output = solve(example_microclimate_problem(; heights, canopy_model, config,
        initial_soil_moisture = fill(0.35, n)))
    dry_output = solve(example_microclimate_problem(; heights, canopy_model, config,
        initial_soil_moisture = fill(0.03, n)))

    @test all(isfinite, ustrip.(u"K", wet_output.canopy.leaf_temperature))
    @test all(isfinite, ustrip.(u"K", dry_output.canopy.leaf_temperature))
    # Different initial soil moisture propagates through Campbell's demand/supply
    # solve into a different leaf_water_potential and hence leaf temperature. The
    # *direction* of that difference over a multi-day run also depends on ground-
    # temperature feedback and how quickly soil moisture re-equilibrates (verified
    # separately, directly on canopy_energy_balance!, that more negative
    # leaf_water_potential warms the leaf — see canopy_energy_balance_picard_test.jl).
    @test wet_output.canopy.leaf_temperature != dry_output.canopy.leaf_temperature
end

# Visual check — run manually (not in CI)
#= using Plots
canopy_model = example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 3.0)
output = solve(example_microclimate_problem(; heights, canopy_model))
layer_heights = sort(heights[heights .<= canopy_model.canopy_height]; rev=true)

# Daytime hour profile: leaf/air temperature and wind speed by height
step = 24 * 5 + 13  # day 6, hour 13 (afternoon)
plot(ustrip.(u"°C", output.canopy.leaf_temperature[step, :]), ustrip.(u"m", layer_heights);
    xlabel="Leaf temperature (°C)", ylabel="Height (m)", label="leaf", marker=:circle) |> display
plot!(ustrip.(u"°C", output.canopy.air_temperature[step, :]), ustrip.(u"m", layer_heights);
    label="in-canopy air", marker=:square) |> display
plot(ustrip.(u"m/s", output.canopy.wind_speed[step, :]), ustrip.(u"m", layer_heights);
    xlabel="Wind speed (m/s)", ylabel="Height (m)", label="canopy wind", marker=:circle) |> display

# Top-of-canopy leaf temperature over time
plot(1:size(output.canopy.leaf_temperature, 1), ustrip.(u"°C", output.canopy.leaf_temperature[:, 1]);
    xlabel="Hour", ylabel="Top-layer leaf temperature (°C)", label=false) |> display

# Wet vs dry soil, full simulation: leaf temperature profile at a daytime hour.
# Effect is modest (~0.1-0.2°C) — Campbell's soil model re-equilibrates the top
# layer to ~0.11 within day 1 regardless of how extreme the initial condition
# is, so this isn't a strong demonstration of the coupling on its own.
config = MicroConfig(; soil_moisture_strategy = DynamicSoilMoisture())
n = length(Microclimate.DEFAULT_DEPTHS)
wet_output = solve(example_microclimate_problem(; heights, canopy_model, config, initial_soil_moisture = fill(0.35, n)))
dry_output = solve(example_microclimate_problem(; heights, canopy_model, config, initial_soil_moisture = fill(0.03, n)))
plot(ustrip.(u"°C", wet_output.canopy.leaf_temperature[step, :]), ustrip.(u"m", layer_heights);
    xlabel="Leaf temperature (°C)", ylabel="Height (m)", label="wet soil", marker=:circle) |> display
plot!(ustrip.(u"°C", dry_output.canopy.leaf_temperature[step, :]), ustrip.(u"m", layer_heights);
    label="dry soil", marker=:square) |> display

# Direct mechanism: canopy_energy_balance! sensitivity to leaf_water_potential
# alone, isolated from the soil model's own moisture dynamics above. This is
# the actual strength of the water-stress -> leaf-temperature coupling
# (order 10+ °C), not muted by how fast Campbell's soil re-equilibrates.
site = example_site()
boundary_layer_model = MoninObukhov()
environment_instant = (;
    atmospheric_pressure = 101325.0u"Pa", reference_temperature = 298.0u"K",
    reference_wind_speed = 2.0u"m/s", reference_humidity = 0.5, zenith_angle = 30.0u"°",
    cloud_emissivity = 1.0, cloud_cover = 0.3, surface_emissivity = 0.95, shade = 0.0,
)
lwp_values = [0.0, -0.5e6, -1.0e6, -1.5e6, -2.0e6, -3.0e6] .* u"J/kg"
lwp_leaf_temps = map(lwp_values) do lwp
    m = example_multilayer_canopy(; canopy_height = 1.0u"m", plant_area_index = 3.0,
        convergence = FixedSoilTemperatureIterations(20))
    b = Microclimate.allocate_canopy(m, heights, boundary_layer_model)
    inputs = Microclimate.CanopyEnergyBalanceInputs(m;
        site, environment_instant, zenith_angle = 30.0u"°",
        direct_horizontal_irradiance = 500.0u"W/m^2", diffuse_horizontal_irradiance = 100.0u"W/m^2",
        ground_reflectance = 0.15, ground_temperature = 295.0u"K", ground_emissivity = 0.95,
        canopy_source_temperature = 298.0u"K", leaf_water_potential = lwp,
    )
    Microclimate.canopy_energy_balance!(b, m, boundary_layer_model, inputs)
    sum(b.leaf.leaf_temperature) / length(b.leaf.leaf_temperature)
end
plot(lwp_values, lwp_leaf_temps;
    xlabel="Leaf water potential", ylabel="Mean leaf temperature", label=false, marker=:circle) |> display
=#
