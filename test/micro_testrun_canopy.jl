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
=#
