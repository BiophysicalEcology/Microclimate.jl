using Test
using Microclimate
using Stencils: Wrap, Remove

mean(x) = sum(x) / length(x)

@testset "Surface Water Flow" begin
    @testset "Conservation" begin
        dem = Float32[10 10 10 10 10; 10 5 5 5 10; 10 5 2 5 10; 10 5 5 5 10; 10 10 10 10 10]
        method = Microclimate.SurfaceWaterFlow(infiltration_rate=0.0, friction=0.9, timesteps=1000)
        water = Microclimate.surface_water_flow(method, dem, 1.0; cellsize=(1.0, 1.0), boundary=Wrap())
        @test sum(water) ≈ 25.0 atol=0.01

        # Flat terrain - no flow
        dem_flat = fill(5.0f0, 10, 10)
        water_flat = Microclimate.surface_water_flow(method, dem_flat, 1.0; cellsize=(1.0, 1.0), boundary=Wrap())
        @test sum(water_flat) ≈ 100.0 atol=0.01
        @test all(w -> w ≈ 1.0, water_flat)

        # Infiltration removes water
        method_infil = Microclimate.SurfaceWaterFlow(infiltration_rate=0.1, friction=0.9, timesteps=10)
        water_infil = Microclimate.surface_water_flow(method_infil, dem_flat[1:5,1:5], 1.0; cellsize=(1.0, 1.0), boundary=Wrap())
        @test sum(water_infil) ≈ 0.0 atol=0.01
    end

    @testset "Lake filling" begin
        dem = Float32[10 10 10 10 10; 10 5 5 5 10; 10 5 2 5 10; 10 5 5 5 10; 10 10 10 10 10]
        method = Microclimate.SurfaceWaterFlow(infiltration_rate=0.0, friction=0.9, timesteps=1000)
        water = Microclimate.surface_water_flow(method, dem, 1.0; cellsize=(1.0, 1.0), boundary=Wrap())

        inner_surfaces = [dem[i,j] + water[i,j] for i in 2:4, j in 2:4]
        @test maximum(inner_surfaces) - minimum(inner_surfaces) < 0.01  # Flat lake
        @test inner_surfaces[2,2] ≈ 7.44 atol=0.1  # Expected equilibrium
    end

    @testset "Flow direction" begin
        # Valley collects water
        dem_valley = Float32[abs(i - 3) for i in 1:5, j in 1:5]
        method = Microclimate.SurfaceWaterFlow(infiltration_rate=0.0, friction=0.9, timesteps=500)
        water_valley = Microclimate.surface_water_flow(method, dem_valley, 1.0; cellsize=(1.0, 1.0), boundary=Wrap())
        @test water_valley[3, 3] > water_valley[1, 3]

        # Slope with Wrap - low side accumulates more
        dem_slope = Float32[i for i in 1:5, j in 1:5]
        water_slope = Microclimate.surface_water_flow(method, dem_slope, 1.0; cellsize=(1.0, 1.0), boundary=Wrap())
        @test mean(water_slope[1, :]) > mean(water_slope[5, :])
    end

    @testset "Boundary conditions" begin
        dem = Float32[i + j for i in 1:5, j in 1:5]
        method = Microclimate.SurfaceWaterFlow(infiltration_rate=0.0, friction=0.5, timesteps=200)

        water_remove = Microclimate.surface_water_flow(method, dem, 1.0; cellsize=(1.0, 1.0), boundary=Remove(0.0f0))
        @test sum(water_remove) < 25.0  # Loses water at edges

        method2 = Microclimate.SurfaceWaterFlow(infiltration_rate=0.0, friction=0.9, timesteps=200)
        water_wrap = Microclimate.surface_water_flow(method2, dem, 1.0; cellsize=(1.0, 1.0), boundary=Wrap())
        @test sum(water_wrap) ≈ 25.0 atol=0.01  # Conserves
    end

    @testset "Numerical stability" begin
        dem = rand(Float32, 20, 20) .* 100
        method = Microclimate.SurfaceWaterFlow(infiltration_rate=0.01, friction=0.5, timesteps=500)
        water = Microclimate.surface_water_flow(method, dem, 0.1; cellsize=(1.0, 1.0), boundary=Remove(0.0f0))
        @test all(w -> w >= 0, water)
        @test !any(isnan, water)
        @test !any(isinf, water)

        # Extreme precipitation
        dem_flat = fill(5.0f0, 5, 5)
        method2 = Microclimate.SurfaceWaterFlow(infiltration_rate=0.0, friction=0.9, timesteps=10)
        water_extreme = Microclimate.surface_water_flow(method2, dem_flat, 1000.0; cellsize=(1.0, 1.0), boundary=Wrap())
        @test sum(water_extreme) ≈ 25000.0 atol=1.0
        @test !any(isnan, water_extreme)
    end
end
