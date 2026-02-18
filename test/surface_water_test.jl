using Microclimate
using Test
using Statistics

@testset "Surface Water Flow" begin
    # DEM with a central depression
    dem = Float32[
        10 10 10 10 10
        10  8  7  8 10
        10  7  5  7 10
        10  8  7  8 10
        10 10 10 10 10
    ]

    @testset "SurfaceWaterFlow basic" begin
        precip = 0.01  # 10mm uniform
        method = SurfaceWaterFlow(infiltration_rate=0.0001, friction=0.1, timesteps=50)
        result = surface_water_flow(method, dem, precip; cellsize=(1.0, 1.0))

        @test size(result) == size(dem)

        # Water should accumulate in depression
        @test result[3, 3] > 0

        # High points should drain
        @test result[1, 1] ≈ 0.0 atol=1e-5
        @test result[5, 5] ≈ 0.0 atol=1e-5

        # Center should have maximum
        @test result[3, 3] == maximum(result)

        # All non-negative
        @test all(result .>= 0)
    end

    @testset "Infiltration effect" begin
        # Test on flat terrain where flow doesn't complicate things
        flat = fill(10.0f0, 5, 5)
        precip = 0.01
        no_infil = surface_water_flow(
            SurfaceWaterFlow(infiltration_rate=0.0, timesteps=10),
            flat, precip; cellsize=(1.0, 1.0)
        )
        with_infil = surface_water_flow(
            SurfaceWaterFlow(infiltration_rate=0.005, timesteps=10),
            flat, precip; cellsize=(1.0, 1.0)
        )

        # With infiltration = less total surface water than without
        @test sum(with_infil) < sum(no_infil)
    end

    @testset "Spatial precipitation" begin
        # Precipitation only on left side
        precip_grid = Float32[
            0.02 0.02 0.0 0.0 0.0
            0.02 0.02 0.0 0.0 0.0
            0.02 0.02 0.0 0.0 0.0
            0.02 0.02 0.0 0.0 0.0
            0.02 0.02 0.0 0.0 0.0
        ]

        result = surface_water_flow(
            SurfaceWaterFlow(timesteps=100),
            dem, precip_grid; cellsize=(1.0, 1.0)
        )

        # Water should still flow to center depression
        @test result[3, 3] > 0
    end

    @testset "surface_water_event convenience" begin
        result = surface_water_event(dem, 0.01; duration=50, cellsize=(1.0, 1.0))
        @test size(result) == size(dem)
        @test result[3, 3] > 0
    end

    @testset "Flat terrain" begin
        flat_dem = fill(10.0f0, 5, 5)
        result = surface_water_flow(
            SurfaceWaterFlow(infiltration_rate=0.0, timesteps=10),
            flat_dem, 0.01; cellsize=(1.0, 1.0)
        )
        # On flat terrain with no infiltration, water stays where it falls
        interior = result[2:4, 2:4]
        @test std(interior) < 0.01  # Should be nearly uniform
    end
end
