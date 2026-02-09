using Microclimate
using Test
using Statistics

@testset "Cold Air Pooling" begin
    # DEM with a central depression
    dem = Float32[
        10 10 10 10 10
        10  8  7  8 10
        10  7  5  7 10
        10  8  7  8 10
        10 10 10 10 10
    ]

    @testset "ColdAirFlow" begin
        method = ColdAirFlow(production=0.1, friction=0.1, timesteps=50)
        result = cold_air_pooling(method, dem; cellsize=(1.0, 1.0))

        @test size(result) == size(dem)

        # Cold air should accumulate in the depression (center)
        @test result[3, 3] > 0

        # High points (corners) should have no cold air - it drains away
        @test result[1, 1] ≈ 0.0 atol=1e-5
        @test result[1, 5] ≈ 0.0 atol=1e-5
        @test result[5, 1] ≈ 0.0 atol=1e-5
        @test result[5, 5] ≈ 0.0 atol=1e-5

        # Center should have maximum cold air
        @test result[3, 3] == maximum(result)

        # All values should be non-negative
        @test all(result .>= 0)
    end

    @testset "ColdAirFlow parameters" begin
        # More production = more cold air
        low_prod = cold_air_pooling(ColdAirFlow(production=0.05, timesteps=50), dem; cellsize=(1.0, 1.0))
        high_prod = cold_air_pooling(ColdAirFlow(production=0.2, timesteps=50), dem; cellsize=(1.0, 1.0))
        @test maximum(high_prod) > maximum(low_prod)

        # More timesteps = more accumulation
        few_steps = cold_air_pooling(ColdAirFlow(timesteps=20), dem; cellsize=(1.0, 1.0))
        many_steps = cold_air_pooling(ColdAirFlow(timesteps=100), dem; cellsize=(1.0, 1.0))
        @test maximum(many_steps) > maximum(few_steps)

        # Higher friction = slower drainage = more pooling on slopes
        low_fric = cold_air_pooling(ColdAirFlow(friction=0.05, timesteps=50), dem; cellsize=(1.0, 1.0))
        high_fric = cold_air_pooling(ColdAirFlow(friction=0.5, timesteps=50), dem; cellsize=(1.0, 1.0))
        # With high friction, intermediate cells retain more air
        @test low_fric[2, 2] <= high_fric[2, 2]
    end

    @testset "Flat terrain" begin
        # On flat terrain, cold air should accumulate uniformly
        flat_dem = fill(10.0f0, 5, 5)
        result = cold_air_pooling(ColdAirFlow(production=0.1, timesteps=10), flat_dem; cellsize=(1.0, 1.0))
        # All cells should have similar amounts (edges may differ due to boundary)
        interior = result[2:4, 2:4]
        @test std(interior) < 0.1 * mean(interior)
    end
end
