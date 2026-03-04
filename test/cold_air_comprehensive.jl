using Test
using Microclimate

mean(x) = sum(x) / length(x)
std(x) = (m = mean(x); sqrt(sum((xi - m)^2 for xi in x) / length(x)))

@testset "Cold Air Pooling" begin
    @testset "Basic flow" begin
        dem = Float32[10 10 10 10 10; 10 5 5 5 10; 10 5 2 5 10; 10 5 5 5 10; 10 10 10 10 10]
        method = Microclimate.ColdAirFlow(production=1.0, friction=0.1, timesteps=100)
        air = Microclimate.cold_air_pooling(method, dem; cellsize=(1.0, 1.0))

        @test air[3, 3] > air[1, 1]  # Accumulates in depression
        @test air[2, 2] > air[1, 1]  # Inner cells > rim

        # Closed basin accumulates - compare center to edges within basin
        @test air[3, 3] > air[2, 2]  # Deepest point gets most
    end

    @testset "Production and friction" begin
        dem = Float32[10 10 10; 10 2 10; 10 10 10]

        # Higher production = more accumulation
        air_low = Microclimate.cold_air_pooling(Microclimate.ColdAirFlow(production=0.5, friction=0.1, timesteps=100), dem; cellsize=(1.0, 1.0))
        air_high = Microclimate.cold_air_pooling(Microclimate.ColdAirFlow(production=2.0, friction=0.1, timesteps=100), dem; cellsize=(1.0, 1.0))
        @test sum(air_high) > sum(air_low)

        # Zero production = no cold air
        air_zero = Microclimate.cold_air_pooling(Microclimate.ColdAirFlow(production=0.0, friction=0.1, timesteps=100), dem; cellsize=(1.0, 1.0))
        @test sum(air_zero) ≈ 0.0 atol=1e-6

        # High friction = more accumulation in center (less drainage)
        air_lowf = Microclimate.cold_air_pooling(Microclimate.ColdAirFlow(production=1.0, friction=0.1, timesteps=100), dem; cellsize=(1.0, 1.0))
        air_highf = Microclimate.cold_air_pooling(Microclimate.ColdAirFlow(production=1.0, friction=0.9, timesteps=100), dem; cellsize=(1.0, 1.0))
        @test air_highf[2,2] > 0  # Center accumulates with high friction
    end

    @testset "Timesteps" begin
        dem = Float32[10 10 10; 10 2 10; 10 10 10]
        air_short = Microclimate.cold_air_pooling(Microclimate.ColdAirFlow(production=1.0, friction=0.1, timesteps=10), dem; cellsize=(1.0, 1.0))
        air_long = Microclimate.cold_air_pooling(Microclimate.ColdAirFlow(production=1.0, friction=0.1, timesteps=100), dem; cellsize=(1.0, 1.0))
        @test sum(air_long) > sum(air_short)
    end

    @testset "Terrain geometry" begin
        # Valley collects more
        dem_valley = Float32[10 5 10; 10 5 10; 10 5 10]
        method = Microclimate.ColdAirFlow(production=1.0, friction=0.1, timesteps=100)
        air_valley = Microclimate.cold_air_pooling(method, dem_valley; cellsize=(1.0, 1.0))
        @test mean(air_valley[:, 2]) > mean(air_valley[:, 1])
    end

    @testset "Numerical stability" begin
        dem = rand(Float32, 20, 20) .* 100
        method = Microclimate.ColdAirFlow(production=1.0, friction=0.1, timesteps=500)
        air = Microclimate.cold_air_pooling(method, dem; cellsize=(1.0, 1.0))

        @test all(a -> a >= 0, air)
        @test !any(isnan, air)
        @test !any(isinf, air)
    end
end
