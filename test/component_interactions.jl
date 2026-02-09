using Test
using Microclimate
using Stencils: Wrap, Remove

mean(x) = sum(x) / length(x)

@testset "Component Interactions" begin
    @testset "Water and cold air share terrain response" begin
        dem = Float32[10 10 10 10 10; 10 5 5 5 10; 10 5 2 5 10; 10 5 5 5 10; 10 10 10 10 10]

        water = Microclimate.surface_water_flow(
            Microclimate.SurfaceWaterFlow(infiltration_rate=0.0, friction=0.9, timesteps=500),
            dem, 1.0; cellsize=(1.0, 1.0), boundary=Wrap()
        )
        air = Microclimate.cold_air_pooling(
            Microclimate.ColdAirFlow(production=1.0, friction=0.1, timesteps=100),
            dem; cellsize=(1.0, 1.0)
        )

        # Both accumulate in depression
        @test water[3, 3] > water[1, 1]
        @test air[3, 3] > air[1, 1]
        @test (water[3, 3] > water[1, 1]) == (air[3, 3] > air[1, 1])

        # Valley collects both
        dem_valley = Float32[abs(i - 5) for i in 1:9, j in 1:5]
        water_v = Microclimate.surface_water_flow(
            Microclimate.SurfaceWaterFlow(infiltration_rate=0.0, friction=0.9, timesteps=500),
            dem_valley, 1.0; cellsize=(1.0, 1.0), boundary=Wrap()
        )
        air_v = Microclimate.cold_air_pooling(
            Microclimate.ColdAirFlow(production=1.0, friction=0.1, timesteps=100),
            dem_valley; cellsize=(1.0, 1.0)
        )
        @test mean(water_v[5, :]) > mean(water_v[1, :])
        @test mean(air_v[5, :]) > mean(air_v[1, :])
    end

    @testset "Infiltration effects" begin
        dem = fill(5.0f0, 5, 5)
        water_low = Microclimate.surface_water_flow(
            Microclimate.SurfaceWaterFlow(infiltration_rate=0.01, friction=0.9, timesteps=100),
            dem, 1.0; cellsize=(1.0, 1.0), boundary=Wrap()
        )
        water_high = Microclimate.surface_water_flow(
            Microclimate.SurfaceWaterFlow(infiltration_rate=0.1, friction=0.9, timesteps=100),
            dem, 1.0; cellsize=(1.0, 1.0), boundary=Wrap()
        )
        @test sum(water_high) < sum(water_low)

        # Ponded water persists in depressions despite infiltration
        dem_pit = Float32[10 10 10; 10 2 10; 10 10 10]
        water_pit = Microclimate.surface_water_flow(
            Microclimate.SurfaceWaterFlow(infiltration_rate=0.01, friction=0.9, timesteps=100),
            dem_pit, 1.0; cellsize=(1.0, 1.0), boundary=Wrap()
        )
        @test water_pit[2, 2] > water_pit[1, 1]
    end

    @testset "Temporal stability" begin
        dem = Float32[10 10 10; 10 2 10; 10 10 10]

        # Water reaches stable equilibrium
        water_500 = Microclimate.surface_water_flow(
            Microclimate.SurfaceWaterFlow(infiltration_rate=0.0, friction=0.9, timesteps=500),
            dem, 1.0; cellsize=(1.0, 1.0), boundary=Wrap()
        )
        water_1000 = Microclimate.surface_water_flow(
            Microclimate.SurfaceWaterFlow(infiltration_rate=0.0, friction=0.9, timesteps=1000),
            dem, 1.0; cellsize=(1.0, 1.0), boundary=Wrap()
        )
        @test water_500 ≈ water_1000 atol=0.1

        # Cold air keeps accumulating
        air_50 = Microclimate.cold_air_pooling(
            Microclimate.ColdAirFlow(production=1.0, friction=0.1, timesteps=50),
            dem; cellsize=(1.0, 1.0)
        )
        air_100 = Microclimate.cold_air_pooling(
            Microclimate.ColdAirFlow(production=1.0, friction=0.1, timesteps=100),
            dem; cellsize=(1.0, 1.0)
        )
        @test sum(air_100) > sum(air_50)
    end

    @testset "Extreme terrain stability" begin
        # Steep slope
        dem_steep = Float32[Float32(i) for i in 1:10, j in 1:10]
        water_steep = Microclimate.surface_water_flow(
            Microclimate.SurfaceWaterFlow(infiltration_rate=0.0, friction=0.5, timesteps=200),
            dem_steep, 1.0; cellsize=(1.0, 1.0), boundary=Wrap()
        )
        air_steep = Microclimate.cold_air_pooling(
            Microclimate.ColdAirFlow(production=1.0, friction=0.1, timesteps=100),
            dem_steep; cellsize=(1.0, 1.0)
        )
        @test !any(isnan, water_steep)
        @test !any(isnan, air_steep)

        # Cliff
        dem_cliff = Float32[100 100 100; 100 0 100; 100 100 100]
        water_cliff = Microclimate.surface_water_flow(
            Microclimate.SurfaceWaterFlow(infiltration_rate=0.0, friction=0.9, timesteps=500),
            dem_cliff, 1.0; cellsize=(1.0, 1.0), boundary=Wrap()
        )
        @test !any(isnan, water_cliff)
        @test sum(water_cliff) ≈ 9.0 atol=0.1
        @test water_cliff[2, 2] > water_cliff[1, 1]
    end

    @testset "Real scenarios" begin
        # Hillslope with closed basin at bottom
        dem_hill = Float32[10 9 8 7 6; 10 9 8 7 5; 10 9 8 7 4; 10 9 8 7 3; 10 9 8 7 2]
        water_hill = Microclimate.surface_water_flow(
            Microclimate.SurfaceWaterFlow(infiltration_rate=0.0, friction=0.5, timesteps=300),
            dem_hill, 1.0; cellsize=(1.0, 1.0), boundary=Wrap()
        )
        @test water_hill[5, 5] > water_hill[1, 1]  # Low corner > high corner

        # Frost hollow
        dem_bowl = Float32[sqrt((i-5)^2 + (j-5)^2) for i in 1:9, j in 1:9]
        air_bowl = Microclimate.cold_air_pooling(
            Microclimate.ColdAirFlow(production=1.0, friction=0.1, timesteps=200),
            dem_bowl; cellsize=(1.0, 1.0)
        )
        @test air_bowl[5, 5] > air_bowl[1, 1]
        @test air_bowl[5, 5] > mean([air_bowl[1,5], air_bowl[9,5], air_bowl[5,1], air_bowl[5,9]])
    end
end
