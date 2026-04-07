using Microclimate
using Microclimate: cold_air_depth, temperature_disturbance, friction_coefficient,
                    effective_cold_air_depth, KLAM21State, KLAM21Coupling,
                    apply_cold_air_coupling!
using Test
using Unitful
using Unitful: m, s, K, J, kg, W, hr, Pa, ustrip

@testset "KLAM21 Cold Air Drainage" begin

    @testset "KLAM21 struct" begin
        # Default construction
        params = KLAM21()
        @test params.temperature_profile.reference_depth == 10.0u"m"
        @test params.temperature_profile.reference_temperature_disturbance == 3.0u"K"
        @test params.effective_fraction == 5/12
        @test params.gravitational_acceleration == 9.81u"m/s^2"
        @test params.heat_capacity == 1005.0u"J/kg/K"
        @test params.reference_density == 1.2u"kg/m^3"
        @test params.reference_temperature == 288.0u"K"
        @test params.conditions isa AlwaysDraining
        @test params.temperature_profile isa ParabolicProfile

        # Custom construction with custom temperature profile
        custom_profile = ParabolicProfile(reference_depth = 15.0u"m")
        custom = KLAM21(temperature_profile = custom_profile, mixing_length = 2.0u"m")
        @test custom.temperature_profile.reference_depth == 15.0u"m"
        @test custom.mixing_length == 2.0u"m"
    end

    @testset "cold_air_depth" begin
        params = KLAM21()

        # Zero heat deficit → zero depth
        @test cold_air_depth(params, 0.0u"J/m^2") == 0.0u"m"
        @test cold_air_depth(params, -10.0u"J/m^2") == 0.0u"m"

        # Positive heat deficit → positive depth
        E = 1000.0u"J/m^2"
        H = cold_air_depth(params, E)
        @test H > 0.0u"m"
        @test unit(H) == u"m"

        # More heat deficit → deeper layer
        E1 = 500.0u"J/m^2"
        E2 = 2000.0u"J/m^2"
        @test cold_air_depth(params, E2) > cold_air_depth(params, E1)

        # Volume reduction factor reduces depth
        H_full = cold_air_depth(params, E)
        H_reduced = cold_air_depth(params, E; volume_reduction_factor=0.5)
        @test H_reduced > H_full  # Less volume → more concentrated → deeper
    end

    @testset "temperature_disturbance" begin
        params = KLAM21()
        profile = params.temperature_profile

        # Zero depth → zero disturbance
        @test temperature_disturbance(params, 0.0u"m") == 0.0u"K"
        @test temperature_disturbance(params, -1.0u"m") == 0.0u"K"

        # At reference depth → reference disturbance
        ΔT = temperature_disturbance(params, profile.reference_depth)
        @test ΔT ≈ profile.reference_temperature_disturbance

        # Deeper layer → larger disturbance
        ΔT1 = temperature_disturbance(params, 5.0u"m")
        ΔT2 = temperature_disturbance(params, 20.0u"m")
        @test ΔT2 > ΔT1
    end

    @testset "effective_cold_air_depth" begin
        params = KLAM21()
        H = 12.0u"m"
        H_eff = effective_cold_air_depth(params, H)
        @test H_eff == params.effective_fraction * H
        @test H_eff == (5/12) * 12.0u"m"
        @test H_eff == 5.0u"m"
    end

    @testset "friction_coefficient" begin
        params = KLAM21()

        # Zero depth → zero friction
        @test friction_coefficient(params, 0.0u"m", 0.05u"m") == 0.0

        # Very shallow layer → maximum friction (1.0)
        @test friction_coefficient(params, 0.1u"m", 0.05u"m") == 1.0

        # Normal conditions
        c_star = friction_coefficient(params, 5.0u"m", 0.05u"m")
        @test 0.0 < c_star < 1.0

        # Higher roughness → more friction
        c_smooth = friction_coefficient(params, 5.0u"m", 0.01u"m")
        c_rough = friction_coefficient(params, 5.0u"m", 0.1u"m")
        @test c_rough > c_smooth
    end

    @testset "KLAM21State" begin
        nx, ny = 10, 8
        state = KLAM21State(nx, ny)

        @test size(state.heat_deficit) == (nx, ny)
        @test size(state.depth) == (nx, ny)
        @test size(state.temperature_disturbance) == (nx, ny)
        @test size(state.velocity_x) == (nx, ny)
        @test size(state.velocity_y) == (nx, ny)

        # Check units
        @test unit(state.heat_deficit[1,1]) == u"J/m^2"
        @test unit(state.depth[1,1]) == u"m"
        @test unit(state.temperature_disturbance[1,1]) == u"K"
        @test unit(state.velocity_x[1,1]) == u"m/s"

        # All initialized to zero
        @test all(ustrip.(u"J/m^2", state.heat_deficit) .== 0.0)
        @test all(ustrip.(u"m", state.depth) .== 0.0)
    end

    @testset "cold_air_step! with SpatialEnvironment" begin
        model = KLAM21()
        nx, ny = 10, 10
        state = KLAM21State(nx, ny)

        # Simple valley DEM (depression in center)
        elevation = [
            100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0
            100.0  90.0  85.0  80.0  75.0  75.0  80.0  85.0  90.0 100.0
            100.0  85.0  75.0  65.0  60.0  60.0  65.0  75.0  85.0 100.0
            100.0  80.0  65.0  55.0  50.0  50.0  55.0  65.0  80.0 100.0
            100.0  75.0  60.0  50.0  45.0  45.0  50.0  60.0  75.0 100.0
            100.0  75.0  60.0  50.0  45.0  45.0  50.0  60.0  75.0 100.0
            100.0  80.0  65.0  55.0  50.0  50.0  55.0  65.0  80.0 100.0
            100.0  85.0  75.0  65.0  60.0  60.0  65.0  75.0  85.0 100.0
            100.0  90.0  85.0  80.0  75.0  75.0  80.0  85.0  90.0 100.0
            100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0
        ]u"m"

        # Create terrain and environment
        terrain = SpatialMicroTerrain(elevation; cellsize = (30.0u"m", 30.0u"m"))
        weather = SpatialWeather(
            air_temperature = 288.0u"K",
            wind_speed = 2.0u"m/s",
            humidity = 0.5,
            pressure = 101325.0u"Pa",
            cloud_cover = 0.0,
        )
        land_surface = LandSurface(roughness_length = 0.05u"m")
        env = SpatialEnvironment(terrain; weather, land_surface)

        # Uniform heat loss rate (nighttime cooling)
        heat_loss_rate = fill(30.0u"W/m^2", nx, ny)

        # Run one timestep
        cold_air_step!(model, state, env;
                       timestep = 1.0u"hr",
                       heat_loss_rate = heat_loss_rate,
                       hours_since_sunset = 2.0)

        # Heat deficit should increase where there's heat loss
        @test any(ustrip.(u"J/m^2", state.heat_deficit) .> 0)

        # Cold air should start forming
        @test any(ustrip.(u"m", state.depth) .> 0)

        # Run more timesteps
        for _ in 1:10
            cold_air_step!(model, state, env;
                           timestep = 1.0u"hr",
                           heat_loss_rate = heat_loss_rate,
                           hours_since_sunset = 2.0)
        end

        # After multiple timesteps, cold air should pool in valley
        center_depth = ustrip(u"m", state.depth[5, 5])
        edge_depth = ustrip(u"m", state.depth[1, 1])
        @test center_depth > edge_depth

        # Velocities should be non-zero (drainage occurring)
        @test any(ustrip.(u"m/s", state.velocity_x) .!= 0) ||
              any(ustrip.(u"m/s", state.velocity_y) .!= 0)
    end

    @testset "SpatialMicroState with KLAM21" begin
        nx, ny, nz = 5, 5, 3
        state = SpatialMicroState(nx, ny, nz)

        # Check nested cold_air state
        @test state.cold_air isa KLAM21State
        @test size(state.cold_air.depth) == (nx, ny)
        @test unit(state.cold_air.heat_deficit[1,1]) == u"J/m^2"

        # Check main state fields have units
        @test unit(state.soil_temperature[1,1,1]) == u"K"
        @test unit(state.surface_water[1,1]) == u"m"
        @test unit(state.soil_moisture[1,1,1]) == u"m^3/m^3"
    end

    @testset "apply_cold_air_coupling!" begin
        nx, ny, nz = 5, 5, 3
        state = SpatialMicroState(nx, ny, nz; initial_soil_temperature = 290.0u"K")

        # Set some cold air depth and temperature disturbance
        state.cold_air.depth[3, 3] = 5.0u"m"
        state.cold_air.temperature_disturbance[3, 3] = 2.0u"K"

        initial_temp = state.soil_temperature[3, 3, 1]

        apply_cold_air_coupling!(state, KLAM21Coupling())

        # Surface temperature should be reduced by temperature disturbance
        @test state.soil_temperature[3, 3, 1] == initial_temp - 2.0u"K"

        # Cells without cold air should be unchanged
        @test state.soil_temperature[1, 1, 1] == 290.0u"K"
    end

    @testset "SpatialMicroTerrain with units" begin
        elevation = [
            100.0 100.0 100.0
             90.0  80.0  90.0
            100.0 100.0 100.0
        ]u"m"

        terrain = SpatialMicroTerrain(elevation; cellsize = (10.0u"m", 10.0u"m"))

        @test size(terrain) == (3, 3)
        @test terrain.cellsize == (10.0u"m", 10.0u"m")
        @test unit(terrain.dem[1,1]) == u"m"
        @test unit(terrain.filled_dem[1,1]) == u"m"
        @test unit(terrain.basin_depth[1,1]) == u"m"

        # Sky view factor should be dimensionless
        @test terrain.sky_view_factor[2,2] isa Real
        @test 0.0 <= terrain.sky_view_factor[2,2] <= 1.0
    end
end
