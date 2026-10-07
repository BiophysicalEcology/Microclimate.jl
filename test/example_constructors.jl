using Microclimate
using Unitful
using Test

@testset "example_* constructors build" begin
    @test example_site() isa Site
    @test example_monthly_weather() isa MonthlyMinMaxEnvironment
    @test example_daily_environment() isa DailyTimeseries
    @test example_hourly_environment() isa HourlyTimeseries
    @test example_soil_properties_model() isa CampbelldeVriesSoilProperties
    @test example_soil_hydraulic_model() isa CampbellSoilHydraulics
    @test example_campbell_hydraulic_profile() isa CampbellHydraulicProfile
    @test example_soil_profile() isa SoilProfile
    @test example_soil_profile().hydraulics isa CampbellHydraulicProfile
    @test example_microclimate_problem() isa MicroProblem
end

@testset "consecutive_days dispatch" begin
    @test Microclimate.consecutive_days(example_monthly_weather()) == false
    daily_minmax = DailyMinMaxEnvironment(; forcings = minmax_forcings(;
        reference_temperature_min = [10, 11, 9]u"°C",
        reference_temperature_max = [20, 22, 19]u"°C",
        reference_wind_speed_min = [1, 1, 1]u"m/s",
        reference_wind_speed_max = [3, 3, 3]u"m/s",
        reference_humidity_min = [0.3, 0.3, 0.3],
        reference_humidity_max = [0.8, 0.8, 0.8],
        cloud_cover_min = [0.1, 0.1, 0.1],
        cloud_cover_max = [0.5, 0.5, 0.5],
    ))
    @test daily_minmax isa DailyMinMaxEnvironment
    @test Microclimate.consecutive_days(daily_minmax) == true
end

@testset "example_microclimate_problem solves end-to-end" begin
    out = solve(example_microclimate_problem())
    @test size(out.soil_temperature, 1) == 12 * 24
    @test size(out.soil_temperature, 2) == 19
    @test !any(isnan, ustrip.(out.soil_temperature))
    @test all(out.soil_temperature .> 200u"K")
    @test all(out.soil_temperature .< 350u"K")
end
