using Microclimate
using Unitful
using Test

@testset "example_* constructors build" begin
    @test example_site() isa Site
    @test example_monthly_weather() isa MonthlyMinMaxEnvironment
    @test example_daily_environment() isa DailyTimeseries
    @test example_hourly_environment() isa HourlyTimeseries
    @test example_soil_thermal_parameters() isa CampbelldeVriesSoilProperties
    @test example_soil_hydraulics() isa CampbellSoilHydraulics
    @test example_microclimate_problem() isa MicroProblem
end

@testset "example_microclimate_problem solves end-to-end" begin
    out = solve(example_microclimate_problem())
    @test size(out.soil_temperature, 1) == 12 * 24
    @test size(out.soil_temperature, 2) == 19
    @test !any(isnan, ustrip.(out.soil_temperature))
    @test all(out.soil_temperature .> 200u"K")
    @test all(out.soil_temperature .< 350u"K")
end
