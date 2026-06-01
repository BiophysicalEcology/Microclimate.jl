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

@testset "example_microclimate_problem solves end-to-end" begin
    out = solve(example_microclimate_problem())
    @test size(out.soil_temperature, 1) == 12 * 24
    @test size(out.soil_temperature, 2) == 19
    @test !any(isnan, ustrip.(out.soil_temperature))
    @test all(out.soil_temperature .> 200u"K")
    @test all(out.soil_temperature .< 350u"K")
end
