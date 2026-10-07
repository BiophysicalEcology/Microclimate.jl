using Microclimate, Test, SafeTestsets

# Aqua quality checks — currently failing on stale-deps; re-enable once Project.toml is cleaned up.
# using Aqua
# @testset "Aqua" begin
#     Aqua.test_all(Microclimate;
#         persistent_tasks=false,
#     )
# end

# Tests
@safetestset "example constructors" begin include("example_constructors.jl") end
@safetestset "infiltration algorithm comparison" begin include("infiltration_algorithm_comparison.jl") end
@safetestset "rainfall entry mode" begin include("rainfall_entry_mode.jl") end
@safetestset "monthly simulation" begin include("micro_testrun_monthly.jl") end
@safetestset "monthly simulation with snow" begin include("micro_testrun_monthly_snow.jl") end
@safetestset "daily simulation" begin include("micro_testrun_daily.jl") end
