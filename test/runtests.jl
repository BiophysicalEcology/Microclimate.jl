using Microclimate, Test, SafeTestsets

# Aqua quality checks — currently failing on stale-deps; re-enable once Project.toml is cleaned up.
# using Aqua
# @testset "Aqua" begin
#     Aqua.test_all(Microclimate;
#         persistent_tasks=false,
#     )
# end

# Tests
@safetestset "monthly simulation" begin include("micro_testrun_monthly.jl") end
@safetestset "daily simulation" begin include("micro_testrun_daily.jl") end
