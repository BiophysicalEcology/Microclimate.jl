using Microclimate, Test, Aqua, SafeTestsets

# Automated quality assurance checks
@testset "Aqua" begin
    Aqua.test_all(Microclimate;
        persistent_tasks=false,
    )
end

# Tests
@safetestset "monthly simulation" begin include("micro_testrun_monthly.jl") end
@safetestset "daily simulation" begin include("micro_testrun_daily.jl") end
