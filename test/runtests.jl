using Microclimate, Test, Aqua, SafeTestsets

# Automated quality assurance checks
@testset "Aqua" begin
    Aqua.test_all(Microclimate;
        persistent_tasks=false,
    )
end

# Tests
@safetestset "solar radiation" begin include("Solrad_tests.jl") end
@safetestset "monthly simulation" begin include("micro_testrun_monthly.jl") end
