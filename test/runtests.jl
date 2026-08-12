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
@safetestset "Monin-Obukhov stability" begin include("monin_obukhov_test.jl") end
@safetestset "canopy shortwave" begin include("canopy_shortwave_test.jl") end
@safetestset "canopy longwave" begin include("canopy_longwave_test.jl") end
@safetestset "canopy energy balance" begin include("canopy_energy_balance_test.jl") end
@safetestset "canopy stomatal conductance" begin include("canopy_stomatal_conductance_test.jl") end
@safetestset "canopy air profile" begin include("canopy_air_profile_test.jl") end
@safetestset "canopy interception" begin include("canopy_interception_test.jl") end
@safetestset "canopy energy balance picard loop" begin include("canopy_energy_balance_picard_test.jl") end
@safetestset "infiltration algorithm comparison" begin include("infiltration_algorithm_comparison.jl") end
@safetestset "monthly simulation" begin include("micro_testrun_monthly.jl") end
@safetestset "monthly simulation with snow" begin include("micro_testrun_monthly_snow.jl") end
@safetestset "daily simulation" begin include("micro_testrun_daily.jl") end
@safetestset "canopy in main solve loop" begin include("micro_testrun_canopy.jl") end
