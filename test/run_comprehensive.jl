# julia --project=test test/run_comprehensive.jl
using Test

@testset "Comprehensive" begin
    include("surface_water_comprehensive.jl")
    include("cold_air_comprehensive.jl")
    include("component_interactions.jl")
end
