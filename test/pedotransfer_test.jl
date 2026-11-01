using Microclimate, Unitful, Test

# Reference values from NicheMapR's pedotransfer.R formulas (as ported in MicroclimateMapper)
@testset "pedotransfer" begin
    loam = pedotransfer(CosbyMultivariate(), 18.7, 45.0, 36.3, 1.2u"Mg/m^3")
    @test loam.saturated_hydraulic_conductivity ≈ 0.0003933972750708022u"kg*s/m^3"
    @test pedotransfer(CosbyUnivariate(), 18.7, 45.0, 36.3, 1.2u"Mg/m^3").saturated_hydraulic_conductivity ≈
        0.0003376924472250286u"kg*s/m^3"
    @test pedotransfer(CampbellParticleSize(), 18.7, 45.0, 36.3, 1.2u"Mg/m^3").saturated_hydraulic_conductivity ≈
        0.00042320609658461354u"kg*s/m^3"
    @test pedotransfer(CosbyMultivariate(), 18.7, 45.0, 36.3, 1.9u"Mg/m^3") == loam
    dense = pedotransfer(CampbellParticleSize(), 46.2, 34.8, 19.0, 1.9u"Mg/m^3")
    @test dense.saturated_hydraulic_conductivity ≈ 1.2180520351835803e-7u"kg*s/m^3"

    profile = CampbellHydraulicProfile(CosbyMultivariate();
        clay_percent = [18.7, 46.2], silt_percent = [45.0, 34.8], sand_percent = [36.3, 19.0],
        bulk_density = [1.2, 1.9]u"Mg/m^3", root_density = [1.0, 0.5]u"m/m^3")
    @test profile.saturated_hydraulic_conductivity[1] == loam.saturated_hydraulic_conductivity
    @test profile.campbell_b_parameter[1] == loam.campbell_b_parameter
end
