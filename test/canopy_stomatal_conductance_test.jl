using Microclimate
using Microclimate: allocate_canopy, stomatal_conductance
using Unitful
using Test

day = 30.0u"°"
night = 100.0u"°"

@testset "PrescribedStomatalConductance: day/night gating" begin
    model = PrescribedStomatalConductance()
    day_conductance = stomatal_conductance(model, day, 0.0u"J/kg")
    night_conductance = stomatal_conductance(model, night, 0.0u"J/kg")

    @test day_conductance.abaxial_vapour_conductance == model.conductance.abaxial_vapour_conductance
    @test night_conductance.abaxial_vapour_conductance == 0.0u"mol/m^2/s"
    @test night_conductance.cuticular_conductance == model.conductance.cuticular_conductance
end

@testset "MoistureResponsiveStomatalConductance: closes as water potential drops" begin
    model = MoistureResponsiveStomatalConductance()
    wet = stomatal_conductance(model, day, 0.0u"J/kg")
    dry = stomatal_conductance(model, day, -1.0e6u"J/kg")
    night_conductance = stomatal_conductance(model, night, 0.0u"J/kg")

    @test wet.abaxial_vapour_conductance ≈ model.conductance.abaxial_vapour_conductance rtol=1e-6
    @test 0.0u"mol/m^2/s" < dry.abaxial_vapour_conductance < wet.abaxial_vapour_conductance
    @test night_conductance.abaxial_vapour_conductance == 0.0u"mol/m^2/s"
    @test night_conductance.cuticular_conductance == model.conductance.cuticular_conductance
end

@testset "MultilayerCanopy composes sub-models" begin
    model = MultilayerCanopy(;
        canopy_height = 1.0u"m", plant_area_index = 3.0,
        stomatal_model = MoistureResponsiveStomatalConductance(),
    )
    @test model.shortwave_model isa TwoStreamRadiation
    @test model.wind_model isa CanopyWindAttenuation
    @test model.stomatal_model isa MoistureResponsiveStomatalConductance
    @test model.leaf_temperature_solver isa LinearizedLeafTemperature

    heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"
    boundary_layer_model = MoninObukhov()
    buffers = Microclimate.allocate_canopy(model, heights, boundary_layer_model)
    # shortwave, wind, and leaf buffers all present after composition, nested one field per sub-model
    @test haskey(buffers.shortwave, :absorbed_shortwave)
    @test haskey(buffers.wind, :wind_speed)
    @test haskey(buffers.leaf, :leaf_body)
    @test length(buffers.leaf.leaf_temperature) == length(buffers.shortwave.absorbed_shortwave)
end
