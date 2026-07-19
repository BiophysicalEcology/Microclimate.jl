using Microclimate
using Unitful
using Test

canopy_height = 1.0u"m"
plant_area_index = 3.0
heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"
boundary_layer_model = MoninObukhov()

@testset "NoInterception: rain passes straight through" begin
    model = NoInterception()
    buffers = Microclimate.allocate_interception(model, canopy_height, plant_area_index, heights, 10, boundary_layer_model)
    @test buffers === nothing
    wind_speed = fill(1.0u"m/s", 10)
    result = Microclimate.canopy_interception!(buffers, model, 1.0; rainfall = 2.0u"kg/m^2", wind_speed)
    @test result.ground_throughfall == 2.0u"kg/m^2"
    @test Microclimate.wet_canopy_fraction(model, buffers, 1) == 0.0
end

@testset "LayeredRainInterception: mass conservation" begin
    model = LayeredRainInterception(; leaf_water_storage_capacity = 0.1u"kg/m^2")
    buffers = Microclimate.allocate_interception(model, canopy_height, plant_area_index, heights, 10, boundary_layer_model)
    wind_speed = fill(1.0u"m/s", 10)
    storage_before = sum(buffers.leaf_surface_water)

    rainfall = 0.5u"kg/m^2"
    result = Microclimate.canopy_interception!(buffers, model, 1.0; rainfall, wind_speed)
    storage_after = sum(buffers.leaf_surface_water)

    @test ustrip(u"kg/m^2", result.ground_throughfall + (storage_after - storage_before)) ≈
          ustrip(u"kg/m^2", rainfall) atol=1e-10
    @test result.ground_throughfall < rainfall  # some was intercepted
    @test all(>=(0.0u"kg/m^2"), buffers.leaf_surface_water)
end

@testset "LayeredRainInterception: heavy rain saturates storage and drips through" begin
    model = LayeredRainInterception(; leaf_water_storage_capacity = 0.01u"kg/m^2")
    buffers = Microclimate.allocate_interception(model, canopy_height, plant_area_index, heights, 10, boundary_layer_model)
    wind_speed = fill(1.0u"m/s", 10)

    result = Microclimate.canopy_interception!(buffers, model, 1.0; rainfall = 50.0u"kg/m^2", wind_speed)

    layer_capacities = model.leaf_water_storage_capacity .* buffers.layer_plant_area_index
    @test all(buffers.leaf_surface_water .<= layer_capacities .+ 1e-12u"kg/m^2")
    @test result.ground_throughfall > 40.0u"kg/m^2"  # most heavy rain still reaches the ground
end

@testset "LayeredRainInterception: no rain leaves storage untouched" begin
    model = LayeredRainInterception()
    buffers = Microclimate.allocate_interception(model, canopy_height, plant_area_index, heights, 10, boundary_layer_model)
    buffers.leaf_surface_water .= 0.02u"kg/m^2"
    wind_speed = fill(1.0u"m/s", 10)

    result = Microclimate.canopy_interception!(buffers, model, 1.0; rainfall = 0.0u"kg/m^2", wind_speed)

    @test result.ground_throughfall == 0.0u"kg/m^2"
    @test all(==(0.02u"kg/m^2"), buffers.leaf_surface_water)
end

@testset "wet_canopy_fraction tracks storage saturation" begin
    model = LayeredRainInterception(; leaf_water_storage_capacity = 0.1u"kg/m^2")
    buffers = Microclimate.allocate_interception(model, canopy_height, plant_area_index, heights, 10, boundary_layer_model)
    @test Microclimate.wet_canopy_fraction(model, buffers, 1) == 0.0
    layer_capacity = model.leaf_water_storage_capacity * buffers.layer_plant_area_index[1]
    buffers.leaf_surface_water[1] = layer_capacity / 2
    @test Microclimate.wet_canopy_fraction(model, buffers, 1) ≈ 0.5 atol=1e-10
    buffers.leaf_surface_water[1] = layer_capacity * 10  # clamps at 1 even if numerically over
    @test Microclimate.wet_canopy_fraction(model, buffers, 1) == 1.0
end

@testset "blend_stomatal_conductance" begin
    dry = LeafEvaporationParameters(;
        abaxial_vapour_conductance = 0.3u"mol/m^2/s", adaxial_vapour_conductance = 0.0u"mol/m^2/s",
        cuticular_conductance = 0.01u"mol/m^2/s",
    )
    dry_result = Microclimate.blend_stomatal_conductance(dry, 0.0)
    @test dry_result.abaxial_vapour_conductance == dry.abaxial_vapour_conductance
    @test dry_result.cuticular_conductance == dry.cuticular_conductance

    wet_result = Microclimate.blend_stomatal_conductance(dry, 1.0)
    @test wet_result.abaxial_vapour_conductance == Microclimate.WET_SURFACE_CONDUCTANCE
    @test wet_result.adaxial_vapour_conductance == Microclimate.WET_SURFACE_CONDUCTANCE
    @test wet_result.cuticular_conductance == dry.cuticular_conductance

    half_result = Microclimate.blend_stomatal_conductance(dry, 0.5)
    @test dry.abaxial_vapour_conductance < half_result.abaxial_vapour_conductance < Microclimate.WET_SURFACE_CONDUCTANCE
end

@testset "WET_SURFACE_CONDUCTANCE saturates leaf evaporation (boundary-layer-limited)" begin
    body = Microclimate.leaf_body(0.05u"m", 0.02u"m")
    conv = Microclimate.leaf_convection(body, 1.0u"m^2", 293.0u"K", 295.0u"K", 1.0u"m/s", 101325.0u"Pa")
    atmos = Microclimate.AtmosphericConditions(0.5, 1.0u"m/s", 101325.0u"Pa")

    function wet_evap(conductance)
        wet = LeafEvaporationParameters(;
            abaxial_vapour_conductance = conductance, adaxial_vapour_conductance = conductance,
            cuticular_conductance = 0.01u"mol/m^2/s",
        )
        Microclimate.HeatExchange.evaporation(wet, conv.mass_transfer_coefficient, atmos, 1.0u"m^2", 295.0u"K", 293.0u"K")
    end

    reference = wet_evap(Microclimate.WET_SURFACE_CONDUCTANCE)
    much_higher = wet_evap(10.0 * Microclimate.WET_SURFACE_CONDUCTANCE)
    relative_difference = abs(ustrip(u"W", much_higher.evaporation_heat_flow - reference.evaporation_heat_flow)) /
        abs(ustrip(u"W", reference.evaporation_heat_flow))
    @test relative_difference < 1e-3
end
