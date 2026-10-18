using Microclimate
using Microclimate: allocate_interception, blend_stomatal_conductance, canopy_interception!,
    evaporation, wet_canopy_fraction
using Unitful
using Test

include("helpers.jl")

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

@testset "VerticalRainInterception: mass conservation, wind-independent" begin
    model = VerticalRainInterception(; leaf_water_storage_capacity = 0.1u"kg/m^2")
    buffers = Microclimate.allocate_interception(model, canopy_height, plant_area_index, heights, 10, boundary_layer_model)
    rainfall = 0.5u"kg/m^2"
    storage_before = sum(buffers.leaf_surface_water)

    result_calm = Microclimate.canopy_interception!(buffers, model, 1.0; rainfall, wind_speed = fill(0.0u"m/s", 10))
    storage_after = sum(buffers.leaf_surface_water)
    @test ustrip(u"kg/m^2", result_calm.ground_throughfall + (storage_after - storage_before)) ≈
          ustrip(u"kg/m^2", rainfall) atol=1e-10
    @test all(>=(0.0u"kg/m^2"), buffers.leaf_surface_water)

    fill!(buffers.leaf_surface_water, 0.0u"kg/m^2")
    result_windy = Microclimate.canopy_interception!(buffers, model, 1.0; rainfall, wind_speed = fill(5.0u"m/s", 10))
    @test result_calm.ground_throughfall == result_windy.ground_throughfall
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
    body = Microclimate.leaf_body(0.05u"m", 0.02u"m", 1.0)
    conv = Microclimate.leaf_convection(Microclimate.ElaborateLeafConvection(), body, 1.0u"m^2", 293.0u"K", 295.0u"K", 1.0u"m/s", 101325.0u"Pa")
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
    @test relative_difference < 2e-3
end

@testset "_cap_wet_surface_evaporation conserves mass and energy" begin
    air_T = 293.0u"K"
    flux = 1.0e-4u"kg/s"  # demands 0.36 kg/m^2 over the hour at wet_fraction=1
    sensible = 50.0u"W/m^2"

    @testset "ample storage: passes through unchanged" begin
        flux2, sensible2, actual = Microclimate._cap_wet_surface_evaporation(flux, sensible, 1.0, 1.0u"kg/m^2", 2.0u"kg/m^2", air_T)
        @test flux2 == flux
        @test sensible2 == sensible
        @test ustrip(u"kg/m^2", actual) ≈
            ustrip(u"kg/m^2", uconvert(u"kg", flux * Microclimate.CANOPY_TIMESTEP) / 1.0u"m^2")
    end

    @testset "scarce storage: withdrawal capped, shortfall moved to sensible" begin
        available = 0.05u"kg/m^2"
        flux2, sensible2, actual = Microclimate._cap_wet_surface_evaporation(flux, sensible, 1.0, available, 0.1u"kg/m^2", air_T)
        @test actual == available
        @test flux2 < flux
        latent_removed = uconvert(u"W", (flux - flux2) * Microclimate.enthalpy_of_vaporisation(air_T))
        sensible_added = (sensible2 - sensible) * 1.0u"m^2"
        @test ustrip(u"W", latent_removed) ≈ ustrip(u"W", sensible_added) rtol=1e-10
        @test ustrip(u"kg/m^2", actual) ≈
            ustrip(u"kg/m^2", uconvert(u"kg", flux2 * Microclimate.CANOPY_TIMESTEP) / 1.0u"m^2") atol=1e-12
    end

    @testset "no storage: wet-surface share fully rerouted, dry (transpiration) share untouched" begin
        wet_fraction = 0.5
        flux2, sensible2, actual = Microclimate._cap_wet_surface_evaporation(flux, sensible, wet_fraction, 0.0u"kg/m^2", 0.1u"kg/m^2", air_T)
        @test actual == 0.0u"kg/m^2"
        @test ustrip(u"kg/s", flux2) ≈ ustrip(u"kg/s", flux * (1.0 - wet_fraction)) rtol=1e-10
    end

    @testset "condensation overflow: storage capped at capacity, shortfall moved to sensible" begin
        condensing_flux = -1.0e-4u"kg/s"
        available = 0.09u"kg/m^2"
        capacity = 0.1u"kg/m^2"
        flux2, sensible2, actual = Microclimate._cap_wet_surface_evaporation(condensing_flux, sensible, 1.0, available, capacity, air_T)
        @test actual == available - capacity  # capped at remaining headroom, not unbounded
        @test flux2 > condensing_flux  # less condensation "happened" than the uncapped flux implied
        latent_delta = uconvert(u"W", (condensing_flux - flux2) * Microclimate.enthalpy_of_vaporisation(air_T))
        sensible_delta = (sensible2 - sensible) * 1.0u"m^2"
        @test ustrip(u"W", latent_delta) ≈ ustrip(u"W", sensible_delta) rtol=1e-10
    end

    @testset "condensation within headroom: passes through unchanged" begin
        condensing_flux = -1.0e-6u"kg/s"  # small dew, well within headroom
        flux2, sensible2, actual = Microclimate._cap_wet_surface_evaporation(condensing_flux, sensible, 1.0, 0.05u"kg/m^2", 0.1u"kg/m^2", air_T)
        @test flux2 == condensing_flux
        @test sensible2 == sensible
        @test actual < 0.0u"kg/m^2"  # storage gains water
    end
end

@testset "canopy_energy_balance!: wet-surface evaporation can't exceed stored water" begin
    run_with_capacity(capacity) = begin
        model = Microclimate.example_multilayer_canopy(;
            canopy_height = 1.0u"m", plant_area_index = 3.0,
            interception_model = LayeredRainInterception(; leaf_water_storage_capacity = capacity),
        )
        buffers = Microclimate.allocate_canopy(model, heights, boundary_layer_model)
        inputs = canopy_inputs(model; zenith_angle = 20.0u"°", direct_horizontal_irradiance = 600.0u"W/m^2",
            reference_temperature = 293.0u"K", ground_temperature = 295.0u"K", canopy_source_temperature = 295.0u"K",
            rainfall = 5.0u"kg/m^2")
        Microclimate.canopy_energy_balance!(buffers, model, boundary_layer_model, inputs; site = Microclimate.example_site())
        buffers
    end

    scarce = run_with_capacity(1.0e-4u"kg/m^2")
    abundant = run_with_capacity(10.0u"kg/m^2")

    @test all(>=(0.0u"kg/m^2"), scarce.interception.leaf_surface_water)
    @test all(isfinite, ustrip.(u"kg/s", scarce.leaf.evaporation_mass_flow))
    @test all(isfinite, ustrip.(u"W/m^2", scarce.leaf.sensible_heat_source))
    # scarce storage forces some latent demand into sensible heat instead
    @test sum(scarce.leaf.evaporation_mass_flow) < sum(abundant.leaf.evaporation_mass_flow)
    @test sum(scarce.leaf.sensible_heat_source) > sum(abundant.leaf.sensible_heat_source)
end
