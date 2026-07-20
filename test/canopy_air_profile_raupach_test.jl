using Microclimate
using Unitful
using Test

canopy_height = 1.0u"m"
plant_area_index = 3.0
heights = vcat(0.1:0.1:1.0, [2.0]) .* u"m"
boundary_layer_model = MoninObukhov()
model = RaupachLTheoryAirProfile()
buffers = Microclimate.allocate_air_profile(model, canopy_height, plant_area_index, heights, 10, boundary_layer_model)
n_layers = length(buffers.air_temperature)

displacement_height = 0.65u"m"
friction_velocity = 0.3u"m/s"
atmospheric_pressure = 101325.0u"Pa"
obukhov_length = Inf * u"m"  # neutral: Φ_h -> 1
zero_evaporation = zeros(typeof(0.0u"g/s"), n_layers)
zero_source = zeros(typeof(0.0u"W/m^2"), n_layers)

@testset "no forcing: uniform boundary temperatures give a uniform profile" begin
    T = 290.0u"K"
    fill!(buffers.air_temperature, T)  # seed the prev-state read at function entry, matching energy_balance!'s own pre-Picard seeding
    result = Microclimate.canopy_air_profile!(buffers, model, boundary_layer_model;
        canopy_height, displacement_height, friction_velocity,
        canopy_top_air_temperature = T, canopy_top_relative_humidity = 0.5, ground_temperature = T,
        ground_relative_humidity = 0.5, sensible_heat_source = zero_source,
        evaporation_mass_flow = zero_evaporation, obukhov_length, atmospheric_pressure,
    )
    @test all(t -> isapprox(ustrip(u"K", t), ustrip(u"K", T); atol=1e-6), result.air_temperature)
end

@testset "a mid-canopy heat source produces a local warm bump" begin
    T = 290.0u"K"
    fill!(buffers.air_temperature, T)
    source = zeros(typeof(0.0u"W/m^2"), n_layers)
    source[n_layers ÷ 2] = 50.0u"W/m^2"
    result = Microclimate.canopy_air_profile!(buffers, model, boundary_layer_model;
        canopy_height, displacement_height, friction_velocity,
        canopy_top_air_temperature = T, canopy_top_relative_humidity = 0.5,
        ground_temperature = T, ground_relative_humidity = 0.5,
        sensible_heat_source = source, evaporation_mass_flow = zero_evaporation, obukhov_length, atmospheric_pressure,
    )
    @test all(isfinite, ustrip.(u"K", result.air_temperature))
    peak = argmax(result.air_temperature)
    @test result.air_temperature[peak] > T
    @test abs(peak - n_layers ÷ 2) <= 1
    # a genuine local bump: every other layer is cooler than the peak. Not
    # asserting strict monotonic decay either side, unlike K-theory's
    # symmetric diffusion — Raupach's σ_w(z) profile is asymmetric (increasing
    # toward canopy top), so the near-field kernel doesn't guarantee unimodal
    # decay step-by-step immediately around the source.
    @test all(t -> t < result.air_temperature[peak], result.air_temperature[1:end .!= peak])
end

@testset "no forcing: uniform boundary RH gives a uniform humidity profile" begin
    T = 290.0u"K"
    fill!(buffers.air_temperature, T)
    fill!(buffers.vapour_density, 0.0u"kg/m^3")  # force the cold-start reseed from this call's own boundary condition
    result = Microclimate.canopy_air_profile!(buffers, model, boundary_layer_model;
        canopy_height, displacement_height, friction_velocity,
        canopy_top_air_temperature = T, canopy_top_relative_humidity = 0.5, ground_temperature = T,
        ground_relative_humidity = 0.5, sensible_heat_source = zero_source,
        evaporation_mass_flow = zero_evaporation, obukhov_length, atmospheric_pressure,
    )
    @test all(rh -> isapprox(rh, 0.5; atol=1e-6), result.relative_humidity)
end

@testset "a mid-canopy evaporative source produces a local humidity bump" begin
    T = 290.0u"K"
    fill!(buffers.air_temperature, T)
    fill!(buffers.vapour_density, 0.0u"kg/m^3")
    source = zeros(typeof(0.0u"g/s"), n_layers)
    source[n_layers ÷ 2] = 1.0e-3u"g/s"
    result = Microclimate.canopy_air_profile!(buffers, model, boundary_layer_model;
        canopy_height, displacement_height, friction_velocity,
        canopy_top_air_temperature = T, canopy_top_relative_humidity = 0.5,
        ground_temperature = T, ground_relative_humidity = 0.5,
        sensible_heat_source = zero_source, evaporation_mass_flow = source, obukhov_length, atmospheric_pressure,
    )
    @test all(isfinite, result.relative_humidity)
    peak = argmax(result.relative_humidity)
    @test result.relative_humidity[peak] > 0.5
    @test abs(peak - n_layers ÷ 2) <= 1
    @test all(rh -> rh < result.relative_humidity[peak], result.relative_humidity[1:end .!= peak])
end

@testset "vapor density <-> relative humidity round-trip (FluidProperties primitives)" begin
    # Dedicated check on the conversion this model's humidity readout leans
    # on throughout (Stage 5's judgment call: diffuse vapor density, not
    # vapor pressure — this confirms the round-trip back to RH is exact).
    for T in (280.0u"K", 290.0u"K", 300.0u"K", 310.0u"K"), rh in (0.1, 0.3, 0.5, 0.7, 0.95)
        ρv = Microclimate.wet_air_properties(T, rh, atmospheric_pressure).vapour_density
        ρv_sat = Microclimate.wet_air_properties(T, 1.0, atmospheric_pressure).vapour_density
        rh_recovered = ustrip(u"kg/m^3", ρv) / ustrip(u"kg/m^3", ρv_sat)
        @test isapprox(rh_recovered, rh; atol=1e-6)
    end
end

@testset "canopy_air_profile! allocation stays bounded" begin
    T = 290.0u"K"
    fill!(buffers.air_temperature, T)
    source = fill(10.0u"W/m^2", n_layers)
    evap_source = fill(1.0e-4u"g/s", n_layers)
    f() = Microclimate.canopy_air_profile!(buffers, model, boundary_layer_model;
        canopy_height, displacement_height, friction_velocity,
        canopy_top_air_temperature = T, canopy_top_relative_humidity = 0.5,
        ground_temperature = 288.0u"K", ground_relative_humidity = 0.5,
        sensible_heat_source = source, evaporation_mass_flow = evap_source, obukhov_length, atmospheric_pressure,
    )
    f() # warm up
    @test (@allocated f()) < 5_000
end
