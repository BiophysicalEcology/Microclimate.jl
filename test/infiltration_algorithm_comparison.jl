using Microclimate
using Unitful
using Test
using Microclimate: allocate_soil_water_balance, infiltration_step!

# infiltration_step!/allocate_soil_water_balance are internal — white-box tested here since
# there's no external (e.g. NicheMapR) reference for MatricFluxPotentialAlgorithm.

function run_infiltration(algorithm; initial_moisture, evapotranspiration, nsteps)
    # example_soil_profile()'s defaults have bulk_density == mineral_density (zero porosity);
    # override with realistic values.
    profile = example_soil_profile(; bulk_density=1.3u"Mg/m^3", mineral_density=2.65u"Mg/m^3")
    num_layers = length(profile.hydraulics.campbell_b_parameter)
    depths = fill(0.15u"m", num_layers) .* (1:num_layers)
    input_soil_temperature = fill(293.0u"K", num_layers)
    model = example_soil_hydraulic_model(; infiltration_algorithm=algorithm)
    buffers = allocate_soil_water_balance(model, num_layers)
    soil_moisture = fill(initial_moisture, num_layers)
    local out
    for _ in 1:nsteps
        out = infiltration_step!(buffers, model;
            soil_profile=profile,
            depths,
            atmospheric_pressure=101325.0u"Pa",
            local_relative_humidity=0.5,
            leaf_area_index=1.0u"Mg/m^3",
            soil_moisture,
            evapotranspiration,
            input_soil_temperature,
            moisture_timestep=360.0u"s",
            moisture_tolerance=1e-6u"kg/m^2/s",
            moisture_max_iterations=200,
        )
        soil_moisture = out.soil_moisture
    end
    return soil_moisture
end

@testset "MatricPotentialAlgorithm vs MatricFluxPotentialAlgorithm" begin
    @testset "moderate moisture: results agree closely (textbook: 'results are identical')" begin
        psi_moisture = run_infiltration(MatricPotentialAlgorithm(); initial_moisture=0.25, evapotranspiration=5e-7u"kg/m^2/s", nsteps=20)
        phi_moisture = run_infiltration(MatricFluxPotentialAlgorithm(); initial_moisture=0.25, evapotranspiration=5e-7u"kg/m^2/s", nsteps=20)
        @test !any(isnan, psi_moisture)
        @test !any(isnan, phi_moisture)
        @test psi_moisture ≈ phi_moisture rtol=1e-2
    end

    @testset "wet moisture: results agree closely" begin
        psi_moisture = run_infiltration(MatricPotentialAlgorithm(); initial_moisture=0.40, evapotranspiration=1e-7u"kg/m^2/s", nsteps=20)
        phi_moisture = run_infiltration(MatricFluxPotentialAlgorithm(); initial_moisture=0.40, evapotranspiration=1e-7u"kg/m^2/s", nsteps=20)
        @test !any(isnan, psi_moisture)
        @test !any(isnan, phi_moisture)
        @test psi_moisture ≈ phi_moisture rtol=1e-2
    end

    @testset "dry moisture: both stay bounded and physically valid (may diverge — Program 8.2 is documented to be less accurate here)" begin
        psi_moisture = run_infiltration(MatricPotentialAlgorithm(); initial_moisture=0.05, evapotranspiration=1e-6u"kg/m^2/s", nsteps=40)
        phi_moisture = run_infiltration(MatricFluxPotentialAlgorithm(); initial_moisture=0.05, evapotranspiration=1e-6u"kg/m^2/s", nsteps=40)
        @test !any(isnan, psi_moisture)
        @test !any(isnan, phi_moisture)
        # saturation water content for bulk_density=1.3, mineral_density=2.65 is 1-1.3/2.65 ≈ 0.509
        @test all(0 .<= psi_moisture .<= 0.51)
        @test all(0 .<= phi_moisture .<= 0.51)
    end
end

@testset "MatricFluxPotentialAlgorithm end-to-end smoke test" begin
    soil_hydraulic_model = example_soil_hydraulic_model(; infiltration_algorithm=MatricFluxPotentialAlgorithm())
    config = MicroConfig(; soil_moisture_strategy=DynamicSoilMoisture())
    soil_profile = example_soil_profile(; bulk_density=1.3u"Mg/m^3", mineral_density=2.65u"Mg/m^3")
    out = solve(example_microclimate_problem(; soil_hydraulic_model, config, soil_profile))
    @test !any(isnan, ustrip.(out.soil_temperature))
    @test all(out.soil_temperature .> 200u"K")
    @test all(out.soil_temperature .< 350u"K")
    @test !any(isnan, out.soil_moisture)
    @test all(0 .<= out.soil_moisture .<= 1)
end
