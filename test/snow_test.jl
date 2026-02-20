using Test
using Microclimate
using Unitful: m, mm, kg, K, °C, W, J, s, hr, d, Pa, @u_str

@testset "Snow Models" begin

    @testset "DegreeDaySnow" begin
        model = DegreeDaySnow()

        # Test precipitation partitioning
        @testset "Precipitation partitioning" begin
            # Cold: all snow
            result = snow_accumulation(model, 0.01u"m", -5.0u"°C")
            @test result.snowfall ≈ 0.01u"m"
            @test result.rainfall ≈ 0.0u"m"

            # Warm: all rain
            result = snow_accumulation(model, 0.01u"m", 5.0u"°C")
            @test result.snowfall ≈ 0.0u"m"
            @test result.rainfall ≈ 0.01u"m"

            # At threshold: all snow
            result = snow_accumulation(model, 0.01u"m", 1.0u"°C")
            @test result.snowfall ≈ 0.01u"m"
        end

        @testset "Melt calculation" begin
            # No melt below threshold
            state = SnowState(water_equivalent=0.1u"m")
            forcing = SnowForcing(air_temperature=-5.0u"°C", precipitation=0.0u"m")
            melt = snow_melt(model, state, forcing, 24.0u"hr")
            @test melt ≈ 0.0u"m"

            # Melt above threshold
            state = SnowState(water_equivalent=0.1u"m")
            forcing = SnowForcing(air_temperature=5.0u"°C", precipitation=0.0u"m")
            melt = snow_melt(model, state, forcing, 1.0u"d")
            expected = 3.0u"mm/K/d" * 5.0u"K" * 1.0u"d"  # temperature difference is in K
            @test melt ≈ expected

            # Can't melt more than available
            state = SnowState(water_equivalent=0.001u"m")
            forcing = SnowForcing(air_temperature=10.0u"°C", precipitation=0.0u"m")
            melt = snow_melt(model, state, forcing, 1.0u"d")
            @test melt ≈ 0.001u"m"

            # No snow, no melt
            state = SnowState(water_equivalent=0.0u"m")
            forcing = SnowForcing(air_temperature=10.0u"°C", precipitation=0.0u"m")
            melt = snow_melt(model, state, forcing, 1.0u"d")
            @test melt ≈ 0.0u"m"
        end
    end

    @testset "Snow17" begin
        model = Snow17()

        @testset "Precipitation partitioning" begin
            # Cold: all snow
            result = snow_accumulation(model, 0.01u"m", -5.0u"°C")
            @test result.snowfall ≈ 0.01u"m"
            @test result.rainfall ≈ 0.0u"m"

            # Warm: all rain
            result = snow_accumulation(model, 0.01u"m", 5.0u"°C")
            @test result.snowfall ≈ 0.0u"m"
            @test result.rainfall ≈ 0.01u"m"

            # Transition zone: mixed
            result = snow_accumulation(model, 0.01u"m", 1.0u"°C")
            @test 0.0u"m" < result.snowfall < 0.01u"m"
            @test 0.0u"m" < result.rainfall < 0.01u"m"
            @test result.snowfall + result.rainfall ≈ 0.01u"m"
        end

        @testset "Seasonal melt factor" begin
            # Winter solstice (day ~355) should have minimum melt factor
            mf_winter = seasonal_melt_factor(model, 355)
            @test mf_winter < (model.maximum_melt_factor + model.minimum_melt_factor) / 2

            # Summer solstice (day ~172) should have maximum melt factor
            mf_summer = seasonal_melt_factor(model, 172)
            @test mf_summer > (model.maximum_melt_factor + model.minimum_melt_factor) / 2
        end
    end

    @testset "UtahEnergyBalance" begin
        model = UtahEnergyBalance()

        @testset "Precipitation partitioning" begin
            # Cold: all snow
            result = snow_accumulation(model, 0.01u"m", -5.0u"°C")
            @test result.snowfall ≈ 0.01u"m"

            # Warm: all rain
            result = snow_accumulation(model, 0.01u"m", 10.0u"°C")
            @test result.rainfall ≈ 0.01u"m"
        end

        @testset "Albedo decay" begin
            formula = model.albedo_formula

            # Fresh snow (age=0) - DickinsonAlbedo averages visible (0.85) and infrared (0.65) bands
            α = snow_albedo(formula, 0.0u"d")
            @test α ≈ 0.75  # (0.85 + 0.65) / 2

            # Aged snow decays
            α = snow_albedo(formula, 1.0u"d")
            @test α < 0.75
            @test α > 0.5
        end

        @testset "Density evolution" begin
            formula = model.density_formula

            # Fresh snow - mixing keeps density near fresh
            state = SnowState(density=100.0u"kg/m^3", water_equivalent=0.01u"m")
            density = snow_density(formula, state, 0.01u"m", 1.0u"hr")
            @test density ≈ 100.0u"kg/m^3" atol=5.0u"kg/m^3"

            # Compaction over time
            state = SnowState(density=100.0u"kg/m^3", water_equivalent=0.1u"m")
            density = snow_density(formula, state, 0.0u"m", 24.0u"hr")
            @test density > 100.0u"kg/m^3"
            @test density < 500.0u"kg/m^3"
        end
    end

    @testset "KearneySnow" begin
        sturm = SturmSnowDensity(
            fresh = 270.0u"kg/m^3",
            maximum = 917.0u"kg/m^3",
        )
        model = KearneySnow(
            accumulation = ThresholdAccumulation(threshold=1.0u"°C"),
            density_formula = sturm,
        )

        @testset "Sturm density formula" begin
            # Fresh shallow snow - some compaction even at 10cm
            state = SnowState(depth=0.1u"m", age=0.0u"d")
            ρ = snow_density(sturm, state, 0.0u"m", 1.0u"hr")
            @test ρ > 270.0u"kg/m^3"
            @test ρ < 400.0u"kg/m^3"

            # Deep old snow approaches maximum
            state = SnowState(depth=3.0u"m", age=100.0u"d")
            ρ = snow_density(sturm, state, 0.0u"m", 1.0u"hr")
            @test ρ > 500.0u"kg/m^3"
            @test ρ < 917.0u"kg/m^3"
        end

        @testset "Thermal conductivity" begin
            djachkova = Djachkova()

            # Low density = low conductivity
            k_low = snow_thermal_conductivity(djachkova, 100.0u"kg/m^3")

            # High density = high conductivity
            k_high = snow_thermal_conductivity(djachkova, 500.0u"kg/m^3")

            @test k_high > k_low
            @test k_low > 0.0u"W/m/K"
            @test k_low < 1.0u"W/m/K"
        end

        @testset "Anderson albedo decay" begin
            anderson = AndersonAlbedo()

            # Fresh snow
            α = snow_albedo(anderson, 0.5u"d")
            @test α ≈ 0.85

            # Aged snow
            α = snow_albedo(anderson, 10.0u"d")
            @test α < 0.85
            @test α > 0.45
        end
    end

    @testset "AndersonDensityEvolution" begin
        formula = AndersonDensityEvolution()
        simple = SimpleMixingDensity()

        # Base state for comparisons
        base_state = SnowState(
            water_equivalent = 100.0u"mm",
            depth = 0.5u"m",
            density = 150.0u"kg/m^3",
            temperature = -5.0u"°C",
            liquid_water = 0.0u"mm",
        )

        @testset "Time-dependent compaction" begin
            # Anderson: density increases over time even without snowfall
            ρ_6h = snow_density(formula, base_state, 0.0u"mm", 6.0u"hr", -5.0u"°C")
            ρ_24h = snow_density(formula, base_state, 0.0u"mm", 24.0u"hr", -5.0u"°C")

            @test ρ_6h > base_state.density
            @test ρ_24h > ρ_6h

            # Simple: density unchanged without snowfall
            ρ_simple = snow_density(simple, base_state, 0.0u"mm", 24.0u"hr")
            @test ρ_simple ≈ base_state.density
        end

        @testset "Wet snow compacts faster" begin
            dry_state = base_state
            wet_state = SnowState(
                water_equivalent = 100.0u"mm",
                depth = 0.5u"m",
                density = 150.0u"kg/m^3",
                temperature = 0.0u"°C",
                liquid_water = 5.0u"mm",
            )

            ρ_dry = snow_density(formula, dry_state, 0.0u"mm", 6.0u"hr", -5.0u"°C")
            ρ_wet = snow_density(formula, wet_state, 0.0u"mm", 6.0u"hr", 0.0u"°C")

            @test ρ_wet > ρ_dry
        end

        @testset "Deeper snowpack compacts faster" begin
            shallow_state = SnowState(
                water_equivalent = 50.0u"mm",
                density = 150.0u"kg/m^3",
                temperature = -5.0u"°C",
            )
            deep_state = SnowState(
                water_equivalent = 500.0u"mm",
                density = 150.0u"kg/m^3",
                temperature = -5.0u"°C",
            )

            ρ_shallow = snow_density(formula, shallow_state, 0.0u"mm", 6.0u"hr", -5.0u"°C")
            ρ_deep = snow_density(formula, deep_state, 0.0u"mm", 6.0u"hr", -5.0u"°C")

            @test ρ_deep > ρ_shallow
        end

        @testset "Warmer snow compacts faster" begin
            cold_state = SnowState(
                water_equivalent = 100.0u"mm",
                density = 150.0u"kg/m^3",
                temperature = -15.0u"°C",
            )
            warm_state = SnowState(
                water_equivalent = 100.0u"mm",
                density = 150.0u"kg/m^3",
                temperature = -2.0u"°C",
            )

            ρ_cold = snow_density(formula, cold_state, 0.0u"mm", 6.0u"hr", -15.0u"°C")
            ρ_warm = snow_density(formula, warm_state, 0.0u"mm", 6.0u"hr", -2.0u"°C")

            @test ρ_warm > ρ_cold
        end

        @testset "Metamorphism slows above threshold density" begin
            low_density_state = SnowState(
                water_equivalent = 100.0u"mm",
                density = 100.0u"kg/m^3",  # below ρd = 150
                temperature = -5.0u"°C",
            )
            high_density_state = SnowState(
                water_equivalent = 100.0u"mm",
                density = 300.0u"kg/m^3",  # above ρd = 150
                temperature = -5.0u"°C",
            )

            ρ_low = snow_density(formula, low_density_state, 0.0u"mm", 6.0u"hr", -5.0u"°C")
            ρ_high = snow_density(formula, high_density_state, 0.0u"mm", 6.0u"hr", -5.0u"°C")

            # Relative increase should be larger for low density snow
            increase_low = (ρ_low - low_density_state.density) / low_density_state.density
            increase_high = (ρ_high - high_density_state.density) / high_density_state.density

            @test increase_low > increase_high
        end

        @testset "Maximum density cap" begin
            near_max_state = SnowState(
                water_equivalent = 100.0u"mm",
                density = 590.0u"kg/m^3",
                temperature = 0.0u"°C",
                liquid_water = 10.0u"mm",  # wet to maximize compaction
            )

            # Even with aggressive compaction, should not exceed 600 kg/m³
            ρ = snow_density(formula, near_max_state, 0.0u"mm", 24.0u"hr", 0.0u"°C")
            @test ρ <= 600.0u"kg/m^3"
        end

        @testset "Fresh snow density is temperature-dependent" begin
            empty_state = SnowState()

            ρ_cold = snow_density(formula, empty_state, 10.0u"mm", 1.0u"hr", -20.0u"°C")
            ρ_warm = snow_density(formula, empty_state, 10.0u"mm", 1.0u"hr", -2.0u"°C")

            # Cold air produces lighter fresh snow
            @test ρ_cold < ρ_warm
            @test ρ_cold ≈ 50.0u"kg/m^3"  # minimum at Ta ≤ -15°C
        end
    end

    @testset "Snow thermal properties" begin
        djachkova = Djachkova()

        # Low density snow
        k = snow_thermal_conductivity(djachkova, 100.0u"kg/m^3")
        @test k > 0.0u"W/m/K"
        @test k < 0.5u"W/m/K"

        # Ice
        k = snow_thermal_conductivity(djachkova, 917.0u"kg/m^3")
        @test k > 1.0u"W/m/K"

        # Specific heat at different temperatures
        c_cold = snow_specific_heat(300.0u"kg/m^3", -10.0u"°C")
        c_melt = snow_specific_heat(300.0u"kg/m^3", 0.0u"°C")

        # Near melting point includes latent heat
        @test c_melt > c_cold
    end

    @testset "SnowState" begin
        state = SnowState()

        @test state.water_equivalent ≈ 0.0u"m"
        @test state.depth ≈ 0.0u"m"
        @test state.density ≈ 100.0u"kg/m^3"
        @test state.albedo ≈ 0.85

        # Update with DegreeDaySnow
        model = DegreeDaySnow()
        forcing = SnowForcing(air_temperature=-5.0u"°C", precipitation=0.02u"m")
        state, water_output = update_snow_state(model, state, forcing, 1.0u"hr")

        @test state.water_equivalent ≈ 0.02u"m"
        @test state.depth > 0.0u"m"
        @test water_output ≈ 0.0u"m"  # no melt, no rain
    end

end

@testset "Model Comparisons" begin
    # All models should agree on basic precipitation partitioning
    @testset "Precipitation partitioning consistency" begin
        # Use consistent accumulation thresholds for comparison
        accum_threshold = ThresholdAccumulation(threshold=1.0u"°C")
        accum_linear = LinearTransitionAccumulation(snow_threshold=0.0u"°C", rain_threshold=2.0u"°C")
        models = [
            DegreeDaySnow(accumulation=accum_threshold),
            Snow17(accumulation=accum_linear),
            UtahEnergyBalance(),  # uses default 0-3°C transition
            KearneySnow(accumulation=accum_threshold),
        ]

        # Cold conditions: all models should produce snow
        for model in models
            result = snow_accumulation(model, 0.01u"m", -10.0u"°C")
            @test result.snowfall ≈ 0.01u"m"
            @test result.rainfall ≈ 0.0u"m"
        end

        # Warm conditions: all models should produce rain
        for model in models
            result = snow_accumulation(model, 0.01u"m", 10.0u"°C")
            @test result.snowfall ≈ 0.0u"m"
            @test result.rainfall ≈ 0.01u"m"
        end
    end

    @testset "Melt behavior ordering" begin
        swe = 0.1u"m"
        T_air = 5.0u"°C"
        timestep = 1.0u"d"

        # Degree-day melt
        dd_model = DegreeDaySnow(melt_factor=3.0u"mm/K/d")
        state = SnowState(water_equivalent=swe)
        forcing = SnowForcing(air_temperature=T_air, precipitation=0.0u"m")
        dd_melt = snow_melt(dd_model, state, forcing, timestep)

        # All models should produce positive melt above freezing
        @test dd_melt > 0.0u"m"

        # Melt should be bounded by available SWE
        @test dd_melt <= swe
    end

    @testset "Thermal conductivity range" begin
        # Djachkova formula should give reasonable conductivity at typical snow densities
        djachkova = Djachkova()

        for ρ in [100.0, 300.0, 500.0]u"kg/m^3"
            k = snow_thermal_conductivity(djachkova, ρ)

            # Should be positive and in reasonable range
            @test k > 0.0u"W/m/K"
            @test k < 3.0u"W/m/K"
        end
    end

    @testset "Albedo decay comparison" begin
        # Both ExponentialSnowAlbedo and AndersonAlbedo have albedo decay
        exponential = ExponentialSnowAlbedo()
        anderson = AndersonAlbedo()

        # Fresh snow: both should be high
        α_exp = snow_albedo(exponential, 0.0u"d")
        α_anderson = snow_albedo(anderson, 0.1u"d")

        @test α_exp > 0.8
        @test α_anderson > 0.8

        # Aged snow: both should decay
        α_exp_aged = snow_albedo(exponential, 10.0u"d")
        α_anderson_aged = snow_albedo(anderson, 10.0u"d")

        @test α_exp_aged < 0.85
        @test α_anderson_aged < 0.85
    end
end

@testset "Spatial Snow" begin

    @testset "Wind redistribution method" begin
        method = WindDrivenSnowRedistribution(
            transport_coefficient = 0.1,
            critical_wind_speed = 5.0,
            timesteps = 5,
        )

        # Create simple DEM
        dem = [100.0 110.0 120.0;
               100.0 110.0 120.0;
               100.0 110.0 120.0]

        # Uniform snow
        snow = fill(0.5, 3, 3)

        # Low wind: no redistribution
        result = snow_redistribution(method, dem, snow, 2.0, 270.0; cellsize=(10.0, 10.0))
        @test result ≈ snow

        # High wind: some redistribution
        result = snow_redistribution(method, dem, snow, 10.0, 270.0; cellsize=(10.0, 10.0))
        @test sum(result) ≈ sum(snow) atol=0.1  # mass conservation (approximate)
    end

    @testset "Gravitational redistribution" begin
        method = GravitationalSnowRedistribution(
            critical_slope = 30.0,
            timesteps = 5,
        )

        # Steep slope DEM
        dem = [100.0 100.0 100.0;
               50.0  50.0  50.0;
               0.0   0.0   0.0]

        snow = fill(0.5, 3, 3)

        result = snow_redistribution(method, dem, snow; cellsize=(10.0, 10.0))

        # Snow redistributes - lower cells should have more snow
        @test result[3, 2] >= result[1, 2]  # bottom row has more than top
        @test all(result .>= 0)  # no negative snow
    end

end
