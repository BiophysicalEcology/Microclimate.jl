function allocate_soil_energy_balance(num_nodes::Int)
    layer_depths = fill(0.0u"cm", num_nodes + 1)
    heat_capacity = fill(1.0u"J/K/m^2", num_nodes)
    thermal_conductance = fill(1.0u"W/K/m^2", num_nodes)
    return (; layer_depths, heat_capacity, thermal_conductance)
end

# This is a 3-parameters OrdinaryDiffEq function
function soil_energy_balance(
    temperature_state::U,  # state
    p::SoilEnergyInputs,   # "parameters"
    t::Quantity,           # timestep
) where U <: SVector{N} where N
    # extract parameters
    (; soil_thermal_model, forcing, buffers, heights, depths, nodes, environment_instant, solar_terrain, micro_terrain, soil_wetness, runmoist) = p
    (; layer_depths, heat_capacity, thermal_conductance) = buffers.soil_energy_balance
    (; soil_moisture, shade) = environment_instant
    # Get environmental data at time t
    (; atmospheric_pressure, air_temperature, wind_speed, zenith_angle, solar_radiation, cloud_cover, relative_humidity, slope_zenith_angle) = interpolate_forcings(forcing, t)
    (; roughness_height, karman_constant, dyer_constant) = micro_terrain
    (; albedo, slope) = solar_terrain

    reference_height = last(heights)
    absorptivity = 1.0 - albedo

    # check for unstable conditions of ground surface temperature
    clamped_temperature = map(t -> clamp(t, (-81.0+273.15)u"K", (85.0+273.15)u"K"), temperature_state)::U

    # get soil properties
    (; bulk_thermal_conductivity, bulk_heat_capacity, bulk_density) = soil_properties!(buffers.soil_properties, soil_thermal_model; soil_temperature=clamped_temperature, soil_moisture, atmospheric_pressure)

    temperature_vector = SVector(MVector(clamped_temperature))

    # set boundary condition of deep soil temperature
    layer_depths[1:N] = depths
    # Compute soil layer properties
    for i in 1:N
        volumetric_heat_capacity = bulk_density[i] * bulk_heat_capacity[i]
        if i == 1
            heat_capacity[i] = volumetric_heat_capacity * layer_depths[2] / 2.0
            thermal_conductance[i] = bulk_thermal_conductivity[1] / layer_depths[2]
        else
            heat_capacity[i] = volumetric_heat_capacity * (layer_depths[i+1] - layer_depths[i-1]) / 2.0
            thermal_conductance[i] = bulk_thermal_conductivity[i] / (layer_depths[i+1] - layer_depths[i])
        end
    end

    # Solar radiation
    Q_solar = absorptivity * solar_radiation * (1.0 - shade)
    if slope > 0 && zenith_angle < 90u"°"
        cos_zenith = cosd(zenith_angle)
        cos_slope_zenith = cosd(slope_zenith_angle)
        Q_solar = (Q_solar / cos_zenith) * cos_slope_zenith
    end

    # Longwave radiation
    longwave_out = longwave_radiation(; micro_terrain, surface_temperature=temperature_vector[1], environment_instant)
    Q_infrared = longwave_out.net_longwave_radiation

    # Conduction
    Q_conduction = thermal_conductance[1] * (temperature_vector[2] - temperature_vector[1])

    # Convection
    log_z_ratio = log(reference_height / roughness_height + 1)
    surface_temperature = temperature_vector[1]
    ΔT = air_temperature - surface_temperature
    mean_temperature = (surface_temperature + air_temperature) / 2
    # TODO call calc_ρ_cp method specific to elevation and RH in final version but do it this way for NicheMapR comparison
    ρ_cp = calc_ρ_cp(mean_temperature)
    if air_temperature ≥ surface_temperature || zenith_angle ≥ 90°
        u_star = calc_u_star(; reference_wind_speed=wind_speed, log_z_ratio, κ=karman_constant)
        convective_heat_flux = calc_convection(; u_star, log_z_ratio, ΔT, ρ_cp, z0=roughness_height)
    else
        # compute ρcpTκg (was a constant in original Fortran version)
        ρcpTκg = 6.003e-8u"cal*minute^2/cm^4"
        Obukhov_out = calc_Obukhov_length(air_temperature, surface_temperature, wind_speed, roughness_height, reference_height, ρcpTκg, karman_constant, log_z_ratio, ΔT, ρ_cp)
        convective_heat_flux = Obukhov_out.convective_heat_flux
    end
    heat_transfer_coefficient = max(abs(convective_heat_flux / (temperature_vector[1] - air_temperature)), 0.5u"W/m^2/K")

    # Evaporation
    wet_air_out = wet_air_properties(u"K"(air_temperature), relative_humidity, atmospheric_pressure)
    air_heat_capacity = wet_air_out.specific_heat
    air_density = wet_air_out.density
    mass_transfer_coefficient = (heat_transfer_coefficient / (air_heat_capacity * air_density)) * (0.71 / 0.60)^0.666
    Q_evaporation, evaporation_mass_flux = evaporation(;
        surface_temperature=u"K"(temperature_state[1]),
        air_temperature=u"K"(air_temperature),
        relative_humidity,
        surface_relative_humidity=1.0,
        mass_transfer_coefficient,
        atmospheric_pressure,
        soil_wetness,
        saturated=false,
    )
    # Construct static vector of change in soil temperature, to return
    # Energy balance at surface
    surface = u"K/minute"((Q_solar + Q_infrared + Q_conduction + convective_heat_flux - Q_evaporation) / heat_capacity[1])
    # Soil conduction for internal nodes
    internal = ntuple(Val{N-2}()) do i
        u"K/minute".((thermal_conductance[i] * (temperature_vector[i] - temperature_vector[i+1]) + thermal_conductance[i+1] * (temperature_vector[i+2] - temperature_vector[i+1])) / heat_capacity[i+1])
    end
    # Lower boundary condition
    lower_boundary = 0.0u"K/minute"
    dT = SVector{N}((surface, internal..., lower_boundary))

    return dT
end

function interpolate_forcings(f, t)
    t_m = ustrip(u"minute", t)
    return (;
        atmospheric_pressure = f.interpolate_pressure(t_m),
        air_temperature = f.interpolate_temperature(t_m),
        wind_speed = max(0.1u"m/s", f.interpolate_wind(t_m)),
        zenith_angle = min(90.0u"°", u"°"(round(f.interpolate_zenith(t_m), digits=3))),
        solar_radiation = max(0.0u"W/m^2", f.interpolate_solar(t_m)),
        cloud_cover = clamp(f.interpolate_cloud(t_m), 0.0, 1.0),
        relative_humidity = clamp(f.interpolate_humidity(t_m), 0.0, 1.0),
        slope_zenith_angle = min(90.0u"°", f.interpolate_slope_zenith(t_m)),
    )
end

function evaporation(;
    surface_temperature,
    air_temperature,
    relative_humidity,
    surface_relative_humidity,
    mass_transfer_coefficient,
    atmospheric_pressure,
    soil_wetness,
    saturated,
)
    # Assumes all units are SI (Kelvin, Pascal, meters, seconds, kg, etc.)

    clamped_surface_temperature = surface_temperature < u"K"(-81.0u"°C") ? u"K"(-81.0u"°C") : surface_temperature

    # surface and air vapor densities
    surface_vapour_density = wet_air_properties(u"K"(clamped_surface_temperature), surface_relative_humidity, atmospheric_pressure).vapour_density
    air_vapour_density = wet_air_properties(u"K"(air_temperature), relative_humidity, atmospheric_pressure).vapour_density

    # Effective surface wetness fraction
    surface_wetness = saturated ? 1.0 : soil_wetness

    # Water evaporated from surface (kg/s/m^2)
    water_flux = surface_wetness * mass_transfer_coefficient * (surface_vapour_density - air_vapour_density)

    # Latent heat of vaporization (J/kg)
    latent_heat_vaporisation = enthalpy_of_vaporisation(clamped_surface_temperature)

    # Energy flux due to evaporation (W/m²)
    Q_evaporation = water_flux * latent_heat_vaporisation

    # Mass flux (g/s/m²)
    evaporation_mass_flux = u"g/s/m^2"(water_flux)

    # No water loss if surface temperature ≤ 0°C (e.g., melting snow only)
    if surface_temperature <= u"K"(0.0u"°C")
        evaporation_mass_flux = 0.0u"g/s/m^2"
    end

    return Q_evaporation, evaporation_mass_flux
end

function allocate_soil_water_balance(num_layers)
    (;
        water_potential = zeros(typeof(0.0u"J/kg"), num_layers+1),
        depth = zeros(typeof(0.0u"m"), num_layers+1),
        layer_water_mass = zeros(typeof(0.0u"kg/m^2"), num_layers+1),
        water_content = zeros(typeof(0.0u"m^3/m^3"), num_layers+1),
        water_content_new = zeros(typeof(0.0u"m^3/m^3"), num_layers+1),
        hydraulic_conductivity = zeros(typeof(0.0u"kg*s/m^3"), num_layers+1),
        soil_humidity = zeros(typeof(0.0), num_layers+1),
        soil_temperature = zeros(typeof(0.0u"K"), num_layers+1),
        root_water_potential = zeros(typeof(0.0u"J/kg"), num_layers+1),
        air_entry_potential = zeros(typeof(0.0u"J/kg"), num_layers+1),
        campbell_b_inverse = zeros(typeof(0.0), num_layers+1),
        campbell_exponent = zeros(typeof(0.0), num_layers+1),
        campbell_exponent_complement = zeros(typeof(0.0), num_layers+1),
        saturation_water_content = zeros(typeof(0.0), num_layers+1),
        root_resistance = zeros(typeof(0.0u"m^4/kg/s"), num_layers+1),
        root_zone_parameter = zeros(typeof(0.0u"m"), num_layers+1),
        vapor_flux = zeros(typeof(0.0u"kg/m^2/s"), num_layers+1),
        vapor_flux_derivative = zeros(typeof(0.0u"kg*s/m^4"), num_layers+1),
        hydraulic_capacitance = zeros(typeof(0.0u"kg*s/m^4"), num_layers+1),
        sub_diagonal = zeros(typeof(0.0u"kg*s/m^4"), num_layers+1),
        diagonal = zeros(typeof(0.0u"kg*s/m^4"), num_layers+1),
        super_diagonal = zeros(typeof(0.0u"kg*s/m^4"), num_layers+1),
        normalized_super_diagonal = zeros(typeof(0.0), num_layers+1),
        mass_balance_residual = zeros(typeof(0.0u"kg/m^2/s"), num_layers+1),
        normalized_residual = zeros(typeof(0.0u"J/kg"), num_layers+1),
        potential_change = zeros(typeof(0.0u"J/kg"), num_layers+1),
        soil_resistance = zeros(typeof(0.0u"m^4/kg/s"), num_layers+1),
        root_water_uptake = zeros(typeof(0.0u"kg/m^2/s"), num_layers+1),
        # Output buffers
        water_potential_out = Vector{typeof(1.0u"J/kg")}(undef, num_layers),
        soil_humidity_out = Vector{Float64}(undef, num_layers),
        root_water_potential_out = Vector{typeof(1.0u"J/kg")}(undef, num_layers),
    )
end

function soil_water_balance!(buffers, smm::SoilMoistureModel;
    depths,
    atmospheric_pressure,
    soil_moisture,
    local_relative_humidity,
    leaf_area_index,
    evapotranspiration,
    input_soil_temperature,
)
    # Local variable names
    θ_soil = soil_moisture
    lai = leaf_area_index
    relative_humidity_local = local_relative_humidity

    dt = smm.moist_step
    saturated_conductivity = smm.saturated_hydraulic_conductivity
    campbell_b = smm.campbell_b_parameter
    bulk_density = smm.soil_bulk_density2
    mineral_density = smm.soil_mineral_density2
    root_density = smm.root_density
    root_resistance_param = smm.root_resistance
    stomatal_closure_potential = smm.stomatal_closure_potential
    leaf_resistance = smm.leaf_resistance
    stomatal_stability = smm.stomatal_stability_parameter
    root_radius = smm.root_radius
    moist_error = smm.moist_error
    moist_count = smm.moist_count

    (; water_potential, depth, layer_water_mass, water_content, water_content_new,
       hydraulic_conductivity, soil_humidity, soil_temperature,
       root_water_potential, air_entry_potential,
       campbell_b_inverse, campbell_exponent, campbell_exponent_complement,
       saturation_water_content, root_resistance, root_zone_parameter,
       vapor_flux, vapor_flux_derivative, hydraulic_capacitance,
       sub_diagonal, diagonal, super_diagonal, normalized_super_diagonal,
       mass_balance_residual, normalized_residual, potential_change,
       soil_resistance, root_water_uptake) = buffers
    num_layers = length(water_potential) - 1

    # Constants
    water_molar_mass = 0.01801528u"kg/mol"
    water_density = 1000.0u"kg/m^3"
    vapor_diffusivity = 2.4e-5u"m^2/s"

    # Convert to negative absolute value; fill positions 1..num_layers, extend boundary
    air_entry_potential[1:num_layers] .= -abs.(smm.air_entry_water_potential)
    air_entry_potential[num_layers+1] = air_entry_potential[num_layers]

    # Saturation water content, m3/m3; extend boundary
    saturation_water_content[1:num_layers] .= 1.0 .- bulk_density ./ mineral_density
    saturation_water_content[num_layers+1] = saturation_water_content[num_layers]

    # Soil hydraulic properties; extend boundary
    campbell_b_inverse[1:num_layers] .= 1.0 ./ campbell_b
    campbell_b_inverse[num_layers+1] = campbell_b_inverse[num_layers]
    campbell_exponent[1:num_layers] .= 2.0 .+ 3.0 ./ campbell_b
    campbell_exponent[num_layers+1] = campbell_exponent[num_layers]
    campbell_exponent_complement[1:num_layers] .= 1.0 .- campbell_exponent[1:num_layers]
    campbell_exponent_complement[num_layers+1] = campbell_exponent_complement[num_layers]

    # Fill depth directly from user-provided fine-resolution depths vector
    for i in 1:num_layers
        depth[i+1] = uconvert(u"m", depths[i])
    end
    depth[1] = 0.0u"m"
    depth[2] = 0.0u"m"

    # Copy soil temperature directly (depths already at fine resolution)
    for i in 1:num_layers
        soil_temperature[i] = input_soil_temperature[i]
    end
    soil_temperature[num_layers+1] = input_soil_temperature[num_layers]

    # Set initial water content and related variables
    for i in 2:num_layers
        water_content_new[i] = θ_soil[i-1]
        water_potential[i] = air_entry_potential[i] * (saturation_water_content[i] / water_content_new[i])^campbell_b[i] # matric water potential, EQ5.9 (note thetas=water_content are inverted so not raised to -campbell_b)
        soil_humidity[i] = exp(water_molar_mass * water_potential[i] / (R * soil_temperature[i-1])) # soil humidity, EQ5.14
        hydraulic_conductivity[i] = saturated_conductivity[i] * (air_entry_potential[i] / water_potential[i])^campbell_exponent[i] # hydraulic conductivity, EQ6.14
        water_content[i] = θ_soil[i-1] # water content
    end

    # Bulk water mass per soil layer
    for i in 2:num_layers
        layer_water_mass[i] = water_density * (depth[i+1] - depth[i-1]) / 2 # bulk density x volume per unit area, kg/m²
    end
    # Lower boundary condition set to saturated (stays constant)
    water_potential[num_layers+1] = air_entry_potential[num_layers] * (saturation_water_content[num_layers+1] / saturation_water_content[num_layers+1])^campbell_b[num_layers] # water potential
    soil_humidity[num_layers+1] = 1.0 # soil humidity
    water_content[num_layers+1] = saturation_water_content[num_layers+1] # water content
    water_content_new[num_layers+1] = saturation_water_content[num_layers+1] # water content
    depth[1] = -1e10u"m" # depth at node 1, m
    depth[num_layers+1] = 1e20u"m" # depth at deepest node, m
    hydraulic_conductivity[num_layers+1] = saturated_conductivity[num_layers] * (air_entry_potential[num_layers] / water_potential[num_layers+1])^campbell_exponent[num_layers+1] # lower boundary conductivity

    # Initialize root water uptake variables
    for i in 2:num_layers
        if root_density[i] > 0.0u"m/m^3"
            root_resistance[i] = root_resistance_param / (root_density[i] * (depth[i+1] - depth[i-1]) / 2.0)
            root_zone_parameter[i] = campbell_exponent_complement[i] * log(π * root_radius^2 * root_density[i]) / (4.0 * π * root_density[i] * (depth[i+1] - depth[i-1]) / 2.0)
        else
            root_resistance[i] = 1e20u"m^4/kg/s"
            root_zone_parameter[i] = 0.0u"m"
        end
    end

    water_potential[1] = water_potential[2]
    hydraulic_conductivity[1] = 0.0u"kg*s/m^3"

    # Evapotranspiration
    evaporation_potential = exp(-0.82 * ustrip(lai)) * evapotranspiration # partition potential evaporation from potential evapotranspiration, EQ12.30
    transpiration_potential = evapotranspiration - evaporation_potential # now get potential transpiration

    # Plant water uptake
    potential_sum = 0.0u"J*s/m^4"  # numerator of first term on left of EQ11.18, J * s / m⁴
    resistance_sum = 0.0u"kg*s/m^4" # weighted mean root-soil resistance, R_bar, m4 /(s kg)
    leaf_water_potential = 0.0u"J/kg"
    for i in 2:num_layers
        soil_resistance[i] = root_zone_parameter[i] / hydraulic_conductivity[i] # soil resistance, simplification of EQ11.14, assuming conductivity constant in the rhizosphere
        potential_sum += water_potential[i] / (soil_resistance[i] + root_resistance[i]) # summing over layers
        resistance_sum += 1.0 / (soil_resistance[i] + root_resistance[i]) # summing over layers
    end
    mean_soil_potential = potential_sum / resistance_sum # final step in evaluating psi_bar, weighted mean soil water potential, first term on right in EQ11.18
    mean_resistance = (1.0 / resistance_sum) # denominator of first and second terms on right in EQ11.18

    # Newton-Raphson to estimate leaf_water_potential
    counter = 0
    while counter < moist_count
        if leaf_water_potential > mean_soil_potential
            # Seems we need to force the units here or leaf_water_potential is type unstable in the loop
            leaf_water_potential = uconvert(u"J/kg", mean_soil_potential - transpiration_potential * (mean_resistance + leaf_resistance)) # variation on EQ11.18
        end
        stomatal_closure_factor = (leaf_water_potential / stomatal_closure_potential)^stomatal_stability # part of EQ12.28 determining stomatal closure
        stomatal_slope = transpiration_potential * (mean_resistance + leaf_resistance) * stomatal_stability * stomatal_closure_factor / (leaf_water_potential * (1.0 + stomatal_closure_factor)^2) - 1.0 # derivative of stomatal function
        transpiration_residual = mean_soil_potential - leaf_water_potential - transpiration_potential * (mean_resistance + leaf_resistance) / (1.0 + stomatal_closure_factor) # transpiration mass balance (variation on EQ11.18)
        leaf_water_potential = uconvert(u"J/kg", leaf_water_potential - (transpiration_residual / stomatal_slope))
        counter += 1
        if abs(transpiration_residual) <= 10.0u"J/kg"
            break
        end
    end
    stomatal_closure_factor = (leaf_water_potential / stomatal_closure_potential)^stomatal_stability
    actual_transpiration = transpiration_potential / (1.0 + stomatal_closure_factor)
    for i in 2:num_layers
        root_water_uptake[i] = (water_potential[i] - leaf_water_potential - leaf_resistance * actual_transpiration) / (root_resistance[i] + soil_resistance[i]) # root water uptake, EQ11.15
    end

    # Convergence loop
    counter = 0
    while counter < moist_count
        mass_balance_error = 0.0u"kg/m^2/s"
        counter += 1
        for i in 2:num_layers
            hydraulic_conductivity[i] = saturated_conductivity[i] * (air_entry_potential[i] / water_potential[i])^campbell_exponent[i]
        end

        vapor_flux[1] = evaporation_potential * (soil_humidity[2] - relative_humidity_local) / (1.0 - relative_humidity_local) # vapour flux at soil surface, EQ9.14
        vapor_flux_derivative[1] = evaporation_potential * water_molar_mass * soil_humidity[2] / (R * soil_temperature[1] * (1.0 - relative_humidity_local)) # derivative of vapour flux at soil surface, combination of EQ9.14 and EQ5.14

        for i in 2:num_layers
            vapor_density = wet_air_properties(u"K"(soil_temperature[i]), 1.0, atmospheric_pressure).vapour_density # vapor_density is vapour density = c'_v in EQ9.7
            vapor_conductivity = 0.66 * vapor_diffusivity * vapor_density * (saturation_water_content[i] - (water_content_new[i] + water_content_new[i+1]) / 2.0) / (depth[i+1] - depth[i]) # vapour conductivity, EQ9.7, assuming epsilon(psi_g) = b*psi_g^m (eq. 3.10) where b = 0.66 and m = 1 (p.99)
            vapor_flux[i] = vapor_conductivity * (soil_humidity[i+1] - soil_humidity[i]) # fluxes of vapour within soil, EQ9.14
            vapor_flux_derivative[i] = water_molar_mass * soil_humidity[i] * vapor_conductivity / (R * soil_temperature[i-1]) # derivatives of vapour fluxes within soil, combination of EQ9.14 and EQ5.14
            hydraulic_capacitance[i] = -1.0 * layer_water_mass[i] * water_content_new[i] / (campbell_b[i] * water_potential[i] * dt) # hydraulic capacity = capacitance, d_theta/d_psi
            # Tridiagonal matrix components
            sub_diagonal[i] = -1.0 * hydraulic_conductivity[i-1] / (depth[i] - depth[i-1]) + Unitful.gn * campbell_exponent[i] * hydraulic_conductivity[i-1] / water_potential[i-1]
            super_diagonal[i] = -1.0 * hydraulic_conductivity[i+1] / (depth[i+1] - depth[i])
            diagonal[i] = hydraulic_conductivity[i] / (depth[i] - depth[i-1]) + hydraulic_conductivity[i] / (depth[i+1] - depth[i]) + hydraulic_capacitance[i] - Unitful.gn * campbell_exponent[i] * hydraulic_conductivity[i] / water_potential[i] + vapor_flux_derivative[i-1] + vapor_flux_derivative[i]
            # mass balance including vapour fluxes and root water uptake
            # version of equation 8.28 that additionally contains vapour fluxes and root water uptake
            mass_balance_residual[i] = ((water_potential[i] * hydraulic_conductivity[i] - water_potential[i-1] * hydraulic_conductivity[i-1]) / (depth[i] - depth[i-1]) - (water_potential[i+1] * hydraulic_conductivity[i+1] - water_potential[i] * hydraulic_conductivity[i]) / (depth[i+1] - depth[i])) / campbell_exponent_complement[i] + layer_water_mass[i] * (water_content_new[i] - water_content[i]) / dt - Unitful.gn * (hydraulic_conductivity[i-1] - hydraulic_conductivity[i]) + vapor_flux[i-1] - vapor_flux[i] + root_water_uptake[i]
            mass_balance_error += abs(mass_balance_residual[i]) # total mass balance error
        end

        # Thomas algorithm (Gauss elimination)
        for i in 2:num_layers-1
            normalized_super_diagonal[i] = super_diagonal[i] / diagonal[i]
            normalized_residual[i] = mass_balance_residual[i] / diagonal[i]
            diagonal[i+1] -= sub_diagonal[i+1] * normalized_super_diagonal[i]
            mass_balance_residual[i+1] -= sub_diagonal[i+1] * normalized_residual[i]
        end

        potential_change[num_layers] = mass_balance_residual[num_layers] / diagonal[num_layers]
        water_potential[num_layers] -= potential_change[num_layers]
        water_potential[num_layers] = min(water_potential[num_layers], air_entry_potential[num_layers])

        for i in (num_layers-1):-1:2
            potential_change[i] = normalized_residual[i] - normalized_super_diagonal[i] * potential_change[i+1] # change in matric potential in an iteration step, J/kg
            water_potential[i] -= potential_change[i] # matric potential, J/kg
            if water_potential[i] > air_entry_potential[i]
                water_potential[i] = (water_potential[i] + potential_change[i] + air_entry_potential[i]) / 2.0
            end
        end

        for i in 2:num_layers
            water_content_new[i] = max(saturation_water_content[i] * (air_entry_potential[i] / water_potential[i])^campbell_b_inverse[i], 1e-7)
            water_potential[i] = air_entry_potential[i] * (saturation_water_content[i] / water_content_new[i])^campbell_b[i]
            soil_humidity[i] = exp(water_molar_mass * water_potential[i] / (R * soil_temperature[i-1]))
        end
        soil_humidity[num_layers+1] = soil_humidity[num_layers]
        if mass_balance_error <= moist_error
            break
        end
    end

    surface_water_flux = ((water_potential[2] * hydraulic_conductivity[2] - water_potential[3] * hydraulic_conductivity[3]) / (campbell_exponent_complement[2] * (depth[3] - depth[2])) + Unitful.gn * hydraulic_conductivity[2] + actual_transpiration) * dt
    water_content .= water_content_new
    for i in 2:num_layers+1
        θ_soil[i-1] = water_content_new[i]
    end

    for i in 2:num_layers
        root_water_potential[i] = -1.0 * (actual_transpiration * soil_resistance[i] - water_potential[i])
    end

    evap = evaporation_potential * (soil_humidity[2] - relative_humidity_local) / (1.0 - relative_humidity_local) * dt
    return (;
        evaporation = evap,
        transpiration = actual_transpiration,
        soil_moisture = θ_soil,
        leaf_water_potential,
        # These need the first value removed. Why?
        soil_water_potential = (buffers.water_potential_out .= view(water_potential, 2:(num_layers+1))),
        root_water_potential = (buffers.root_water_potential_out .= view(root_water_potential, 2:(num_layers+1))),
        soil_humidity = (buffers.soil_humidity_out .= view(soil_humidity, 2:(num_layers+1))),
        surface_water_flux,
        drainage = Unitful.gn * hydraulic_conductivity[num_layers]
    )
end


get_soil_water_balance(soil_moisture_model; num_layers=18, kw...) =
    get_soil_water_balance!(allocate_soil_water_balance(num_layers), soil_moisture_model; kw...)

function get_soil_water_balance!(buffers, soil_moisture_model::SoilMoistureModel;
    depths,
    micro_terrain,
    environment_instant,
    T0,
    pool,
    niter_moist,
    soil_wetness,
    soil_moisture,
)
    air_temperature = environment_instant.reference_temperature
    atmospheric_pressure = environment_instant.atmospheric_pressure
    relative_humidity = environment_instant.reference_humidity
    leaf_area_index = environment_instant.leaf_area_index

    (; maxpool, moist_step, soil_bulk_density2, soil_mineral_density2) = soil_moisture_model

    bulk_density = soil_bulk_density2
    mineral_density = soil_mineral_density2
    θ_soil = soil_moisture
    surface_temperature = T0[1]

    # compute scalar profiles
    profile_out = atmospheric_surface_profile!(buffers.profile;
        micro_terrain, environment_instant, surface_temperature,
    )

    # convection
    convective_heat_flux = profile_out.convective_heat_flux

    # evaporation
    wet_air_out_ref = wet_air_properties(u"K"(last(profile_out.air_temperature)), last(profile_out.relative_humidity), atmospheric_pressure)
    wet_air_out_loc = wet_air_properties(u"K"(profile_out.air_temperature[1]), 1.0, atmospheric_pressure)
    local_relative_humidity = clamp(wet_air_out_ref.vapour_pressure / wet_air_out_loc.vapour_pressure, 0.0, 0.99)
    heat_transfer_coefficient = max(abs(convective_heat_flux / (surface_temperature - air_temperature)), 0.5u"W/m^2/K")
    wet_air_out = wet_air_properties(air_temperature, relative_humidity, atmospheric_pressure)
    air_heat_capacity = wet_air_out.specific_heat
    air_density = wet_air_out.density
    mass_transfer_coefficient = (heat_transfer_coefficient / (air_heat_capacity * air_density)) * (0.71 / 0.60)^0.666
    Q_evaporation, evaporation_mass_flux = evaporation(;
        surface_temperature,
        air_temperature,
        relative_humidity,
        surface_relative_humidity=1.0,
        mass_transfer_coefficient,
        atmospheric_pressure,
        soil_wetness,
        saturated=true,
    )
    latent_heat_vaporisation = enthalpy_of_vaporisation(surface_temperature)
    evaporation_potential = max(1e-7u"kg/m^2/s", Q_evaporation / latent_heat_vaporisation)

    # run infiltration algorithm
    infil_out = soil_water_balance!(buffers.soil_water_balance, soil_moisture_model;
        depths,
        atmospheric_pressure,
        local_relative_humidity,
        leaf_area_index,
        soil_moisture,
        evapotranspiration=evaporation_potential,
        input_soil_temperature=T0,
    )
    soil_moisture = infil_out.soil_moisture
    surf_evap = max(0.0u"kg/m^2", infil_out.evaporation)
    water_flux = max(0.0u"kg/m^2", infil_out.surface_water_flux)
    pool = clamp(pool - water_flux - surf_evap, 0.0u"kg/m^2", maxpool) # pooling surface water
    if pool > 0.0u"kg/m^2" # surface is wet - saturate it for infiltration
        soil_moisture[1] = 1 - bulk_density[1] / mineral_density[1]
    end
    for _ in 1:(niter_moist-1)
        infil_out = soil_water_balance!(buffers.soil_water_balance, soil_moisture_model;
            depths,
            atmospheric_pressure,
            local_relative_humidity,
            soil_moisture,
            leaf_area_index,
            evapotranspiration=evaporation_potential,
            input_soil_temperature=T0,
        )
        soil_moisture = infil_out.soil_moisture
        surf_evap = max(0.0u"kg/m^2", infil_out.evaporation)
        water_flux = max(0.0u"kg/m^2", infil_out.surface_water_flux)
        pool = clamp(pool - water_flux - surf_evap, 0.0u"kg/m^2", maxpool)
        if pool > 0.0u"kg/m^2"
            soil_moisture[1] = 1 - bulk_density[1] / mineral_density[1]
        end
    end
    soil_wetness = clamp(abs(surf_evap / (evaporation_potential * moist_step) * 1.0), 0, 1.0)

    return (; infil_out, soil_wetness, pool, soil_moisture)
end

function allocate_phase_transition(num_nodes)
    layer_mass = zeros(Float64, num_nodes)u"kg"
    phase_change_heat = zeros(Float64, num_nodes)u"J"
    return (; layer_mass, phase_change_heat)
end

phase_transition(; depths, kw...) =
    phase_transition!(allocate_phase_transition(length(depths)); depths, kw...)

function phase_transition!(
    buffers::NamedTuple;
    temperatures::AbstractVector,
    temperatures_past::AbstractVector,
    accumulated_latent_heat::AbstractVector,
    soil_moisture::AbstractVector,
    depths::AbstractVector,
)
    (; layer_mass, phase_change_heat) = buffers
    latent_heat_fusion = 333550.0u"J/kg"
    specific_heat_water = 4184.0u"J/kg/K"
    num_nodes = length(depths)
    mean_temperature = similar(temperatures)
    mean_temperature_past = similar(temperatures)
    temperature = MVector(temperatures)
    tolerance = 1.0e-4u"°C"

    for i in 1:num_nodes
        # --- Always compute mean layer temperatures ---
        if i < num_nodes
            mean_temperature[i] = 0.5 * (temperature[i] + temperature[i+1])
            mean_temperature_past[i] = 0.5 * (temperatures_past[i] + temperatures_past[i+1])
        else
            mean_temperature[i] = temperature[i]
            mean_temperature_past[i] = temperatures_past[i]
        end
        # --- Compute layer mass (kg), handle dry layers ---
        if soil_moisture[i] > 0
            if i < num_nodes
                layer_mass[i] = (u"m"(depths[i+1] - depths[i])) * 1000.0u"kg/m" * soil_moisture[i]
            else
                layer_mass[i] = (u"m"(depths[i] + 100.0u"cm" - depths[i])) * 1000.0u"kg/m" * soil_moisture[i]
            end
        else
            layer_mass[i] = 0.0u"kg"
        end
        max_latent_heat = latent_heat_fusion * layer_mass[i]

        # --- If no water, reset and skip ---
        if soil_moisture[i] <= 0.0 || layer_mass[i] <= 0.0u"kg"
            accumulated_latent_heat[i] = 0.0u"J"
            phase_change_heat[i] = 0.0u"J"
        end

        # ==============================
        #  PHASE CHANGE CALCULATIONS
        # ==============================

        # --- FREEZING (above → below 0°C) ---
        if (mean_temperature_past[i] > tolerance) && (mean_temperature[i] <= -tolerance)
            phase_change_heat[i] = (mean_temperature_past[i] - mean_temperature[i]) * layer_mass[i] * specific_heat_water
            accumulated_latent_heat[i] += phase_change_heat[i]
            if accumulated_latent_heat[i] >= max_latent_heat
                accumulated_latent_heat[i] = max_latent_heat
                phase_change_heat[i] = 0.0u"J"
            end

            temperature[i] = 0.0u"°C"
            if i < num_nodes
                temperature[i+1] = 0.0u"°C"
            end

        # --- THAWING (below → above 0°C) ---
        elseif (mean_temperature_past[i] < -tolerance) && (mean_temperature[i] >= tolerance)
            phase_change_heat[i] = (mean_temperature[i] - mean_temperature_past[i]) * layer_mass[i] * specific_heat_water
            accumulated_latent_heat[i] -= phase_change_heat[i]

            if accumulated_latent_heat[i] <= 0.0u"J"
                accumulated_latent_heat[i] = 0.0u"J"
                phase_change_heat[i] = 0.0u"J"
            else
                temperature[i] = 0.0u"°C"
                if i < num_nodes
                    temperature[i+1] = 0.0u"°C"
                end
            end

        else
            phase_change_heat[i] = 0.0u"J"
        end
    end
    return (; accumulated_latent_heat, phase_change_heat, temperature=SVector(temperature))
end
