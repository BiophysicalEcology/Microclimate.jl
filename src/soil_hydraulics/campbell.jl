"""
    CampbellSoilHydraulics(; ...)

Soil hydraulic parameters for Campbell's (1985) soil water balance model.
Holds parameters only — the moisture-solver strategy (`PrescribedSoilMoisture` /
`DynamicSoilMoisture`) and tuning (`moisture_tolerance`, `moisture_max_iterations`,
`moisture_timestep`, `max_surface_pool`) live on `MicroConfig`, since they are
solver/strategy choices, not Campbell-specific.

# References
Campbell, G. S. (1985). Soil Physics with BASIC. Elsevier.
"""
@kwdef struct CampbellSoilHydraulics{AEWP,SHC,CBP,BD,MD,RDen,RRes,SCP,LRes,SSP,RRad} <: AbstractSoilHydraulicsModel
    air_entry_water_potential::AEWP
    saturated_hydraulic_conductivity::SHC
    campbell_b_parameter::CBP
    bulk_density::BD       # per-depth profile (Vector); also consumed by the soil thermal model
    mineral_density::MD    # per-depth profile (Vector); also consumed by the soil thermal model
    root_density::RDen
    root_resistance::RRes
    stomatal_closure_potential::SCP
    leaf_resistance::LRes
    stomatal_stability_parameter::SSP
    root_radius::RRad
end

function allocate_soil_water_balance(::CampbellSoilHydraulics, num_layers)
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

function infiltration_step!(buffers, smm::CampbellSoilHydraulics;
    depths,
    atmospheric_pressure,
    soil_moisture,
    local_relative_humidity,
    leaf_area_index,
    evapotranspiration,
    input_soil_temperature,
    moisture_timestep,    # solver tuning, lives on MicroConfig
    moisture_tolerance,
    moisture_max_iterations,
    vapour_pressure_equation=GoffGratch(),
)
    # Local variable names
    θ_soil = soil_moisture
    lai = leaf_area_index
    relative_humidity_local = local_relative_humidity

    dt = moisture_timestep
    saturated_conductivity = smm.saturated_hydraulic_conductivity
    campbell_b = smm.campbell_b_parameter
    bulk_density = smm.bulk_density
    mineral_density = smm.mineral_density
    root_density = smm.root_density
    root_resistance_param = smm.root_resistance
    stomatal_closure_potential = smm.stomatal_closure_potential
    leaf_resistance = smm.leaf_resistance
    stomatal_stability = smm.stomatal_stability_parameter
    root_radius = smm.root_radius

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
    while counter < moisture_max_iterations
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
    while counter < moisture_max_iterations
        mass_balance_error = 0.0u"kg/m^2/s"
        counter += 1
        for i in 2:num_layers
            hydraulic_conductivity[i] = saturated_conductivity[i] * (air_entry_potential[i] / water_potential[i])^campbell_exponent[i]
        end

        vapor_flux[1] = evaporation_potential * (soil_humidity[2] - relative_humidity_local) / (1.0 - relative_humidity_local) # vapour flux at soil surface, EQ9.14
        vapor_flux_derivative[1] = evaporation_potential * water_molar_mass * soil_humidity[2] / (R * soil_temperature[1] * (1.0 - relative_humidity_local)) # derivative of vapour flux at soil surface, combination of EQ9.14 and EQ5.14

        # Campbell (1985) Soil Physics with BASIC, Eq. 3.10: ε(ψ_g) = b·ψ_g^m with b=0.66, m=1 (p. 99).
        # Matches Fortran INFIL.f:236 literal `KV = 0.66 * DV * ...` — keep the two in sync if m ≠ 1 is ever used.
        tortuosity_b = 0.66
        for i in 2:num_layers
            vapor_density = wet_air_properties(u"K"(soil_temperature[i]), 1.0, atmospheric_pressure; vapour_pressure_equation).vapour_density # vapor_density is vapour density = c'_v in EQ9.7
            vapor_conductivity = tortuosity_b * vapor_diffusivity * vapor_density * (saturation_water_content[i] - (water_content_new[i] + water_content_new[i+1]) / 2.0) / (depth[i+1] - depth[i]) # vapour conductivity, EQ9.7
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
        # Fortran INFIL.f lines 255-256: clamp tiny normalized super-diagonal values to zero
        for i in 2:num_layers-1
            normalized_super_diagonal[i] = super_diagonal[i] / diagonal[i]
            if abs(normalized_super_diagonal[i]) < 1e-8
                normalized_super_diagonal[i] = zero(normalized_super_diagonal[i])
            end
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
        if mass_balance_error <= moisture_tolerance
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


soil_water_balance(soil_hydraulics::CampbellSoilHydraulics; num_layers=18, kw...) =
    soil_water_balance!(allocate_soil_water_balance(soil_hydraulics, num_layers), soil_hydraulics; kw...)

function soil_water_balance!(buffers, soil_hydraulics::CampbellSoilHydraulics;
    depths,
    site,
    boundary_layer_model,
    environment_instant,
    T0,
    pool,
    niter_moist,
    soil_wetness,
    soil_moisture,
    moisture_timestep,    # solver tuning, lives on MicroConfig
    moisture_tolerance,
    moisture_max_iterations,
    max_surface_pool,
    evaporation_model::AbstractEvaporationModel=KearneyEvaporation(),
    vapour_pressure_equation=GoffGratch(),
    snow_present=false,
)
    air_temperature = environment_instant.reference_temperature
    atmospheric_pressure = environment_instant.atmospheric_pressure
    relative_humidity = environment_instant.reference_humidity
    leaf_area_index = environment_instant.leaf_area_index

    (; bulk_density, mineral_density) = soil_hydraulics

    θ_soil = soil_moisture
    surface_temperature = T0[1]

    # compute scalar profiles
    profile_out = atmospheric_surface_profile!(boundary_layer_model, buffers.profile;
        site, environment_instant, surface_temperature, vapour_pressure_equation,
    )

    # convection
    convective_heat_flux = profile_out.convective_heat_flux

    # evaporation
    wet_air_out_ref = wet_air_properties(u"K"(last(profile_out.air_temperature)), last(profile_out.relative_humidity), atmospheric_pressure; vapour_pressure_equation)
    wet_air_out_loc = wet_air_properties(u"K"(profile_out.air_temperature[1]), 1.0, atmospheric_pressure; vapour_pressure_equation)
    local_relative_humidity = clamp(wet_air_out_ref.vapour_pressure / wet_air_out_loc.vapour_pressure, 0.0, 0.99)
    heat_transfer_coefficient = max(abs(convective_heat_flux / (surface_temperature - air_temperature)), 0.5u"W/m^2/K")
    wet_air_out = wet_air_properties(air_temperature, relative_humidity, atmospheric_pressure; vapour_pressure_equation)
    air_heat_capacity = wet_air_out.specific_heat
    air_density = wet_air_out.density
    mass_transfer_coefficient = (heat_transfer_coefficient / (air_heat_capacity * air_density)) * (0.71 / 0.60)^0.666
    Q_evaporation, evaporation_mass_flux = evaporation(evaporation_model;
        surface_temperature,
        air_temperature,
        relative_humidity,
        surface_relative_humidity=1.0,
        mass_transfer_coefficient,
        atmospheric_pressure,
        soil_wetness,
        saturated=true,
        vapour_pressure_equation,
    )
    latent_heat_vaporisation = enthalpy_of_vaporisation(surface_temperature)
    evaporation_potential = max(1e-7u"kg/m^2/s", Q_evaporation / latent_heat_vaporisation)
    # Fortran OSUB.f lines 1188-1196: suppress soil evaporation when snow covers the ground
    if snow_present
        evaporation_potential = 1e-7u"kg/m^2/s"
    end

    # When pool > 0, top boundary is saturated. Establish that BC before the
    # first infil call so iter 1 uses the same top BC as iters 2..N. Matches
    # Fortran OSUB.f:1198-1206 (commit fbf9a91): refill mass is debited from
    # pool, soil node 1 set to saturation. Half-thickness uses the coarse
    # spacing (depths[3]-depths[1])/2 to match Fortran's (dep(2)-dep(1))/2,
    # since `depths` here is the fine 19-node grid in which depths[3] is the
    # first user-specified depth below the surface.
    sat = 1 - bulk_density[1] / mineral_density[1]
    half_thickness = (depths[3] - depths[1]) / 2
    if pool > 0.0u"kg/m^2"
        refill = max(0.0u"kg/m^2", uconvert(u"kg/m^2", (sat - soil_moisture[1]) * half_thickness * 1000.0u"kg/m^3"))
        pool = max(0.0u"kg/m^2", pool - refill)
        soil_moisture[1] = sat
    end

    # run infiltration algorithm
    infil_out = infiltration_step!(buffers.soil_water_balance, soil_hydraulics;
        depths,
        atmospheric_pressure,
        local_relative_humidity,
        leaf_area_index,
        soil_moisture,
        evapotranspiration=evaporation_potential,
        input_soil_temperature=T0,
        moisture_timestep, moisture_tolerance, moisture_max_iterations,
        vapour_pressure_equation,
    )
    soil_moisture = infil_out.soil_moisture
    surf_evap = max(0.0u"kg/m^2", infil_out.evaporation)
    water_flux = max(0.0u"kg/m^2", infil_out.surface_water_flux)
    pool = clamp(pool - water_flux - surf_evap, 0.0u"kg/m^2", max_surface_pool) # pooling surface water
    if pool > 0.0u"kg/m^2"
        refill = max(0.0u"kg/m^2", uconvert(u"kg/m^2", (sat - soil_moisture[1]) * half_thickness * 1000.0u"kg/m^3"))
        pool = max(0.0u"kg/m^2", pool - refill)
        soil_moisture[1] = sat
    end
    for _ in 1:(niter_moist-1)
        infil_out = infiltration_step!(buffers.soil_water_balance, soil_hydraulics;
            depths,
            atmospheric_pressure,
            local_relative_humidity,
            soil_moisture,
            leaf_area_index,
            evapotranspiration=evaporation_potential,
            input_soil_temperature=T0,
            moisture_timestep, moisture_tolerance, moisture_max_iterations,
            vapour_pressure_equation,
        )
        soil_moisture = infil_out.soil_moisture
        surf_evap = max(0.0u"kg/m^2", infil_out.evaporation)
        water_flux = max(0.0u"kg/m^2", infil_out.surface_water_flux)
        pool = clamp(pool - water_flux - surf_evap, 0.0u"kg/m^2", max_surface_pool)
        if pool > 0.0u"kg/m^2"
            refill = max(0.0u"kg/m^2", uconvert(u"kg/m^2", (sat - soil_moisture[1]) * half_thickness * 1000.0u"kg/m^3"))
            pool = max(0.0u"kg/m^2", pool - refill)
            soil_moisture[1] = sat
        end
    end
    # Fortran OSUB.f line 1239: ptwet = surflux / (ep * timestep) * 100
    # Note Fortran's `surflux` is INFIL's FL output (the humidity-gradient
    # evaporation flux EP*(H(2)-HA)/(1-HA)*DT, INFIL.f:300), NOT the Darcy
    # surface_water_flux (INFIL's SW). In Julia these map to `infil_out.evaporation`
    # and `infil_out.surface_water_flux` respectively.
    # Fortran OSUB.f:1223-1225 clamps `surflux >= 0` before this expression.
    raw_surflux = max(0.0u"kg/m^2", infil_out.evaporation)
    soil_wetness = clamp(raw_surflux / (evaporation_potential * moisture_timestep), 0, 1.0)

    return (; infil_out, soil_wetness, pool, soil_moisture)
end
