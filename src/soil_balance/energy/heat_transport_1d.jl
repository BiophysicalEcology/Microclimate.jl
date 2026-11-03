"""
    SoilHeatTransport1D(; ode_solver=Tsit5(),
                          ode_kwargs=(; reltol=1e-6, abstol=1e-3u"K"),
                          maximum_surface_temperature=85.0u"°C",
                          freezing_model=PhaseTransitionLatentHeat())

One dimension soil column energy balance (as in NicheMapR). The ODE right hand side
computes per-node dT/dt from heat conduction between layers plus a surface
boundary condition combining net radiation (clear-sky + cloud + view
factor), convection (Monin–Obukhov), evaporation latent flux, and a snow
melt freeze correction. The deepest node is pinned at the user-supplied
deep soil temperature (Dirichlet BC).

Solver tuning lives here, not on `MicroConfig`, because each
`SoilHeatTransportModel` chooses its own integration strategy:

- `ode_solver` — any SciML algorithm (`Tsit5()`, `RK4()`, …), explicit or implicit.
- `ode_kwargs` — NamedTuple of kwargs forwarded to the integrator, such as `dt=60u"minute"`
- `maximum_surface_temperature` — surface-temperature safety clamp
  (Fortran microinput[74])
- `freezing_model` — soil ice/water phase-transition correction applied to
  the ODE result each hour
"""
@kwdef struct SoilHeatTransport1D{SOS,SOK,MSF,FM} <: SoilHeatTransportModel
    ode_solver::SOS = Tsit5()
    ode_kwargs::SOK = (; reltol=1e-6, abstol=1e-3u"K")
    maximum_surface_temperature::MSF = 85.0u"°C"
    freezing_model::FM = PhaseTransitionLatentHeat()
end

"""
    HourlySoilEnergyInputs

Forcing and state held fixed over one hour of the soil energy ODE.
"""
@kwdef struct HourlySoilEnergyInputs{D<:Vector{<:Number},EI,SW,LW,QF,SNS,GST,GIL,GA}
    # Node depths, including any active snow nodes
    depths::D
    environment_instant::EI
    soil_wetness::SW
    longwave_sky::LW
    albedo::Float64
    Q_freeze::QF = 0.0u"W/m^2"
    snow_state::SNS = nothing
    # Canopy ground-boundary overrides (nothing under NoCanopy; see apply_canopy_overrides)
    ground_shortwave_transmission::GST = nothing
    ground_incoming_longwave::GIL = nothing
    # The canopy's GroundAirExchange, or nothing without a canopy
    ground_air::GA = nothing
end

"""
    SoilEnergyInputs(model, hourly, forcing, buffers, soil_moisture)

The ODE integrator's `p`: the [`MicroModel`](@ref) parameter tree, the current
[`HourlySoilEnergyInputs`](@ref), the day's forcing interpolants, scratch
buffers and the soil moisture state. Replaced each hour with
[`set_hourly_inputs!`](@ref).
"""
struct SoilEnergyInputs{M,HI<:HourlySoilEnergyInputs,F,B,SM}
    model::M
    hourly::HI
    forcing::F
    buffers::B
    soil_moisture::SM
end

"""
    set_hourly_inputs!(integrator, hourly::HourlySoilEnergyInputs)

Replace the hourly forcing and state in the soil energy ODE `integrator`.
"""
function set_hourly_inputs!(integrator, hourly::HourlySoilEnergyInputs)
    (; model, forcing, buffers, soil_moisture) = integrator.p
    integrator.p = SoilEnergyInputs(model, hourly, forcing, buffers, soil_moisture)
    return integrator
end

# The surface: the first node with any heat capacity
first_active_node(heat_capacity::SVector{N}) where N = something(findfirst(!iszero, heat_capacity), N)

# The thermal properties of soil node `node` at `temperature`
function soil_node_properties(p, node, temperature, atmospheric_pressure, vapour_pressure_equation)
    (; properties_model, profile) = p.model.soil
    properties = soil_properties(maybegetindex(properties_model, node);
        atmospheric_pressure, soil_temperature=temperature, soil_moisture=p.soil_moisture[node],
        bulk_density=profile.bulk_density[node], mineral_density=profile.mineral_density[node],
        mineral_conductivity=profile.mineral_conductivity[node], mineral_heat_capacity=profile.mineral_heat_capacity[node],
        de_vries_shape_factor=profile.de_vries_shape_factor[node], vapour_pressure_equation,
    )
    return (;
        bulk_thermal_conductivity=uconvert(u"W/m/K", properties.bulk_thermal_conductivity),
        bulk_heat_capacity=uconvert(u"J/kg/K", properties.bulk_heat_capacity),
        bulk_density=uconvert(unit(eltype(profile.mineral_density)), properties.bulk_density),
    )
end

# Every soil node's thermal properties at `temperature`
soil_column_properties(p, temperature::SVector{N}, atmospheric_pressure, vapour_pressure_equation) where N =
    SVector(ntuple(node -> soil_node_properties(p, node, temperature[node], atmospheric_pressure, vapour_pressure_equation), Val(N)))

# Every node's thermal properties at `temperature`. The SnowModel method is in snow/snow_model.jl.
node_properties(::Nothing, p, temperature, atmospheric_pressure, vapour_pressure_equation) =
    soil_column_properties(p, temperature, atmospheric_pressure, vapour_pressure_equation)

# Each node's heat capacity per area, and its thermal conductance to the node below
function soil_layers(properties::SVector{N}, depths, n_snow) where N
    depth(i) = i <= N ? uconvert(u"cm", depths[i]) : 0.0u"cm"
    heat_capacity = SVector(ntuple(Val(N)) do i
        volumetric_heat_capacity = properties[i].bulk_density * properties[i].bulk_heat_capacity
        capacity = if i == 1
            volumetric_heat_capacity * depth(2) / 2.0
        # A snow node at depth 0, with snow below it, is the surface and holds the whole
        # distance to the next node
        elseif i <= n_snow && iszero(depth(i)) && !iszero(depth(i + 1))
            volumetric_heat_capacity * depth(i + 1)
        else
            volumetric_heat_capacity * (depth(i + 1) - depth(i - 1)) / 2.0
        end
        uconvert(u"J/K/m^2", capacity)
    end)
    thermal_conductance = SVector(ntuple(Val(N)) do i
        conductivity = properties[i].bulk_thermal_conductivity
        conductance = if i == 1
            iszero(depth(2)) ? zero(conductivity / oneunit(depth(2))) : conductivity / depth(2)
        else
            spacing = depth(i + 1) - depth(i)
            previous_spacing = depth(i) - depth(i - 1)
            # Nodes at the same depth, at the snow-soil interface, use the spacing above
            # them, so the snow stays coupled to the soil
            if iszero(spacing)
                iszero(previous_spacing) ? zero(conductivity / oneunit(spacing)) : conductivity / previous_spacing
            else
                conductivity / spacing
            end
        end
        uconvert(u"W/K/m^2", conductance)
    end)
    return (; heat_capacity, thermal_conductance)
end

# Compile-time snow node count for an ODE's `p.model.snow`. The `::SnowModel{N}`
# method is in snow/snow_model.jl; both ::Nothing and ::NoSnow return 0.
@inline _n_snow(::Nothing) = 0

# Ground shortwave fraction: NoCanopy uses `shade`; MultilayerCanopy substitutes
# its hourly transmission ratio (see apply_canopy_overrides).
ground_shortwave_fraction(shade, ::Nothing) = 1.0 - shade
ground_shortwave_fraction(shade, transmission) = transmission

# Ground longwave terms: NoCanopy passes the shade-based split through;
# MultilayerCanopy supersedes it with its own converged incoming flux.
ground_longwave_terms(incoming, outgoing_coeff, shade_term, surface_emissivity, ::Nothing) =
    (incoming, outgoing_coeff, shade_term)
ground_longwave_terms(incoming, outgoing_coeff, shade_term, surface_emissivity, override) =
    (override, σ * surface_emissivity, 0.0u"W/m^2")

# Ground convection and evaporation. Without a canopy, from the reference air.
ground_heat_vapor_flux(evaporation_model, ::Nothing, ctx) =
    surface_convection_evaporation(evaporation_model; ctx...)
# Under a canopy, with the air and conductances the canopy gives the ground
function ground_heat_vapor_flux(evaporation_model, ground_air::GroundAirExchange, ctx)
    convective_heat_flux = uconvert(u"W/m^2", ground_air.heat_conductance * (ground_air.air_temperature - ctx.surface_temperature))
    Q_evaporation, _ = evaporation(evaporation_model;
        surface_temperature=u"K"(ctx.surface_temperature), air_temperature=u"K"(ground_air.air_temperature),
        relative_humidity=ground_air.air_relative_humidity, surface_relative_humidity=1.0,
        mass_transfer_coefficient=ground_air.vapor_conductance,
        atmospheric_pressure=ctx.atmospheric_pressure, soil_wetness=ctx.soil_wetness,
        saturated=false, vapour_pressure_equation=ctx.vapour_pressure_equation,
    )
    return Q_evaporation, convective_heat_flux
end

function allocate_soil_energy_balance(::SoilHeatTransport1D, num_nodes::Int)
    obukhov_length_prev = Ref(INITIAL_OBUKHOV_LENGTH)  # warm-start across hourly ODE solves
    return (; obukhov_length_prev)
end

# The surface node of state `temperature` at time `t`
function first_active_node(p::SoilEnergyInputs, temperature, t)
    (; atmospheric_pressure) = interpolate_forcings(p.forcing, t)
    snow_model = p.model.snow
    properties = node_properties(snow_model, p, temperature, atmospheric_pressure, p.model.vapour_pressure_equation)
    return first_active_node(soil_layers(properties, p.hourly.depths, _n_snow(snow_model)).heat_capacity)
end

# SciML ODE rhs. Dispatches the actual physics on `p.model.soil.energy_model`.
soil_energy_balance(u, p::SoilEnergyInputs, t) = soil_energy_balance(p.model.soil.energy_model, u, p, t)

function soil_energy_balance(model::SoilHeatTransport1D,
    temperature_state::U,  # state
    p::SoilEnergyInputs,   # "parameters"
    t::Quantity,           # timestep
) where U <: SVector{N} where N
    # extract parameters
    (; hourly, forcing, buffers) = p
    (; site, vapour_pressure_equation) = p.model
    (; heights, boundary_layer_model) = p.model.atmosphere
    atmospheric_radiation_model = p.model.atmosphere.radiation.longwave_model.atmospheric_radiation_model
    evaporation_model = p.model.surface.evaporation_model
    snow_model = p.model.snow
    (; depths, environment_instant, soil_wetness, longwave_sky, albedo, Q_freeze) = hourly
    maximum_surface_temperature = model.maximum_surface_temperature
    (; shade) = environment_instant
    # Get environmental data at time t
    (; atmospheric_pressure, air_temperature, wind_speed, zenith_angle, solar_radiation, cloud_cover, relative_humidity, slope_zenith_angle, diffuse_fraction) = interpolate_forcings(forcing, t)
    (; roughness_height, slope) = site

    reference_height = last(heights)
    absorptivity = 1.0 - albedo

    # check for unstable conditions of ground surface temperature.
    # Fortran microinput(74)=maxsurf (MICROCLIMATE.f:482); −81°C lower bound matches
    # DSUB.f:551.
    T_lo = u"K"(-81.0u"°C")
    T_hi = u"K"(maximum_surface_temperature)
    soil_temperature = clamped_temperature = map(t -> clamp(t, T_lo, T_hi), temperature_state)::U

    # Soil and snow properties at the current temperatures, on every evaluation, as the Fortran
    # (DSUB.f lines 274-337) does
    properties = node_properties(snow_model, p, soil_temperature, atmospheric_pressure, vapour_pressure_equation)
    (; heat_capacity, thermal_conductance) = soil_layers(properties, depths, _n_snow(snow_model))
    j = first_active_node(heat_capacity)

    # Solar radiation. Split into direct + diffuse so that self-shadowed
    # slopes still receive the diffuse sky component. Scaling the entire
    # global flux by `cos_slope_zenith / cos_zenith` (the old behaviour)
    # zeroes both components when the sun is behind the slope, which
    # produced sharp per-pixel sun/shadow streaks at fine DEM resolution.
    Q_solar = absorptivity * solar_radiation * ground_shortwave_fraction(shade, hourly.ground_shortwave_transmission)
    if slope > 0 && zenith_angle < 90u"°"
        cos_zenith = cosd(zenith_angle)
        cos_slope_zenith = cosd(slope_zenith_angle)
        direct_share  = (1.0 - diffuse_fraction) * Q_solar
        diffuse_share = diffuse_fraction * Q_solar
        Q_direct  = max(0.0u"W/m^2", (direct_share / cos_zenith) * cos_slope_zenith)
        Q_diffuse = diffuse_share * site.sky_view_fraction
        Q_solar = Q_direct + Q_diffuse
    end

    # Longwave radiation — recompute from interpolated forcings at each ODE sub-step,
    # matching Fortran DSUB.f lines 396-465 which recomputes at every call.
    # Use interpolated cloud_cover (from line 20), not the hourly snapshot.
    (; surface_emissivity, cloud_emissivity) = environment_instant
    sky_view_fraction = site.sky_view_fraction
    wet_air_out = wet_air_properties(u"K"(air_temperature), relative_humidity, atmospheric_pressure; vapour_pressure_equation)
    atmospheric_longwave = atmospheric_radiation(atmospheric_radiation_model, wet_air_out.vapour_pressure, air_temperature)
    cloud_radiation = σ * cloud_emissivity * (u"K"(air_temperature) - 2.0u"K")^4
    hillshade_radiation = σ * (u"K"(air_temperature))^4 # TODO this is assuming unit emissivity, but should have a vegetation emissivity term
    clear_sky_fraction = 1.0 - cloud_cover
    longwave_radiation_sky = (atmospheric_longwave * clear_sky_fraction + cloud_radiation * cloud_cover) * (1.0 - shade)
    longwave_radiation_vegetation = shade * hillshade_radiation
    incoming_longwave = (longwave_radiation_sky + longwave_radiation_vegetation) * sky_view_fraction +
                        hillshade_radiation * (1.0 - sky_view_fraction)
    outgoing_coeff = (1.0 - shade) * σ * surface_emissivity
    ground_shade_term = shade * hillshade_radiation
    (incoming_longwave, outgoing_coeff, ground_shade_term) = ground_longwave_terms(
        incoming_longwave, outgoing_coeff, ground_shade_term, surface_emissivity, hourly.ground_incoming_longwave)
    Q_infrared = incoming_longwave - outgoing_coeff * (u"K"(soil_temperature[j]))^4 - ground_shade_term

    # Conduction from first active node to the node below it
    # Fortran DSUB.f lines 483-487: on day 1 hour 1, uses T(1) instead of T(j).
    # t == 0 minute corresponds to the very first ODE evaluation of the simulation.
    conduction_ref = (t == 0.0u"minute" && j > 1) ? soil_temperature[1] : soil_temperature[j]
    Q_conduction = j < N ? thermal_conductance[j] * (soil_temperature[j+1] - conduction_ref) : 0.0u"W/m^2"

    # Convection and evaporation using first active node as surface
    # Fortran DSUB.f lines 551-561: temperature safety guards before convection
    surface_temperature = soil_temperature[j]
    if abs(u"°C"(surface_temperature)) > 100u"°C"
        surface_temperature = oftype(surface_temperature, air_temperature + 1.0u"K")  # DSUB.f line 557
    end
    if surface_temperature == air_temperature
        surface_temperature = surface_temperature + 0.1u"K"  # DSUB.f line 561: prevent divide-by-zero
    end
    # Fortran DSUB.f line 565: CALL EVAP(T(1),TAIR,RH,100.0D0,...) — the 100.0 is
    # surface relative humidity (RHSURF), not wetness. PTWET comes from COMMON block.
    convection_ctx = (; boundary_layer_model,
        surface_temperature, air_temperature, wind_speed, relative_humidity,
        atmospheric_pressure, roughness_height, reference_height,
        zenith_angle, soil_wetness, vapour_pressure_equation,
        obukhov_length_prev=buffers.soil_energy_balance.obukhov_length_prev,
    )
    Q_evaporation, convective_heat_flux = ground_heat_vapor_flux(evaporation_model, hourly.ground_air, convection_ctx)

    # Surface energy derivative for the first active node
    # QFREZE: snow melt energy correction, matching Fortran DSUB.f line 567
    surface_dtdt = u"K/minute"((Q_solar + Q_infrared + Q_conduction + convective_heat_flux + Q_freeze - Q_evaporation) / heat_capacity[j])

    # Build derivative vector: inactive nodes follow surface, active nodes use conduction.
    # Fortran DSUB.f:591-606 runs DO 10 I=2,N. T(N) is pinned to TDS at every RHS call
    # (DSUB.f:579), so any DTDT(N) the integrator computes is overwritten next step —
    # effectively a Dirichlet BC. Julia matches that by holding `dt[N] = 0`; T(N) is set
    # once per hour to `environment_instant.deep_soil_temperature` and stays there for
    # the integration interval.
    dt = SVector{N}(ntuple(Val{N}()) do i
        if i <= j
            # Inactive nodes (i < j) follow surface; node j IS the surface
            surface_dtdt
        elseif i == N
            # Lower boundary: pinned at TDS (Dirichlet). Matches Fortran DSUB.f:579.
            zero(surface_dtdt)
        else
            # Internal conduction
            iszero(heat_capacity[i]) ? zero(surface_dtdt) :
                u"K/minute"((thermal_conductance[i-1] * (soil_temperature[i-1] - soil_temperature[i]) + thermal_conductance[i] * (soil_temperature[i+1] - soil_temperature[i])) / heat_capacity[i])
        end
    end)

    return dt
end


