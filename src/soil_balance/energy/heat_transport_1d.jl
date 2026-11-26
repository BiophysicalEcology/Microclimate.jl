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
@kwdef struct HourlySoilEnergyInputs{D<:Vector{<:Number},EI,G<:Ground,LW,QF,SNS,GST,GDL,GA}
    # Node depths, including any active snow nodes
    depths::D
    environment_instant::EI
    ground::G
    longwave_sky::LW
    Q_freeze::QF = 0.0u"W/m^2"
    snow_state::SNS = nothing
    # Canopy ground-boundary overrides (nothing under NoCanopy; see apply_canopy_overrides)
    ground_shortwave_transmission::GST = nothing
    ground_downward_longwave::GDL = nothing
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

# The surface node of the column at `temperatures` and `time`. Snow nodes above the snow
# surface have no heat capacity.
function first_active_node(p::SoilEnergyInputs, temperatures, time)
    (; atmospheric_pressure) = interpolate_forcings(p.forcing, time)
    properties = node_properties(p.model.snow, p, temperatures, atmospheric_pressure, p.model.vapour_pressure_equation)
    return first_active_node(soil_layers(properties, p.hourly.depths, _n_snow(p.model.snow)).heat_capacity)
end

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

# Direct and diffuse shares of `shortwave` reaching the ground. NoCanopy shades it and splits
# by the forcing's diffuse fraction; MultilayerCanopy gives the direct and diffuse fractions
# of the global shortwave that reach the ground (see apply_canopy_overrides).
function ground_shortwave_shares(shortwave, diffuse_fraction, shade, ::Nothing)
    unshaded = shortwave * (1.0 - shade)
    return ((1.0 - diffuse_fraction) * unshaded, diffuse_fraction * unshaded)
end
ground_shortwave_shares(shortwave, diffuse_fraction, shade, transmission::NamedTuple) =
    (shortwave * transmission.direct, shortwave * transmission.diffuse)

# The longwave reaching the ground, and the fraction of the ground that exchanges it. Without a
# canopy, the sky and hillside longwave over the unshaded ground, the shaded ground being in
# longwave balance with its shade (DSUB.f lines 396-465). Under a canopy, the canopy's over all of it.
ground_longwave_terms(open_downward_longwave, shade, ::Nothing) = (open_downward_longwave, 1.0 - shade)
ground_longwave_terms(open_downward_longwave, shade, canopy_downward_longwave) = (canopy_downward_longwave, 1.0)

# The air the ground exchanges with: the reference air without a canopy, the canopy's under one
ground_air_exchange(::Nothing, reference_air) = reference_air
ground_air_exchange(ground_air::GroundAirExchange, reference_air) = ground_air

"""
    ground_conditions(p::SoilEnergyInputs, instant_forcing)

The [`GroundConditions`](@ref) of the hour's soil energy ODE, from its forcing at one instant
and the hour's inputs. Recomputed at every ODE evaluation, as DSUB.f recomputes radiation at
every call.
"""
function ground_conditions(p::SoilEnergyInputs, instant_forcing)
    (; hourly, buffers) = p
    (; site, vapour_pressure_equation) = p.model
    (; heights, boundary_layer_model) = p.model.atmosphere
    (; atmospheric_radiation_model) = p.model.atmosphere.radiation.longwave_model
    (; environment_instant) = hourly
    (; shade, ground_roughness_length, cloud_emissivity) = environment_instant
    (; atmospheric_pressure, air_temperature, wind_speed, zenith_angle, solar_radiation, cloud_cover,
       relative_humidity, slope_zenith_angle, diffuse_fraction) = instant_forcing
    (; slope, sky_view_fraction) = site

    # Solar radiation, split into direct and diffuse so that self-shadowed slopes still receive the diffuse sky
    (direct_shortwave, diffuse_shortwave) = ground_shortwave_shares(solar_radiation, diffuse_fraction, shade,
        hourly.ground_shortwave_transmission)
    sun_elevation_sine = cosd(zenith_angle)
    if slope > 0 && zenith_angle < 90u"°"
        cos_slope_zenith = cosd(slope_zenith_angle)
        direct_shortwave = max(0.0u"W/m^2", (direct_shortwave / cosd(zenith_angle)) * cos_slope_zenith)
        diffuse_shortwave = diffuse_shortwave * sky_view_fraction
        sun_elevation_sine = cos_slope_zenith
    end

    # Longwave from the forcing interpolated to this instant (the interpolated cloud_cover, not
    # the hourly snapshot), or the supplied hourly longwave
    open_downward_longwave = downward_longwave(atmospheric_radiation_model; air_temperature, relative_humidity,
        atmospheric_pressure, cloud_cover, cloud_emissivity, sky_view_fraction,
        supplied_longwave=get(environment_instant, :longwave_radiation, nothing), vapour_pressure_equation).open_ground
    (longwave_reaching_ground, exposed_fraction) = ground_longwave_terms(open_downward_longwave, shade,
        hourly.ground_downward_longwave)

    reference_air = ReferenceAirExchange(; boundary_layer_model, air_temperature, air_relative_humidity=relative_humidity,
        wind_speed, height=last(heights), roughness_length=ground_roughness_length, zenith_angle,
        obukhov_length_prev=buffers.soil_energy_balance.obukhov_length_prev)
    # The downward longwave is held at what the sky, or the canopy's solve, sends down
    return GroundConditions(; direct_shortwave, diffuse_shortwave, sun_elevation_sine, downward_longwave=longwave_reaching_ground,
        longwave_return_fraction=0.0, exposed_fraction, atmospheric_pressure, vapour_pressure_equation,
        air=ground_air_exchange(hourly.ground_air, reference_air))
end

function allocate_soil_energy_balance(::SoilHeatTransport1D, num_nodes::Int)
    obukhov_length_prev = Ref(INITIAL_OBUKHOV_LENGTH)  # warm-start across hourly ODE solves
    return (; obukhov_length_prev)
end

# SciML ODE rhs. Dispatches the actual physics on `p.model.soil.energy_model`.
soil_energy_balance(u, p::SoilEnergyInputs, t) = soil_energy_balance(p.model.soil.energy_model, u, p, t)

function soil_energy_balance(model::SoilHeatTransport1D,
    temperature_state::U,  # state
    p::SoilEnergyInputs,   # "parameters"
    t::Quantity,           # timestep
) where U <: SVector{N} where N
    (; hourly, forcing) = p
    (; vapour_pressure_equation) = p.model
    (; evaporation_model) = p.model.surface
    snow_model = p.model.snow
    (; depths, Q_freeze) = hourly
    maximum_surface_temperature = model.maximum_surface_temperature
    instant_forcing = interpolate_forcings(forcing, t)
    (; atmospheric_pressure) = instant_forcing

    # check for unstable conditions of ground surface temperature.
    # Fortran microinput(74)=maxsurf (MICROCLIMATE.f:482); −81°C lower bound matches
    # DSUB.f:551.
    T_lo = u"K"(-81.0u"°C")
    T_hi = u"K"(maximum_surface_temperature)
    soil_temperature = clamped_temperature = map(t -> clamp(t, T_lo, T_hi), temperature_state)::U

    # Soil and snow properties at the current temperatures, on every evaluation, as the Fortran
    # (DSUB.f lines 274-337) does
    properties = node_properties(snow_model, p, soil_temperature, atmospheric_pressure, vapour_pressure_equation)
    n_snow = _n_snow(snow_model)
    (; heat_capacity, thermal_conductance) = soil_layers(properties, depths, n_snow)
    j = first_active_node(heat_capacity)

    # Radiation, sensible and latent heat into the first active node from above. Snow buries
    # any litter, so a snow surface exchanges with the air directly.
    conditions = ground_conditions(p, instant_forcing)
    (; ground) = hourly
    surface_heat_flux = soil_surface_heat_flux(exposed_litter(ground), ground, conditions, evaporation_model, soil_temperature[j])

    # Conduction from first active node to the node below it
    # Fortran DSUB.f lines 483-487 uses T(1) here on the first call. We always use the
    # surface node, so that only neighbouring nodes enter each derivative.
    Q_conduction = j < N ? thermal_conductance[j] * (soil_temperature[j+1] - soil_temperature[j]) : 0.0u"W/m^2"

    # Surface energy derivative for the first active node
    # QFREZE: snow melt energy correction, matching Fortran DSUB.f line 567
    surface_dtdt = u"K/minute"((surface_heat_flux + Q_conduction + Q_freeze) / heat_capacity[j])

    # Build derivative vector: the surface node takes the surface balance, active nodes use
    # conduction. Inactive snow nodes above the surface hold still, so that each row only
    # depends on its neighbours, as the solver's tridiagonal Jacobian needs. step_to_hour!
    # gives them the surface temperature afterwards. Fortran DSUB.f:591-606 runs DO 10 I=2,N. T(N) is pinned to TDS at every RHS call
    # (DSUB.f:579), so any DTDT(N) the integrator computes is overwritten next step —
    # effectively a Dirichlet BC. Julia matches that by holding `dt[N] = 0`; T(N) is set
    # once per hour to `environment_instant.deep_soil_temperature` and stays there for
    # the integration interval.
    dt = SVector{N}(ntuple(Val{N}()) do i
        if i < j
            zero(surface_dtdt)
        elseif i == j
            surface_dtdt
        elseif i == N
            # Lower boundary: pinned at TDS (Dirichlet). Matches Fortran DSUB.f:579.
            zero(surface_dtdt)
        elseif iszero(heat_capacity[i])
            zero(surface_dtdt)
        else
            # Internal conduction
            u"K/minute"((thermal_conductance[i-1] * (soil_temperature[i-1] - soil_temperature[i]) + thermal_conductance[i] * (soil_temperature[i+1] - soil_temperature[i])) / heat_capacity[i])
        end
    end)

    return dt
end


