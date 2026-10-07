"""
    SoilHeatTransport1D(; ode_solver=Tsit5(),
                          ode_kwargs=(; reltol=1e-6u"K", abstol=1e-3u"K"),
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

- `ode_solver` — any SciML algorithm (`Tsit5()`, `RK4()`, …)
- `ode_kwargs` — NamedTuple of kwargs forwarded to the integrator
- `maximum_surface_temperature` — surface-temperature safety clamp
  (Fortran microinput[74])
- `freezing_model` — soil ice/water phase-transition correction applied to
  the ODE result each hour
"""
@kwdef struct SoilHeatTransport1D{SOS,SOK,MSF,FM} <: SoilHeatTransportModel
    ode_solver::SOS = Tsit5()
    ode_kwargs::SOK = (; reltol=1e-6u"K", abstol=1e-3u"K")
    maximum_surface_temperature::MSF = 85.0u"°C"
    freezing_model::FM = PhaseTransitionLatentHeat()
end

"""
    SoilEnergyInputs

Mutable parameters passed to the ODE integrator via `integrator.p`. Mutable so
the per-hour-changing fields (`environment_instant`, `longwave_sky`, `albedo`,
`Q_freeze`, `snow_state`, `depths`, `soil_wetness`) can be updated in place
instead of allocating a fresh struct each hour. All fields are concretely
typed via the parametric `{F,B,…}` slots; the underlying heap object is the
same across hours.
"""
mutable struct SoilEnergyInputs{M<:SoilHeatTransportModel,EVM<:AbstractEvaporationModel,ARM<:AbstractAtmosphericRadiationModel,F,B,SP,SPR,D<:Vector{<:Number},H<:Vector{<:Number},S,BLM,EI,SW,VP,LW,QF,SNM,SNS,SM,GST,GIL,GWS,GAT,GRH,GRFH,GHC,GVC}
    model::M
    evaporation_model::EVM
    atmospheric_radiation_model::ARM
    forcing::F
    buffers::B
    soil_properties_model::SP
    soil_profile::SPR
    depths::D
    heights::H
    site::S
    boundary_layer_model::BLM
    environment_instant::EI
    soil_wetness::SW
    vapour_pressure_equation::VP
    longwave_sky::LW
    albedo::Float64
    Q_freeze::QF
    # Snow/soil property recomputation inside ODE (matching Fortran DSUB calling SOILPROPS)
    snow_model::SNM
    snow_state::SNS
    soil_moisture::SM
    # Canopy ground-boundary overrides (nothing under NoCanopy; see apply_canopy_overrides)
    ground_shortwave_transmission::GST
    ground_incoming_longwave::GIL
    # Canopy-resolved sub-canopy boundary layer for convection/evaporation
    # (nothing under NoCanopy): the canopy's ground-most-layer wind/air
    # temperature/humidity and that layer's height, used as the reference
    # level in place of the free-atmosphere forcing.
    ground_wind_speed::GWS
    ground_air_temperature::GAT
    ground_air_relative_humidity::GRH
    ground_reference_height::GRFH
    # Air-profile model's own ground-to-lowest-layer conductances (nothing
    # under NoCanopy) -- used instead of surface_fluxes's log-law inversion,
    # which is invalid when ground_reference_height is too close to
    # site.roughness_height (see ground_heat_vapor_flux).
    ground_heat_conductance::GHC
    ground_vapor_conductance::GVC
end

function SoilEnergyInputs(; model=SoilHeatTransport1D(), evaporation_model=BulkTransferEvaporation(),
    atmospheric_radiation_model=CampbellNormanAtmosphericRadiation(),
    forcing, buffers, soil_properties_model, soil_profile, depths, heights, site,
    boundary_layer_model, environment_instant, soil_wetness,
    vapour_pressure_equation=GoffGratch(),
    longwave_sky, albedo, Q_freeze=0.0u"W/m^2",
    snow_model=nothing, snow_state=nothing,
    soil_moisture=nothing,
    ground_shortwave_transmission=nothing, ground_incoming_longwave=nothing,
    ground_wind_speed=nothing, ground_air_temperature=nothing,
    ground_air_relative_humidity=nothing, ground_reference_height=nothing,
    ground_heat_conductance=nothing, ground_vapor_conductance=nothing,
)
    SoilEnergyInputs(model, evaporation_model, atmospheric_radiation_model,
        forcing, buffers, soil_properties_model, soil_profile,
        depths, heights, site, boundary_layer_model, environment_instant,
        soil_wetness, vapour_pressure_equation, longwave_sky, albedo, Q_freeze,
        snow_model, snow_state, soil_moisture,
        ground_shortwave_transmission, ground_incoming_longwave,
        ground_wind_speed, ground_air_temperature, ground_air_relative_humidity, ground_reference_height,
        ground_heat_conductance, ground_vapor_conductance)
end

# In-place per-hour update — no allocation. Only the fields that change per
# hour are written; static fields stay put.
@inline function update_soil_energy_inputs!(p::SoilEnergyInputs;
    depths, environment_instant, soil_wetness, longwave_sky, albedo, Q_freeze,
    snow_state, ground_shortwave_transmission=nothing, ground_incoming_longwave=nothing,
    ground_wind_speed=nothing, ground_air_temperature=nothing,
    ground_air_relative_humidity=nothing, ground_reference_height=nothing,
    ground_heat_conductance=nothing, ground_vapor_conductance=nothing,
)
    p.depths = depths
    p.environment_instant = environment_instant
    p.soil_wetness = soil_wetness
    p.longwave_sky = longwave_sky
    p.albedo = albedo
    p.Q_freeze = Q_freeze
    p.snow_state = snow_state
    p.ground_shortwave_transmission = ground_shortwave_transmission
    p.ground_incoming_longwave = ground_incoming_longwave
    p.ground_wind_speed = ground_wind_speed
    p.ground_air_temperature = ground_air_temperature
    p.ground_air_relative_humidity = ground_air_relative_humidity
    p.ground_reference_height = ground_reference_height
    p.ground_heat_conductance = ground_heat_conductance
    p.ground_vapor_conductance = ground_vapor_conductance
    return p
end

@inline function _first_active_node(heat_capacity, N)
    j = 1
    while j < N && iszero(heat_capacity[j])
        j += 1
    end
    return j
end

# Dispatched soil-properties recompute. Methods for the no-snow case here; the
# SnowModel{N} method lives in snow/snow_model.jl (after SnowModel is defined).
# Both `::Nothing` and `::NoSnow` use the same fall-through; the SnowModel branch
# uses Val(N) so SVector sizes are compile-time.
@inline function _recompute_soil_snow_properties_no_snow!(p, soil_temperature,
    atmospheric_pressure, vapour_pressure_equation)
    soil_properties!(p.buffers.soil_properties, p.soil_properties_model;
        soil_temperature, soil_moisture=p.soil_moisture,
        bulk_density=p.soil_profile.bulk_density, mineral_density=p.soil_profile.mineral_density,
        mineral_conductivity=p.soil_profile.mineral_conductivity, mineral_heat_capacity=p.soil_profile.mineral_heat_capacity,
        atmospheric_pressure, vapour_pressure_equation,
    )
    return nothing
end
@inline _recompute_soil_snow_properties!(::Nothing, p, st, ap, vpe) =
    _recompute_soil_snow_properties_no_snow!(p, st, ap, vpe)

# Compile-time snow node count for an ODE's `p.snow_model`. The `::SnowModel{N}`
# method is in snow/snow_model.jl; both ::Nothing and ::NoSnow return 0.
@inline _n_snow(::Nothing) = 0

# Ground shortwave fraction: NoCanopy uses `shade`; MultilayerCanopy substitutes
# its hourly transmission ratio (see apply_canopy_overrides).
@inline ground_shortwave_fraction(shade, ::Nothing) = 1.0 - shade
@inline ground_shortwave_fraction(shade, transmission) = transmission

# Ground longwave terms: NoCanopy passes the shade-based split through;
# MultilayerCanopy supersedes it with its own converged incoming flux.
@inline ground_longwave_terms(incoming, outgoing_coeff, shade_term, surface_emissivity, ::Nothing) =
    (incoming, outgoing_coeff, shade_term)
@inline ground_longwave_terms(incoming, outgoing_coeff, shade_term, surface_emissivity, override) =
    (override, σ * surface_emissivity, 0.0u"W/m^2")

# Sub-canopy convection/evaporation boundary condition: NoCanopy uses the
# free-atmosphere forcing at reference_height directly; MultilayerCanopy
# substitutes its resolved ground-most-layer wind/air temperature/humidity
# and that layer's height as the reference level. roughness_height is
# untouched in both cases — it's the soil surface's own micro-roughness,
# unaffected by canopy above it.
@inline ground_convection_conditions(air_temperature, wind_speed, relative_humidity, reference_height,
    ::Nothing, ::Nothing, ::Nothing, ::Nothing) =
    (air_temperature, wind_speed, relative_humidity, reference_height)
@inline ground_convection_conditions(air_temperature, wind_speed, relative_humidity, reference_height,
    canopy_air_temperature, canopy_wind_speed, canopy_relative_humidity, canopy_reference_height) =
    (canopy_air_temperature, canopy_wind_speed, canopy_relative_humidity, canopy_reference_height)

# Ground sensible/latent heat exchange: without canopy-supplied conductances,
# use surface_fluxes's Monin-Obukhov log-law inversion (valid for the
# free-atmosphere reference height). MultilayerCanopy instead supplies its
# air-profile model's own ground-to-lowest-layer conductances directly --
# ground_reference_height (the canopy's ground-most layer height) is often
# too close to roughness_height for the log law to invert validly there
# (see _MIN_LOGLAW_RATIO in monin_obukhov.jl).
@inline ground_heat_vapor_flux(evaporation_model, ::Nothing, ::Nothing, ctx) =
    surface_convection_evaporation(evaporation_model; ctx...)
@inline function ground_heat_vapor_flux(evaporation_model, ground_heat_conductance, ground_vapor_conductance, ctx)
    convective_heat_flux = uconvert(u"W/m^2", ground_heat_conductance * (ctx.air_temperature - ctx.surface_temperature))
    Q_evaporation, _ = evaporation(evaporation_model;
        surface_temperature=u"K"(ctx.surface_temperature), air_temperature=u"K"(ctx.air_temperature),
        relative_humidity=ctx.relative_humidity, surface_relative_humidity=1.0,
        mass_transfer_coefficient=ground_vapor_conductance,
        atmospheric_pressure=ctx.atmospheric_pressure, soil_wetness=ctx.soil_wetness,
        saturated=false, vapour_pressure_equation=ctx.vapour_pressure_equation,
    )
    return Q_evaporation, convective_heat_flux
end

function allocate_soil_energy_balance(::SoilHeatTransport1D, num_nodes::Int)
    layer_depths = fill(0.0u"cm", num_nodes + 1)
    heat_capacity = fill(1.0u"J/K/m^2", num_nodes)
    thermal_conductance = fill(1.0u"W/K/m^2", num_nodes)
    obukhov_length_prev = Ref(-0.3u"m")  # warm-start across hourly ODE solves
    return (; layer_depths, heat_capacity, thermal_conductance, obukhov_length_prev)
end

# SciML ODE rhs. Dispatches the actual physics on `p.model`.
soil_energy_balance(u, p::SoilEnergyInputs, t) = soil_energy_balance(p.model, u, p, t)

function soil_energy_balance(model::SoilHeatTransport1D,
    temperature_state::U,  # state
    p::SoilEnergyInputs,   # "parameters"
    t::Quantity,           # timestep
) where U <: SVector{N} where N
    # extract parameters
    (; forcing, buffers, heights, depths, environment_instant, site, boundary_layer_model, soil_wetness, vapour_pressure_equation, longwave_sky, Q_freeze, evaporation_model) = p
    maximum_surface_temperature = model.maximum_surface_temperature
    (; layer_depths, heat_capacity, thermal_conductance) = buffers.soil_energy_balance
    (; shade) = environment_instant
    # Get environmental data at time t
    (; atmospheric_pressure, air_temperature, wind_speed, zenith_angle, solar_radiation, cloud_cover, relative_humidity, slope_zenith_angle, diffuse_fraction) = interpolate_forcings(forcing, t)
    (; roughness_height, slope) = site
    (; karman_constant, dyer_constant) = boundary_layer_model
    albedo = p.albedo

    reference_height = last(heights)
    absorptivity = 1.0 - albedo

    # check for unstable conditions of ground surface temperature.
    # Fortran microinput(74)=maxsurf (MICROCLIMATE.f:482); −81°C lower bound matches
    # DSUB.f:551.
    T_lo = u"K"(-81.0u"°C")
    T_hi = u"K"(maximum_surface_temperature)
    soil_temperature = clamped_temperature = map(t -> clamp(t, T_lo, T_hi), temperature_state)::U

    # Recompute soil/snow properties at current temperatures on every ODE sub-step,
    # matching Fortran DSUB.f lines 274-337 which calls SOILPROPS on every evaluation.
    # Dispatch on snow_model type so N_snow is compile-time (no runtime ntuple).
    _recompute_soil_snow_properties!(p.snow_model, p, soil_temperature,
        atmospheric_pressure, vapour_pressure_equation)
    n_snow = _n_snow(p.snow_model)
    (; bulk_thermal_conductivity, bulk_heat_capacity, bulk_density) = buffers.soil_properties


    # set boundary condition of deep soil temperature
    layer_depths[1:N] = depths
    # Compute soil layer properties (guard zero depths for inactive snow nodes)
    for i in 1:N
        volumetric_heat_capacity = bulk_density[i] * bulk_heat_capacity[i]
        if i == 1
            d2 = layer_depths[2]
            heat_capacity[i] = volumetric_heat_capacity * d2 / 2.0
            thermal_conductance[i] = iszero(d2) ? zero(bulk_thermal_conductivity[1] / oneunit(d2)) : bulk_thermal_conductivity[1] / d2
        else
            Δd = layer_depths[i+1] - layer_depths[i]
            # Fortran DSUB.f special case: when a snow node sits at depth 0 with snow
            # active (depp(i) < 1e-8 and maxsnode1 > 0), the heat capacity uses the
            # full distance to the next node rather than the centered half-cell.
            # This is the first active snow node acting as the surface.
            # Guard: only applies to snow nodes (i <= n_snow). When snow model is active
            # but no snow exists, the first soil node is at depth 0 and must NOT trigger
            # this — its heat capacity uses the normal half-cell formula.
            if i <= n_snow && iszero(layer_depths[i]) && !iszero(layer_depths[i+1])
                heat_capacity[i] = volumetric_heat_capacity * layer_depths[i+1]
            else
                heat_capacity[i] = volumetric_heat_capacity * (layer_depths[i+1] - layer_depths[i-1]) / 2.0
            end
            # Co-located nodes (Δd = 0) occur at the snow-soil interface: the last active
            # snow node and the soil surface node share the same depth (= snow depth).
            # Fortran DSUB.f computes C = k/0 = ∞, which forces the stiff ODE solver to
            # keep the two co-located nodes at the same temperature (perfect thermal coupling).
            # Julia guards against /0 with zero conductance, which completely decouples
            # the snow bottom from the soil surface — breaking the snow insulation effect.
            # Fix: when Δd = 0, use the conductance of the preceding interval so the
            # snow-soil interface maintains thermal coupling (approximates Fortran's ∞).
            prev_Δd = layer_depths[i] - layer_depths[i-1]
            thermal_conductance[i] = if iszero(Δd)
                iszero(prev_Δd) ? zero(bulk_thermal_conductivity[i] / oneunit(Δd)) : bulk_thermal_conductivity[i] / prev_Δd
            else
                bulk_thermal_conductivity[i] / Δd
            end
        end
    end

    # Find first active node (first with nonzero heat capacity).
    # Extracted to a helper so `j` is not both mutated and captured by the
    # ntuple do-block below — that combination boxes `j` and propagates
    # ::Any through every downstream computation.
    j = _first_active_node(heat_capacity, N)

    # Solar radiation. Split into direct + diffuse so that self-shadowed
    # slopes still receive the diffuse sky component. Scaling the entire
    # global flux by `cos_slope_zenith / cos_zenith` (the old behaviour)
    # zeroes both components when the sun is behind the slope, which
    # produced sharp per-pixel sun/shadow streaks at fine DEM resolution.
    Q_solar = absorptivity * solar_radiation * ground_shortwave_fraction(shade, p.ground_shortwave_transmission)
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
    atmospheric_longwave = atmospheric_radiation(p.atmospheric_radiation_model, wet_air_out.vapour_pressure, air_temperature)
    cloud_radiation = σ * cloud_emissivity * (u"K"(air_temperature) - 2.0u"K")^4
    hillshade_radiation = σ * cloud_emissivity * (u"K"(air_temperature))^4
    clear_sky_fraction = 1.0 - cloud_cover
    longwave_radiation_sky = (atmospheric_longwave * clear_sky_fraction + cloud_radiation * cloud_cover) * (1.0 - shade)
    longwave_radiation_vegetation = shade * hillshade_radiation
    incoming_longwave = (longwave_radiation_sky + longwave_radiation_vegetation) * sky_view_fraction +
                        hillshade_radiation * (1.0 - sky_view_fraction)
    outgoing_coeff = (1.0 - shade) * σ * surface_emissivity
    ground_shade_term = shade * hillshade_radiation
    (incoming_longwave, outgoing_coeff, ground_shade_term) = ground_longwave_terms(
        incoming_longwave, outgoing_coeff, ground_shade_term, surface_emissivity, p.ground_incoming_longwave)
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
        surface_temperature = air_temperature + 1.0u"K"  # DSUB.f line 557
    end
    if surface_temperature == air_temperature
        surface_temperature = surface_temperature + 0.1u"K"  # DSUB.f line 561: prevent divide-by-zero
    end
    # Fortran DSUB.f line 565: CALL EVAP(T(1),TAIR,RH,100.0D0,...) — the 100.0 is
    # surface relative humidity (RHSURF), not wetness. PTWET comes from COMMON block.
    # Distinct names from the longwave section above, which must keep using the
    # free-atmosphere air_temperature/relative_humidity, not the canopy-attenuated ones.
    (convection_air_temperature, convection_wind_speed, convection_relative_humidity, convection_reference_height) =
        ground_convection_conditions(air_temperature, wind_speed, relative_humidity, reference_height,
            p.ground_air_temperature, p.ground_wind_speed, p.ground_air_relative_humidity, p.ground_reference_height)
    convection_ctx = (; boundary_layer_model,
        surface_temperature, air_temperature=convection_air_temperature, wind_speed=convection_wind_speed,
        relative_humidity=convection_relative_humidity,
        atmospheric_pressure, roughness_height, reference_height=convection_reference_height,
        zenith_angle, soil_wetness, vapour_pressure_equation,
        obukhov_length_prev=buffers.soil_energy_balance.obukhov_length_prev,
    )
    Q_evaporation, convective_heat_flux = ground_heat_vapor_flux(evaporation_model,
        p.ground_heat_conductance, p.ground_vapor_conductance, convection_ctx)

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
            0.0u"K/minute"
        else
            # Internal conduction
            iszero(heat_capacity[i]) ? 0.0u"K/minute" :
                u"K/minute"((thermal_conductance[i-1] * (soil_temperature[i-1] - soil_temperature[i]) + thermal_conductance[i] * (soil_temperature[i+1] - soil_temperature[i])) / heat_capacity[i])
        end
    end)

    return dt
end


