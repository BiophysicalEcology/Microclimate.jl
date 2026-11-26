"""
    SoilHeatTransportModel

Supertype for soil column energy balance models. [`SoilHeatTransport1D`](@ref) is the
only implementation.
"""
abstract type SoilHeatTransportModel end

"""
    soil_energy_balance(model, temperature_state, p, t)

ODE right hand side for the soil temperature column: returns the per-node `dT/dt`
SVector. Variants supply the energy budget formulation (conduction +
surface boundary condition: net radiation + convection + evaporation + phase corrections).
"""
function soil_energy_balance end

"""
    allocate_soil_energy_balance(model, num_nodes)

Allocate the buffers for the energy ODE.
"""
function allocate_soil_energy_balance end

"""
    AbstractGroundTemperatureModel

How the soil gives the ground temperature at the end of an hour for a surface boundary.
Concrete models implement [`ground_temperature`](@ref):
[`IntegratedGroundTemperature`](@ref), [`PredictedGroundTemperature`](@ref) and [`HourStartGroundTemperature`](@ref).
"""
abstract type AbstractGroundTemperatureModel end

"""
    HourStartGroundTemperature()

Keeps the ground at its temperature at the start of the hour. The canopy is solved at that
temperature, and the soil is integrated once, with the canopy's solved boundary.
"""
struct HourStartGroundTemperature <: AbstractGroundTemperatureModel end

"""
    IntegratedGroundTemperature()

Integrates the whole soil column over the hour.
"""
struct IntegratedGroundTemperature <: AbstractGroundTemperatureModel end

"""
    PredictedGroundTemperature(; nodes=3, solver=NewtonRaphson(), tolerance=1e-3u"K")

One backward Euler step over the hour for the surface node and the `nodes - 1` nodes below it.
The deeper nodes keep their temperatures from the start of the hour. `solver`, a NonlinearSolve
algorithm, solves the step to within `tolerance`.
"""
struct PredictedGroundTemperature{N,S,T} <: AbstractGroundTemperatureModel
    solver::S
    tolerance::T
end
function PredictedGroundTemperature(; nodes::Int=3, solver=NewtonRaphson(), tolerance=1e-3u"K")
    nodes >= 1 || throw(ArgumentError("PredictedGroundTemperature needs at least 1 node, got $nodes"))
    return PredictedGroundTemperature{nodes,typeof(solver),typeof(tolerance)}(solver, tolerance)
end

"""
    ground_temperature(model::AbstractGroundTemperatureModel, solver, integrator, hour_state, hour_index; surface_node)

The temperature of node `surface_node` at the end of hour `hour_index`, starting from
`hour_state`. The surface boundary is the one in `integrator`'s inputs. `solver` is
`model`'s own solver, from `allocate_ground_temperature`.
"""
function ground_temperature end

"""
    GroundAirExchange(; air_temperature, air_relative_humidity, wind_speed, height,
        heat_conductance, vapor_conductance)

The air at `height` just above the ground, and how the ground exchanges heat and water
vapor with it. A canopy gives this to the soil, and the soil uses these conductances in place of its own surface exchange.

- Sensible heat flux up from the ground = `heat_conductance * (ground temperature - air_temperature)`.
- Water vapor flux up from the ground = `vapor_conductance * (ground vapor density - air vapor density)`.
"""
@kwdef struct GroundAirExchange{AT,RH,WS,H,HC,VC}
    air_temperature::AT
    air_relative_humidity::RH
    wind_speed::WS
    height::H
    heat_conductance::HC
    vapor_conductance::VC
end

"""
    ReferenceAirExchange(; boundary_layer_model, air_temperature, air_relative_humidity, wind_speed,
        height, roughness_length, zenith_angle, obukhov_length_prev)

The free air at the reference `height`, with no canopy. The ground exchanges heat and water
vapor with it through `boundary_layer_model`.
"""
@kwdef struct ReferenceAirExchange{BLM,AT,RH,WS,H,Z0,ZA,OL}
    boundary_layer_model::BLM
    air_temperature::AT
    air_relative_humidity::RH
    wind_speed::WS
    height::H
    roughness_length::Z0
    zenith_angle::ZA
    obukhov_length_prev::OL
end

"""
    ground_air_transfer(air, surface_temperature, atmospheric_pressure, vapour_pressure_equation)

Heat and water vapor exchange between a surface at `surface_temperature` and the `air`
([`ReferenceAirExchange`](@ref) or [`GroundAirExchange`](@ref)), as
`(; sensible_heat_flux, mass_transfer_coefficient)`.
`sensible_heat_flux` is into the surface.
"""
function ground_air_transfer(air::ReferenceAirExchange, surface_temperature, atmospheric_pressure, vapour_pressure_equation)
    (; boundary_layer_model, air_temperature, air_relative_humidity, wind_speed, height, roughness_length,
       zenith_angle, obukhov_length_prev) = air
    (; convective_heat_flux, temperature_difference) = surface_fluxes(boundary_layer_model;
        surface_temperature, air_temperature, wind_speed, zenith_angle,
        roughness_length, reference_height=height, atmospheric_pressure, obukhov_length_prev,
    )
    heat_transfer_coefficient = calc_heat_transfer_coefficient(convective_heat_flux, temperature_difference)
    wet_air = wet_air_properties(u"K"(air_temperature), air_relative_humidity, atmospheric_pressure; vapour_pressure_equation)
    mass_transfer_coefficient = calc_mass_transfer_coefficient(heat_transfer_coefficient, wet_air.specific_heat, wet_air.density)
    return (; sensible_heat_flux=convective_heat_flux, mass_transfer_coefficient)
end
function ground_air_transfer(air::GroundAirExchange, surface_temperature, atmospheric_pressure, vapour_pressure_equation)
    sensible_heat_flux = uconvert(u"W/m^2", air.heat_conductance * (air.air_temperature - surface_temperature))
    return (; sensible_heat_flux, mass_transfer_coefficient=air.vapor_conductance)
end

"""
    GroundConditions(; direct_shortwave, diffuse_shortwave, sun_elevation_sine, downward_longwave,
        longwave_return_fraction, exposed_fraction, atmospheric_pressure, vapour_pressure_equation, air)

Radiation and air above the ground at one instant.

- `direct_shortwave`, `diffuse_shortwave`: shortwave reaching the ground
- `sun_elevation_sine`: sine of the sun's elevation above the ground
- `downward_longwave`: longwave reaching the ground from above when the ground sends none up
- `longwave_return_fraction`: fraction of the longwave the ground sends up that comes back down
- `exposed_fraction`: fraction of the ground that emits longwave. NicheMapR's `shade` leaves
  the rest with no net longwave.
- `air`: the [`ReferenceAirExchange`](@ref) or [`GroundAirExchange`](@ref)
"""
@kwdef struct GroundConditions{SW,LW,P,VPE,A}
    direct_shortwave::SW
    diffuse_shortwave::SW
    sun_elevation_sine::Float64
    downward_longwave::LW
    longwave_return_fraction::Float64
    exposed_fraction::Float64
    atmospheric_pressure::P
    vapour_pressure_equation::VPE
    air::A
end

"""
    Ground(; soil_albedo, soil_emissivity, soil_longwave_absorption=EmissivityAbsorption(), soil_wetness,
        litter, litter_water_content=0.0, snow_covered)

The ground: the soil surface, and the `litter` ([`AbstractLitterModel`](@ref)) lying on it.
Snow buries the litter while `snow_covered`. Whatever is above the ground uses it only through
[`ground_reflectance`](@ref), [`ground_surface`](@ref) and [`soil_surface_heat_flux`](@ref).

- `soil_albedo`, `soil_emissivity`: of the soil surface, or of snow covering it
- `soil_longwave_absorption`: how much of the longwave reaching it the soil absorbs, an
  [`AbstractGroundLongwaveAbsorption`](@ref)
- `soil_wetness`: the evaporating fraction of the soil surface, which evaporates as if
  saturated at its temperature
- `litter_water_content`: water per dry mass the litter holds, for [`WetLitter`](@ref)
"""
@kwdef struct Ground{A<:AbstractGroundLongwaveAbsorption,L}
    soil_albedo::Float64
    soil_emissivity::Float64
    soil_longwave_absorption::A = EmissivityAbsorption()
    soil_wetness::Float64
    litter::L
    litter_water_content::Float64 = 0.0
    snow_covered::Bool
end

"""
    GroundSurface(; temperature, longwave_reflectance, emitted_longwave, reflectance, vapor_density)

The ground as what is above it sees it.

- `temperature`: of the surface exchanging sensible heat with the air
- `longwave_reflectance`, `emitted_longwave`: the ground reflects `longwave_reflectance` of the
  longwave reaching it, absorbs the rest, and emits `emitted_longwave`
- `reflectance`: shortwave reflectance, one number or a [`GroundShortwaveReflectance`](@ref)
  for the unscattered beam and diffuse light
- `vapor_density`: the water vapor density that, across the air's own resistance, gives the
  ground's evaporation into the air above it. It depends on that air's vapor density.
"""
@kwdef struct GroundSurface{T,LR,LW,R,VD}
    temperature::T
    longwave_reflectance::LR
    emitted_longwave::LW
    reflectance::R
    vapor_density::VD
end

# The surface temperature convection is evaluated at. The reference air's heat transfer coefficient
# divides by the air-surface difference, so it takes DSUB.f's guards against a runaway surface
# (line 557) and a zero difference (line 561); the canopy's ground air exchange is linear in it.
convection_surface_temperature(::GroundAirExchange, surface_temperature) = surface_temperature
function convection_surface_temperature(air::ReferenceAirExchange, surface_temperature)
    (; air_temperature) = air
    abs(u"°C"(surface_temperature)) > 100u"°C" && return oftype(surface_temperature, air_temperature + 1.0u"K")
    surface_temperature == air_temperature && return surface_temperature + 0.1u"K"
    return surface_temperature
end

# Latent heat flux of evaporation from a soil surface at `surface_temperature` (W/m², away
# from the surface) into air at `air_relative_humidity`. DSUB.f line 565 passes 100% surface
# relative humidity.
function ground_latent_heat_flux(evaporation_model, ground::Ground, conditions::GroundConditions, surface_temperature,
    mass_transfer_coefficient; air_relative_humidity=conditions.air.air_relative_humidity,
)
    (; air, atmospheric_pressure, vapour_pressure_equation) = conditions
    latent_heat_flux, _ = evaporation(evaporation_model;
        surface_temperature=u"K"(surface_temperature), air_temperature=u"K"(air.air_temperature),
        relative_humidity=air_relative_humidity, surface_relative_humidity=1.0, mass_transfer_coefficient,
        atmospheric_pressure, soil_wetness=ground.soil_wetness, saturated=false, vapour_pressure_equation,
    )
    return latent_heat_flux
end
