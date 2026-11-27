"""
    AbstractLitterMoisture

How a litter layer holds water: [`DryLitter`](@ref) or [`WetLitter`](@ref).
"""
abstract type AbstractLitterMoisture end

"""
    DryLitter()

The litter holds no water: soil water vapor passes through it and nothing evaporates from it.
"""
struct DryLitter <: AbstractLitterMoisture end

"""
    WetLitter(; element_vapor_resistance, water_potential_coefficient=-53.72u"m",
        water_potential_exponent=1.32, maximum_relative_humidity=0.999, initial_water_content=0.0)

The litter holds water, as SHAW's residue does (Flerchinger 2017, eqs. 81-82 and "Interception
by Canopy and Residue"). Its gravimetric water content `w`, water per dry mass, sets the water
potential of the litter elements, `ψ = g·a·w^-b` with `a` the `water_potential_coefficient`
and `b` the `water_potential_exponent` (Myrold et al. 1981, for wheat straw), and so their
humidity `exp(M_w·ψ / (R·T))`. Water evaporates from the elements, or condenses onto them,
through `element_vapor_resistance` over the ground the litter covers in series with the air's
resistance, and takes its latent heat from the litter. Over the ground it covers the litter
catches rain until it holds the water in equilibrium with `maximum_relative_humidity`; the
rest falls through to the soil. The soil's and the elements' vapor meet in the litter's air
spaces ([`LitterVapor`](@ref)).

- `element_vapor_resistance`: between the litter elements and the air around them, per
  ground area the litter covers. SHAW suggests 1000-50,000 s/m.
- `water_potential_coefficient`, `water_potential_exponent`: `a` and `b` above
- `maximum_relative_humidity`: element humidity at the water the litter can hold. The default,
  just below saturation, is an estimate, not from SHAW.
- `initial_water_content`: water per dry mass at the start of the simulation

# References
- Flerchinger, G. N. (2017). The Simultaneous Heat and Water (SHAW) Model: Technical
  Documentation, Version 3.0. Technical Report NWRC 2017-02, USDA-ARS Northwest Watershed
  Research Center.
- Myrold, D. D., Elliott, L. F., Papendick, R. I. and Campbell, G. S. (1981). Water potential-
  water content characteristics of wheat straw. Soil Science Society of America Journal 45,
  329-333.
"""
@kwdef struct WetLitter{R,A} <: AbstractLitterMoisture
    element_vapor_resistance::R
    water_potential_coefficient::A = -53.72u"m"
    water_potential_exponent::Float64 = 1.32
    maximum_relative_humidity::Float64 = 0.999
    initial_water_content::Float64 = 0.0
end

initial_litter_water_content(::DryLitter) = 0.0
initial_litter_water_content(moisture::WetLitter) = moisture.initial_water_content

# Volume fraction of the litter layer that is water, at `water_content` water per dry mass
litter_water_fraction(::DryLitter, litter, water_content) = 0.0
litter_water_fraction(::WetLitter, litter, water_content) =
    uconvert(NoUnits, water_content * litter_bulk_density(litter) / water_density)

# Relative humidity at the surface of litter elements holding `water_content` at `temperature`
function litter_element_humidity(moisture::WetLitter, water_content, temperature)
    # Litter holding no water has no humidity at any temperature: its water potential is -Inf. The zero is in the temperature's number type
    iszero(water_content) && return zero(ustrip(u"K", temperature))
    (; water_potential_coefficient, water_potential_exponent) = moisture
    water_potential = g_n * water_potential_coefficient * water_content^(-water_potential_exponent)
    return exp(uconvert(NoUnits, water_molar_mass * water_potential / (R * u"K"(temperature))))
end

# The water content, per dry mass, of litter elements at `maximum_relative_humidity`
function litter_maximum_water_content(moisture::WetLitter, temperature)
    (; water_potential_coefficient, water_potential_exponent, maximum_relative_humidity) = moisture
    water_potential = R * u"K"(temperature) * log(maximum_relative_humidity) / water_molar_mass
    return uconvert(NoUnits, water_potential / (g_n * water_potential_coefficient))^(-1 / water_potential_exponent)
end

"""
    AbstractLitterConductance

Heat and water vapor transfer through a litter layer: [`ShawLitterConductance`](@ref).
"""
abstract type AbstractLitterConductance end

"""
    ShawLitterConductance(; still_air_conductivity=0.025u"W/m/K", temperature_coefficient=0.007u"K^-1",
        wind_coefficient=4.0u"s/m", material_conductivity=0.05u"W/m/K", water_conductivity=0.57u"W/m/K",
        material_density=170.0u"kg/m^3")

SHAW's litter conductivity (Flerchinger 2000, eqs. 57-58): convection through the air in the
litter, enhanced by wind within it, plus conduction through the litter material and water,
each weighted by its volume fraction:

`k = k_air(1 + temperature_coefficient·T)(1 + wind_coefficient·u)·f_air + material_conductivity·f_material + water_conductivity·f_water`

with `T` in °C and `f_material` the litter's bulk density over `material_density`. The
convective part over the air's volumetric heat capacity is the water vapor diffusivity.
`wind_coefficient` is 4.0 s/m for wheat straw and 8.5 s/m for larger elements such as corn
stalks. The defaults are SHAW 3.0.3's constants.
"""
@kwdef struct ShawLitterConductance{KA,TC,WC,KM,KW,DM} <: AbstractLitterConductance
    still_air_conductivity::KA = 0.025u"W/m/K"
    temperature_coefficient::TC = 0.007u"K^-1"
    wind_coefficient::WC = 4.0u"s/m"
    material_conductivity::KM = 0.05u"W/m/K"
    water_conductivity::KW = 0.57u"W/m/K"
    material_density::DM = 170.0u"kg/m^3"
end

"""
    AbstractLitterWind

Wind within a litter layer: [`ShawLitterWind`](@ref).
"""
abstract type AbstractLitterWind end

"""
    ShawLitterWind(; attenuation_per_density=0.1u"m^3/kg", sparse_density=20.0u"kg/m^3")

SHAW's wind within litter: it falls exponentially with depth below the litter top at a rate
`attenuation_per_density·(ρ - sparse_density)` per layer depth for litter bulk density `ρ`
(none for sparser litter). It is taken at a quarter of the layer's depth, as SHAW does for a
single litter node.
"""
@kwdef struct ShawLitterWind{A,D} <: AbstractLitterWind
    attenuation_per_density::A = 0.1u"m^3/kg"
    sparse_density::D = 20.0u"kg/m^3"
end

"""
    AbstractLitterRadiation

Radiation transfer through a litter layer: [`GapLitterRadiation`](@ref) or
[`ShawLitterRadiation`](@ref). Both solve the layer as two identical sublayers over the soil,
as SHAW does for a single litter layer. Litter elements lie flat, so they scatter only
backwards.
"""
abstract type AbstractLitterRadiation end

"""
    GapLitterRadiation()

Gaps between litter elements pass all light, direct and diffuse, so the soil between them is
lit as bare ground: the layer transmits `1 - cover_fraction` of each.
"""
struct GapLitterRadiation <: AbstractLitterRadiation end

"""
    ShawLitterRadiation(; diffuse_transmission_factor=0.667)

SHAW's litter radiation. Gaps between litter elements pass `1 - cover_fraction` of light from
overhead, and less of slanting light, which crosses more of the layer: the direct beam passes
`(1 - cover_fraction)^(1/sin β)` at sun elevation `β` (Flerchinger 2000, eq. 18), and diffuse light
`(1 - cover_fraction)^(1/diffuse_transmission_factor)`, as if it all came from an elevation with
that sine. Zero cover passes all light.
"""
@kwdef struct ShawLitterRadiation{F} <: AbstractLitterRadiation
    diffuse_transmission_factor::F = 0.667
end

"""
    LitterLayer(; dry_mass, thickness, cover_fraction, albedo=0.25, emissivity=0.95,
        moisture=DryLitter(), conductance=ShawLitterConductance(), wind=ShawLitterWind(),
        radiation=GapLitterRadiation(), solver=ITP(), tolerance=1e-4u"K",
        temperature_bounds=(-81.0u"°C", 100.0u"°C"))

A layer of plant litter over part of the soil. Where it covers the ground, the litter absorbs
and emits radiation, exchanges sensible heat with the air at its top, and passes heat and soil
water vapor between its top and the soil by conduction and wind-driven convection. The soil
between litter elements exchanges heat and water vapor with the air directly, so zero cover is
bare ground. The litter stores no heat, so its temperature balances its fluxes at every instant.
The `Shaw` components follow SHAW's formulations (Flerchinger 2000).

- `dry_mass`: dry litter per ground area
- `thickness`: depth of the layer
- `cover_fraction`: fraction of the ground the litter covers
- `albedo`, `emissivity`: of the litter elements

The formulations are components: `moisture` ([`AbstractLitterMoisture`](@ref)), `conductance`
([`AbstractLitterConductance`](@ref)), `wind` ([`AbstractLitterWind`](@ref)) and `radiation`
([`AbstractLitterRadiation`](@ref)). `solver`, a bracketing NonlinearSolve algorithm, solves the
litter temperature to within `tolerance`, between `temperature_bounds`: by default the soil's
-81 °C floor and the 100 °C ceiling NicheMapR's DSUB.f puts on the surface.
Snow buries the litter: while snow lies, the snow surface exchanges with the air directly.

# References
- Flerchinger, G. N. (2000). The Simultaneous Heat and Water (SHAW) Model: Technical
  Documentation. Technical Report NWRC 2000-09, USDA-ARS Northwest Watershed Research Center.
"""
@kwdef struct LitterLayer{M,T,MO<:AbstractLitterMoisture,
    C<:AbstractLitterConductance,W<:AbstractLitterWind,R<:AbstractLitterRadiation,S,TO,TB} <: AbstractLitterModel
    dry_mass::M
    thickness::T
    cover_fraction::Float64
    albedo::Float64 = 0.25
    emissivity::Float64 = 0.95
    moisture::MO = DryLitter()
    conductance::C = ShawLitterConductance()
    wind::W = ShawLitterWind()
    radiation::R = GapLitterRadiation()
    solver::S = ITP()
    tolerance::TO = 1e-4u"K"
    temperature_bounds::TB = (-81.0u"°C", 100.0u"°C")
end

litter_bulk_density(litter::LitterLayer) = litter.dry_mass / litter.thickness

initial_litter_water_content(litter::LitterLayer) = initial_litter_water_content(litter.moisture)

"""
    litter_top_wind_speed(litter::LitterLayer, air)

Wind speed at the top of the litter. Under a canopy, the wind of the canopy's lowest layer, as
SHAW takes the wind at the canopy's lowest node. Without one, the neutral log law from the
reference height down to the litter top.
"""
litter_top_wind_speed(litter::LitterLayer, air::GroundAirExchange) = air.wind_speed
function litter_top_wind_speed(litter::LitterLayer, air::ReferenceAirExchange)
    (; wind_speed, roughness_length, height) = air
    return wind_speed * log((litter.thickness + roughness_length) / roughness_length) /
        log((height + roughness_length) / roughness_length)
end

"""
    litter_wind_speed(wind::AbstractLitterWind, litter, top_wind_speed)

Wind speed within the litter, which drives convection through it.
"""
function litter_wind_speed(wind::ShawLitterWind, litter, top_wind_speed)
    attenuation = max(uconvert(NoUnits, wind.attenuation_per_density * (litter_bulk_density(litter) - wind.sparse_density)), 0.0)
    return top_wind_speed * exp(-attenuation / 4)
end
litter_wind_speed(litter::LitterLayer, air) = litter_wind_speed(litter.wind, litter, litter_top_wind_speed(litter, air))

"""
    litter_conductivity(conductance::AbstractLitterConductance, litter, water_content, wind_speed, temperature)

Thermal conductivity of the litter layer holding `water_content` per dry mass at `temperature`
with `wind_speed` within it, as `(; heat, convective)`: the total and the part carried by air
convection.
"""
function litter_conductivity(conductance::ShawLitterConductance, litter, water_content, wind_speed, temperature)
    (; still_air_conductivity, temperature_coefficient, wind_coefficient, material_conductivity,
       water_conductivity, material_density) = conductance
    material_fraction = uconvert(NoUnits, litter_bulk_density(litter) / material_density)
    material_fraction <= 1 || throw(ArgumentError("litter bulk density $(litter_bulk_density(litter)) exceeds its material density $material_density"))
    water_fraction = litter_water_fraction(litter.moisture, litter, water_content)
    # Water beyond the pore space leaves no air in it
    air_fraction = max(1.0 - material_fraction - water_fraction, 0.0)
    temperature_factor = 1.0 + temperature_coefficient * (u"K"(temperature) - 273.15u"K")
    convective = still_air_conductivity * (1.0 + wind_coefficient * wind_speed) * air_fraction * temperature_factor
    conductive = material_conductivity * material_fraction + water_conductivity * water_fraction
    return (; heat=convective + conductive, convective)
end

"""
    litter_transfer(litter::LitterLayer, water_content, wind_speed, temperature)

Heat conductance and water vapor resistance across the litter layer holding `water_content`
per dry mass at `temperature` with `wind_speed` within it, as
`(; heat_conductance, vapor_resistance)`.
"""
function litter_transfer(litter::LitterLayer, water_content, wind_speed, temperature)
    conductivity = litter_conductivity(litter.conductance, litter, water_content, wind_speed, temperature)
    vapor_diffusivity = conductivity.convective / calc_air_volumetric_heat_capacity(u"K"(temperature))
    return (; heat_conductance=uconvert(u"W/m^2/K", conductivity.heat / litter.thickness),
        vapor_resistance=uconvert(u"s/m", litter.thickness / vapor_diffusivity))
end

"""
    litter_gap_fraction(litter::LitterLayer)

Fraction of the ground between litter elements, open to the sky and the air. SHAW caps the
cover's extinction at exp(-6.9).
"""
litter_gap_fraction(litter::LitterLayer) = max(1.0 - litter.cover_fraction, exp(-6.9))

"""
    litter_sublayer_optics(radiation::AbstractLitterRadiation, litter, sun_elevation_sine)

Direct and diffuse transmission of each of the litter's two radiation sublayers, as
`(; direct_transmission, diffuse_transmission)`.
"""
function litter_sublayer_optics(::GapLitterRadiation, litter, sun_elevation_sine)
    transmission = sqrt(litter_gap_fraction(litter))
    return (; direct_transmission=transmission, diffuse_transmission=transmission)
end
function litter_sublayer_optics(radiation::ShawLitterRadiation, litter, sun_elevation_sine)
    gap_fraction = litter_gap_fraction(litter)
    # Each sublayer has half the path. No beam below the horizon.
    direct_transmission = sun_elevation_sine > 0 ? gap_fraction^(1 / (2 * sun_elevation_sine)) : 0.0
    diffuse_transmission = gap_fraction^(1 / (2 * radiation.diffuse_transmission_factor))
    return (; direct_transmission, diffuse_transmission)
end

# Diffuse two-stream fluxes through two identical sublayers over a reflecting bottom; returns each interface's fluxes and the top's reflectance and source
function litter_two_stream(reflectance, transmittance, upward_sources, downward_sources,
    bottom_reflectance, bottom_source, top_downward,
)
    # The upward flux at each interface is its reflectance times the downward flux there plus its source
    middle_reflectance = reflectance + transmittance^2 * bottom_reflectance / (1.0 - reflectance * bottom_reflectance)
    middle_source = transmittance * (bottom_source + bottom_reflectance * downward_sources[2]) /
        (1.0 - reflectance * bottom_reflectance) + upward_sources[2]
    top_reflectance = reflectance + transmittance^2 * middle_reflectance / (1.0 - reflectance * middle_reflectance)
    top_source = transmittance * (middle_source + middle_reflectance * downward_sources[1]) /
        (1.0 - reflectance * middle_reflectance) + upward_sources[1]
    middle_downward = (transmittance * top_downward + reflectance * middle_source + downward_sources[1]) /
        (1.0 - reflectance * middle_reflectance)
    bottom_downward = (transmittance * middle_downward + reflectance * bottom_source + downward_sources[2]) /
        (1.0 - reflectance * bottom_reflectance)
    return (; downward=(top_downward, middle_downward, bottom_downward),
        upward=(top_reflectance * top_downward + top_source, middle_reflectance * middle_downward + middle_source,
            bottom_reflectance * bottom_downward + bottom_source),
        top_reflectance, top_source)
end

"""
    litter_shortwave(litter::LitterLayer, direct_shortwave, diffuse_shortwave, sun_elevation_sine, soil_albedo)

Shortwave absorbed by the litter and by the soil below it, and reflected to the air, as
`(; litter_absorbed, soil_absorbed, reflected)`.
"""
function litter_shortwave(litter::LitterLayer, direct_shortwave, diffuse_shortwave, sun_elevation_sine, soil_albedo)
    (; direct_transmission, diffuse_transmission) = litter_sublayer_optics(litter.radiation, litter, sun_elevation_sine)
    (; albedo) = litter
    direct = (direct_shortwave, direct_shortwave * direct_transmission, direct_shortwave * direct_transmission^2)
    # Litter elements backscatter what they intercept of the beam as diffuse light
    upward_sources = (albedo * (1.0 - direct_transmission) * direct[1], albedo * (1.0 - direct_transmission) * direct[2])
    no_source = zero(direct_shortwave)
    (; downward, upward) = litter_two_stream(albedo * (1.0 - diffuse_transmission), diffuse_transmission,
        upward_sources, (no_source, no_source), soil_albedo, soil_albedo * direct[3], diffuse_shortwave)
    litter_absorbed = (1.0 - albedo) * ((1.0 - direct_transmission) * (direct[1] + direct[2]) +
        (1.0 - diffuse_transmission) * (downward[1] + upward[2] + downward[2] + upward[3]))
    soil_absorbed = (1.0 - soil_albedo) * (downward[3] + direct[3])
    return (; litter_absorbed, soil_absorbed, reflected=upward[1])
end

# Longwave through the litter's two sublayers with `top_downward` entering from above. The
# soil absorbs the longwave reaching it as its `soil_longwave_absorption` says and reflects the rest.
function litter_longwave_streams(litter::LitterLayer, ground::Ground, top_downward,
    litter_temperature, soil_surface_temperature,
)
    (; diffuse_transmission) = litter_sublayer_optics(litter.radiation, litter, 0.0)
    (; emissivity) = litter
    # Emitted from each face of each sublayer, through the fraction of it the litter fills
    emission = (1.0 - diffuse_transmission) * emissivity * σ * u"K"(litter_temperature)^4
    soil_emission = ground.soil_emissivity * σ * u"K"(soil_surface_temperature)^4
    streams = litter_two_stream((1.0 - emissivity) * (1.0 - diffuse_transmission), diffuse_transmission,
        (emission, emission), (emission, emission), soil_longwave_reflectance(ground), soil_emission, top_downward)
    return (; streams..., emission, soil_emission, diffuse_transmission)
end

"""
    litter_longwave(litter::LitterLayer, ground::Ground, conditions::GroundConditions, litter_temperature, soil_surface_temperature)

Net longwave into the litter and into the soil below it, and the longwave leaving the ground
upward, as `(; litter_net, soil_net, upward)`. The longwave reaching the ground includes the
`longwave_return_fraction` of what leaves it. Only the `exposed_fraction` of the ground has a net
exchange, as for bare soil.
"""
function litter_longwave(litter::LitterLayer, ground::Ground, conditions::GroundConditions, litter_temperature,
    soil_surface_temperature,
)
    (; downward_longwave, longwave_return_fraction, exposed_fraction) = conditions
    # The upward longwave is linear in the downward, so the returned share solves directly
    unlit = litter_longwave_streams(litter, ground, zero(downward_longwave), litter_temperature,
        soil_surface_temperature)
    top_downward = (downward_longwave + longwave_return_fraction * unlit.top_source) /
        (1.0 - longwave_return_fraction * unlit.top_reflectance)
    (; downward, upward, emission, soil_emission, diffuse_transmission) = litter_longwave_streams(litter, ground, top_downward,
        litter_temperature, soil_surface_temperature)
    litter_absorbed = litter.emissivity * (1.0 - diffuse_transmission) * (downward[1] + upward[2] + downward[2] + upward[3])
    soil_absorbed = absorbed_ground_longwave(ground.soil_longwave_absorption, ground.soil_emissivity, downward[3])
    return (; litter_net=exposed_fraction * (litter_absorbed - 4 * emission),
        soil_net=exposed_fraction * (soil_absorbed - soil_emission), upward=upward[1])
end

"""
    LitterVapor(; covered_fraction, half_resistance, element_conductance, element_vapor_density, soil_wetness)

The litter's water vapor network at one instant, as SHAW's residue has it (Flerchinger 2017,
eq. 81). Over the ground the litter covers, the air in its spaces exchanges vapor with the
soil below across half the layer, `half_resistance`, with the air above across the other half
and the air's own resistance, and with the litter elements through `element_conductance`
from their `element_vapor_density`. The soil between the elements exchanges with the air
directly. Its vapor storage is neglected. Only the `soil_wetness` fraction of the soil
evaporates, as if saturated at its temperature.
"""
@kwdef struct LitterVapor{R,G,D}
    covered_fraction::Float64
    half_resistance::R
    element_conductance::G
    element_vapor_density::D
    soil_wetness::Float64
end

# No litter: the soil exchanges with the air directly
no_litter_vapor() = LitterVapor(; covered_fraction=0.0, half_resistance=0.0u"s/m", element_conductance=0.0u"m/s",
    element_vapor_density=0.0u"kg/m^3", soil_wetness=0.0)

litter_element_conductance(::DryLitter) = 0.0u"m/s"
litter_element_conductance(moisture::WetLitter) = 1 / moisture.element_vapor_resistance

# Water vapor density at the surface of the litter elements at `temperature`
litter_element_vapor_density(::DryLitter, ground::Ground, conditions::GroundConditions, temperature) = 0.0u"kg/m^3"
function litter_element_vapor_density(moisture::WetLitter, ground::Ground, conditions::GroundConditions, temperature)
    (; atmospheric_pressure, vapour_pressure_equation) = conditions
    humidity = litter_element_humidity(moisture, ground.litter_water_content, temperature)
    return wet_air_properties(u"K"(temperature), humidity, atmospheric_pressure; vapour_pressure_equation).vapour_density
end

# The litter's vapor network with the litter at `temperature` and `transfer` across it
litter_vapor(litter::LitterLayer, ground::Ground, conditions::GroundConditions, temperature, transfer) =
    LitterVapor(; covered_fraction=1.0 - litter_gap_fraction(litter), half_resistance=transfer.vapor_resistance / 2,
        element_conductance=litter_element_conductance(litter.moisture),
        element_vapor_density=litter_element_vapor_density(litter.moisture, ground, conditions, temperature),
        soil_wetness=ground.soil_wetness)

# Conductances from the air in the litter's spaces to the soil's wet fraction, the air above,
# with the air's own `air_mass_transfer_coefficient`, and the litter elements, per ground area covered
function pore_conductances(vapor::LitterVapor, air_mass_transfer_coefficient)
    (; half_resistance, element_conductance, soil_wetness) = vapor
    return (; soil=soil_wetness / half_resistance, air=1 / (half_resistance + 1 / air_mass_transfer_coefficient),
        elements=element_conductance)
end

# Vapor density of the air in the litter's spaces, where the soil's, the elements' and the air
# above's fluxes into it balance
function pore_vapor_density(vapor::LitterVapor, air_mass_transfer_coefficient, soil_vapor_density, air_vapor_density)
    (; soil, air, elements) = pore_conductances(vapor, air_mass_transfer_coefficient)
    return (soil * soil_vapor_density + elements * vapor.element_vapor_density + air * air_vapor_density) /
        (soil + elements + air)
end

"""
    soil_vapor_exchange(vapor::LitterVapor, air_mass_transfer_coefficient, air_temperature, air_relative_humidity,
        atmospheric_pressure, vapour_pressure_equation)

What the soil surface evaporates into through the litter, over all the ground, as
`(; mass_transfer_coefficient, relative_humidity)`: an equivalent air, at `air_temperature`,
behind one conductance, which the soil's wet fraction evaporates into. Under the litter the
soil's vapor reaches its spaces, which the elements and the air above also feed; between its
elements the soil meets the air directly.
"""
function soil_vapor_exchange(vapor::LitterVapor, air_mass_transfer_coefficient, air_temperature, air_relative_humidity,
    atmospheric_pressure, vapour_pressure_equation,
)
    if iszero(vapor.covered_fraction)
        return (; mass_transfer_coefficient=air_mass_transfer_coefficient, relative_humidity=air_relative_humidity)
    end
    (; soil, air, elements) = pore_conductances(vapor, air_mass_transfer_coefficient)
    saturated_vapor_density = wet_air_properties(u"K"(air_temperature), 1.0, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    air_vapor_density = air_relative_humidity * saturated_vapor_density
    # Per unit of wet fraction, as the soil's evaporation multiplies by its wetness
    covered_conductance = (elements + air) / (vapor.half_resistance * (soil + elements + air))
    covered_vapor_density = (elements * vapor.element_vapor_density + air * air_vapor_density) / (elements + air)
    covered = vapor.covered_fraction
    gap = 1.0 - covered
    mass_transfer_coefficient = covered * covered_conductance + gap * air_mass_transfer_coefficient
    vapor_density = (covered * covered_conductance * covered_vapor_density + gap * air_mass_transfer_coefficient * air_vapor_density) /
        mass_transfer_coefficient
    return (; mass_transfer_coefficient, relative_humidity=uconvert(NoUnits, vapor_density / saturated_vapor_density))
end

"""
    litter_element_evaporation(vapor::LitterVapor, air_mass_transfer_coefficient, soil_vapor_density, air_vapor_density)

Water evaporating from the litter elements into the air in the litter's spaces, per ground
area: negative when vapor condenses onto them.
"""
function litter_element_evaporation(vapor::LitterVapor, air_mass_transfer_coefficient, soil_vapor_density, air_vapor_density)
    pore = pore_vapor_density(vapor, air_mass_transfer_coefficient, soil_vapor_density, air_vapor_density)
    return uconvert(u"kg/m^2/s", vapor.covered_fraction * vapor.element_conductance * (vapor.element_vapor_density - pore))
end

# Vapor densities at the soil surface's wet fraction, saturated, and of the air above
function soil_and_air_vapor_density(conditions::GroundConditions, soil_surface_temperature)
    (; atmospheric_pressure, vapour_pressure_equation, air) = conditions
    soil = wet_air_properties(u"K"(soil_surface_temperature), 1.0, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    air_vapor_density = wet_air_properties(u"K"(air.air_temperature), air.air_relative_humidity, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    return (; soil, air=air_vapor_density)
end

# Latent heat of the water evaporating from the litter elements at `temperature`, per ground area
litter_latent_heat_flux(::DryLitter, litter, ground::Ground, conditions::GroundConditions, temperature, transfer,
    soil_surface_temperature, air_mass_transfer_coefficient) = 0.0u"W/m^2"
function litter_latent_heat_flux(::WetLitter, litter, ground::Ground, conditions::GroundConditions, temperature, transfer,
    soil_surface_temperature, air_mass_transfer_coefficient,
)
    vapor = litter_vapor(litter, ground, conditions, temperature, transfer)
    densities = soil_and_air_vapor_density(conditions, soil_surface_temperature)
    evaporation =litter_element_evaporation(vapor, air_mass_transfer_coefficient, densities.soil, densities.air)
    return uconvert(u"W/m^2", evaporation * enthalpy_of_vaporisation(u"K"(temperature)))
end

# Net heat flux into the litter at `temperature`: absorbed shortwave, net longwave, sensible
# heat from the air and from the soil over the ground it covers, less the latent heat of the
# water evaporating from it
function litter_net_heat_flux(litter::LitterLayer, ground::Ground, conditions::GroundConditions, absorbed_shortwave,
    wind_speed, temperature, soil_surface_temperature,
)
    (; atmospheric_pressure, vapour_pressure_equation, air) = conditions
    longwave = litter_longwave(litter, ground, conditions, temperature, soil_surface_temperature)
    exchange_temperature = convection_surface_temperature(air, temperature)
    (; sensible_heat_flux, mass_transfer_coefficient) = ground_air_transfer(air, exchange_temperature,
        atmospheric_pressure, vapour_pressure_equation)
    transfer = litter_transfer(litter, ground.litter_water_content, wind_speed, temperature)
    latent_heat_flux = litter_latent_heat_flux(litter.moisture, litter, ground, conditions, temperature, transfer,
        soil_surface_temperature, mass_transfer_coefficient)
    covered_fraction = 1.0 - litter_gap_fraction(litter)
    return absorbed_shortwave + longwave.litter_net +
        covered_fraction * (sensible_heat_flux + transfer.heat_conductance * (soil_surface_temperature - temperature)) -
        latent_heat_flux
end

# The litter's net heat flux in W/m² at a litter temperature in kelvin, with the soil surface temperature in kelvin as the problem's parameter
struct LitterHeatBalance{L,G,C,A,W}
    litter::L
    ground::G
    conditions::C
    absorbed_shortwave::A
    wind_speed::W
end
function (balance::LitterHeatBalance)(temperature_in_kelvin, soil_surface_temperature_in_kelvin)
    (; litter, ground, conditions, absorbed_shortwave, wind_speed) = balance
    flux = litter_net_heat_flux(litter, ground, conditions, absorbed_shortwave, wind_speed,
        temperature_in_kelvin * u"K", soil_surface_temperature_in_kelvin * u"K")
    return ustrip(u"W/m^2", flux)
end

"""
    litter_temperature(litter::LitterLayer, ground::Ground, conditions::GroundConditions, absorbed_shortwave,
        soil_surface_temperature)

Temperature of the litter, as `(; temperature, transfer)` with the [`litter_transfer`](@ref)
across it, from its energy balance with `absorbed_shortwave` over a soil surface at
`soil_surface_temperature`.
"""
function litter_temperature(litter::LitterLayer, ground::Ground, conditions::GroundConditions, absorbed_shortwave,
    soil_surface_temperature,
)
    wind_speed = litter_wind_speed(litter, conditions.air)
    balance = LitterHeatBalance(litter, ground, conditions, absorbed_shortwave, wind_speed)
    bounds = map(bound -> ustrip(u"K", bound), litter.temperature_bounds)
    problem = IntervalNonlinearProblem{false}(balance, bounds, ustrip(u"K", soil_surface_temperature))
    solution = solve(problem, litter.solver; abstol=ustrip(u"K", litter.tolerance))
    if !SciMLBase.successful_retcode(solution)
        error("Litter temperature solve failed with $(solution.retcode) within $(litter.temperature_bounds)")
    end
    temperature = solution.u * u"K"
    return (; temperature, transfer=litter_transfer(litter, ground.litter_water_content, wind_speed, temperature))
end

# The shortwave, litter temperature and transfer across it, its vapor network, and the soil's
# own exchange with the air between litter elements, over a soil surface at `soil_surface_temperature`
function litter_state(litter::LitterLayer, ground::Ground, conditions::GroundConditions, soil_surface_temperature)
    (; direct_shortwave, diffuse_shortwave, sun_elevation_sine, atmospheric_pressure, vapour_pressure_equation, air) = conditions
    shortwave = litter_shortwave(litter, direct_shortwave, diffuse_shortwave, sun_elevation_sine, ground.soil_albedo)
    (; temperature, transfer) = litter_temperature(litter, ground, conditions, shortwave.litter_absorbed, soil_surface_temperature)
    exchange_temperature = convection_surface_temperature(air, soil_surface_temperature)
    soil_air = ground_air_transfer(air, exchange_temperature, atmospheric_pressure, vapour_pressure_equation)
    vapor = litter_vapor(litter, ground, conditions, temperature, transfer)
    return (; shortwave, temperature, transfer, exchange_temperature, soil_air, vapor)
end

function soil_surface_heat_flux(litter::LitterLayer, ground::Ground, conditions::GroundConditions, evaporation_model,
    soil_surface_temperature,
)
    (; atmospheric_pressure, vapour_pressure_equation, air) = conditions
    (; shortwave, temperature, transfer, exchange_temperature, soil_air, vapor) =
        litter_state(litter, ground, conditions, soil_surface_temperature)
    longwave = litter_longwave(litter, ground, conditions, temperature, soil_surface_temperature)
    gap_fraction = litter_gap_fraction(litter)
    (; mass_transfer_coefficient, relative_humidity) = soil_vapor_exchange(vapor, soil_air.mass_transfer_coefficient,
        air.air_temperature, air.air_relative_humidity, atmospheric_pressure, vapour_pressure_equation)
    latent_heat_flux = ground_latent_heat_flux(evaporation_model, ground, conditions, exchange_temperature, mass_transfer_coefficient;
        air_relative_humidity=relative_humidity)
    # Under the litter the soil exchanges heat with the litter; between its elements, with the air
    return shortwave.soil_absorbed + longwave.soil_net +
        (1.0 - gap_fraction) * transfer.heat_conductance * (temperature - soil_surface_temperature) +
        gap_fraction * soil_air.sensible_heat_flux - latent_heat_flux
end

function ground_reflectance(litter::LitterLayer, ground::Ground, sun_elevation_sine)
    # The share of a unit of direct, or of diffuse, shortwave that the litter and soil reflect
    irradiance = 1.0u"W/m^2"
    reflected(direct, diffuse) = litter_shortwave(litter, direct, diffuse, sun_elevation_sine, ground.soil_albedo).reflected / irradiance
    return GroundShortwaveReflectance(; direct=reflected(irradiance, zero(irradiance)), diffuse=reflected(zero(irradiance), irradiance))
end

function ground_surface(litter::LitterLayer, ground::Ground, conditions::GroundConditions, soil_surface_temperature)
    (; downward_longwave, sun_elevation_sine) = conditions
    (; temperature, soil_air, vapor) = litter_state(litter, ground, conditions, soil_surface_temperature)
    # The ground's upward longwave is its emission plus its reflectance times the downward
    unlit = litter_longwave_streams(litter, ground, zero(downward_longwave), temperature,
        soil_surface_temperature)
    densities = soil_and_air_vapor_density(conditions, soil_surface_temperature)
    # The air sees litter over the ground it covers and soil between its elements
    gap_fraction = litter_gap_fraction(litter)
    return GroundSurface(; temperature=(1.0 - gap_fraction) * temperature + gap_fraction * u"K"(soil_surface_temperature),
        longwave_reflectance=unlit.top_reflectance, emitted_longwave=unlit.top_source,
        reflectance=ground_reflectance(litter, ground, sun_elevation_sine),
        vapor_density=ground_vapor_source(vapor, soil_air.mass_transfer_coefficient, densities.soil, densities.air))
end

# The ground's water vapor source, its soil's wet fraction at `soil_vapor_density` and its litter
# elements, as the vapor density that gives the same flux into air at `air_vapor_density` across
# the air's own `air_mass_transfer_coefficient`
function ground_vapor_source(vapor::LitterVapor, air_mass_transfer_coefficient, soil_vapor_density, air_vapor_density)
    (; air) = pore_conductances(vapor, air_mass_transfer_coefficient)
    pore = pore_vapor_density(vapor, air_mass_transfer_coefficient, soil_vapor_density, air_vapor_density)
    covered = vapor.covered_fraction
    gap = 1.0 - covered
    to_air = covered * air * (pore - air_vapor_density) +
        gap * vapor.soil_wetness * air_mass_transfer_coefficient * (soil_vapor_density - air_vapor_density)
    return air_vapor_density + to_air / air_mass_transfer_coefficient
end

"""
    litter_water_exchange(ground::Ground, conditions::GroundConditions, soil_surface_temperature)

Water evaporating from the exposed litter's elements per ground area over a soil surface at
`soil_surface_temperature`, the litter's temperature, and its [`LitterVapor`](@ref) network,
as `(; vapor_flux, temperature, vapor)`.
"""
litter_water_exchange(ground::Ground, conditions::GroundConditions, soil_surface_temperature) =
    litter_water_exchange(exposed_litter(ground), ground, conditions, soil_surface_temperature)
litter_water_exchange(::NoLitter, ground::Ground, conditions::GroundConditions, soil_surface_temperature) =
    (; vapor_flux=0.0u"kg/m^2/s", temperature=u"K"(soil_surface_temperature), vapor=no_litter_vapor())
function litter_water_exchange(litter::LitterLayer, ground::Ground, conditions::GroundConditions, soil_surface_temperature)
    (; atmospheric_pressure, vapour_pressure_equation, air) = conditions
    (; temperature, vapor) = litter_state(litter, ground, conditions, soil_surface_temperature)
    exchange_temperature = convection_surface_temperature(air, temperature)
    (; mass_transfer_coefficient) = ground_air_transfer(air, exchange_temperature, atmospheric_pressure, vapour_pressure_equation)
    densities = soil_and_air_vapor_density(conditions, soil_surface_temperature)
    vapor_flux =litter_element_evaporation(vapor, mass_transfer_coefficient, densities.soil, densities.air)
    return (; vapor_flux, temperature, vapor)
end

"""
    litter_water_step(litter_model::AbstractLitterModel, water_content; vapor_flux, temperature, rain, duration, buried)

The litter's water content per dry mass after `duration` of evaporating `vapor_flux` at
`temperature` (from [`litter_water_exchange`](@ref)) and catching `rain`, and the water
passing it to the soil, as `(; water_content, throughfall)`. Water condensing beyond what the
litter holds drips to the soil. Litter buried under snow neither catches rain nor exchanges
water.
"""
litter_water_step(::NoLitter, water_content; rain, kw...) = (; water_content, throughfall=rain)
litter_water_step(litter::LitterLayer, water_content; kw...) = litter_water_step(litter.moisture, litter, water_content; kw...)
litter_water_step(::DryLitter, litter, water_content; rain, kw...) = (; water_content, throughfall=rain)
function litter_water_step(moisture::WetLitter, litter, water_content; vapor_flux, temperature, rain, duration, buried)
    buried && return (; water_content, throughfall=rain)
    (; dry_mass) = litter
    capacity = litter_maximum_water_content(moisture, temperature) * dry_mass
    held = max(water_content * dry_mass - uconvert(unit(rain), vapor_flux * duration), zero(rain))
    drip = max(held - capacity, zero(rain))
    caught = min((1.0 - litter_gap_fraction(litter)) * rain, capacity - (held - drip))
    return (; water_content=uconvert(NoUnits, (held - drip + caught) / dry_mass), throughfall=rain - caught + drip)
end
