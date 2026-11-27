"""
    AbstractLitterModel

Plant litter lying on the soil, for `SurfaceModel.litter_model`: [`NoLitter`](@ref) or
[`LitterLayer`](@ref).
"""
abstract type AbstractLitterModel end

"""
    NoLitter()

Bare soil: the soil surface exchanges radiation, heat and water vapor with the air directly.
"""
struct NoLitter <: AbstractLitterModel end

# The litter the air sees: none while snow buries it
exposed_litter(ground::Ground) = ground.snow_covered ? NoLitter() : ground.litter

# The water per dry mass the litter holds when the simulation starts
initial_litter_water_content(::NoLitter) = 0.0

"""
    soil_surface_heat_flux(litter::AbstractLitterModel, ground::Ground, conditions::GroundConditions,
        evaporation_model, soil_surface_temperature)

Heat flux into the soil surface from above (W/m²): net radiation, sensible heat and the
latent heat of soil evaporation, through `litter`.
"""
function soil_surface_heat_flux(::NoLitter, ground::Ground, conditions::GroundConditions, evaporation_model,
    soil_surface_temperature,
)
    (; direct_shortwave, diffuse_shortwave, downward_longwave, exposed_fraction, atmospheric_pressure,
       vapour_pressure_equation, air) = conditions
    absorbed_shortwave = (1.0 - ground.soil_albedo) * (direct_shortwave + diffuse_shortwave)
    (; soil_emissivity, soil_longwave_absorption) = ground
    net_longwave = exposed_fraction * (absorbed_ground_longwave(soil_longwave_absorption, soil_emissivity, downward_longwave) -
        σ * soil_emissivity * u"K"(soil_surface_temperature)^4)
    exchange_temperature = convection_surface_temperature(air, soil_surface_temperature)
    (; sensible_heat_flux, mass_transfer_coefficient) = ground_air_transfer(air, exchange_temperature,
        atmospheric_pressure, vapour_pressure_equation)
    latent_heat_flux = ground_latent_heat_flux(evaporation_model, ground, conditions, exchange_temperature, mass_transfer_coefficient)
    return absorbed_shortwave + net_longwave + sensible_heat_flux - latent_heat_flux
end

"""
    ground_reflectance(ground::Ground, sun_elevation_sine)

Shortwave reflectance of the ground: one number, or a [`GroundShortwaveReflectance`](@ref) for
the unscattered beam and diffuse light.
"""
ground_reflectance(ground::Ground, sun_elevation_sine) = ground_reflectance(exposed_litter(ground), ground, sun_elevation_sine)
ground_reflectance(::NoLitter, ground::Ground, sun_elevation_sine) = ground.soil_albedo

"""
    ground_surface(ground::Ground, conditions::GroundConditions, soil_surface_temperature)

The [`GroundSurface`](@ref) over a soil surface at `soil_surface_temperature`.
"""
ground_surface(ground::Ground, conditions::GroundConditions, soil_surface_temperature) =
    ground_surface(exposed_litter(ground), ground, conditions, soil_surface_temperature)
function ground_surface(::NoLitter, ground::Ground, conditions::GroundConditions, soil_surface_temperature)
    (; soil_emissivity, soil_wetness) = ground
    (; atmospheric_pressure, vapour_pressure_equation, sun_elevation_sine, air) = conditions
    saturated = wet_air_properties(u"K"(soil_surface_temperature), 1.0, atmospheric_pressure; vapour_pressure_equation).vapour_density
    air_vapor_density = wet_air_properties(u"K"(air.air_temperature), air.air_relative_humidity, atmospheric_pressure;
        vapour_pressure_equation).vapour_density
    # The wet fraction of the soil evaporates as if saturated, as ground_latent_heat_flux has it
    vapor_density = air_vapor_density + soil_wetness * (saturated - air_vapor_density)
    return GroundSurface(; temperature=soil_surface_temperature, longwave_reflectance=soil_longwave_reflectance(ground),
        emitted_longwave=soil_emissivity * σ * u"K"(soil_surface_temperature)^4,
        reflectance=ground_reflectance(NoLitter(), ground, sun_elevation_sine), vapor_density)
end

# The fraction of the longwave reaching the soil that it reflects
soil_longwave_reflectance(ground::Ground) =
    1.0 - absorbed_ground_longwave(ground.soil_longwave_absorption, ground.soil_emissivity, 1.0)
