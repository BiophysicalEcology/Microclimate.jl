"""
    AbstractGroundLongwaveAbsorption

How much of the incoming longwave the ground absorbs without a canopy:
[`EmissivityAbsorption`](@ref) or [`FullAbsorption`](@ref).
"""
abstract type AbstractGroundLongwaveAbsorption end

"""
    EmissivityAbsorption()

The ground absorbs its emissivity of the incoming longwave and reflects the rest.
"""
struct EmissivityAbsorption <: AbstractGroundLongwaveAbsorption end

"""
    FullAbsorption()

The ground absorbs all the incoming longwave, as in NicheMapR.
"""
struct FullAbsorption <: AbstractGroundLongwaveAbsorption end

"""
    absorbed_ground_longwave(absorption::AbstractGroundLongwaveAbsorption, emissivity, incoming_longwave)

The longwave the ground absorbs from `incoming_longwave`, at `emissivity`.
"""
absorbed_ground_longwave(::EmissivityAbsorption, emissivity, incoming_longwave) = emissivity * incoming_longwave
absorbed_ground_longwave(::FullAbsorption, emissivity, incoming_longwave) = incoming_longwave

"""
    ViewFactorLongwave(; atmospheric_radiation_model=CampbellNormanAtmosphericRadiation(),
        ground_absorption=EmissivityAbsorption())

View factor-based surface-longwave budget (as in NicheMapR). Combines the clear-sky
downwelling longwave from `atmospheric_radiation_model` with cloud
emissivity, sky view fraction, vegetation shade, and hillshade contributions
to give the net longwave at the surface. `ground_absorption`, an
[`AbstractGroundLongwaveAbsorption`](@ref), is how much of it the ground absorbs.

The `atmospheric_radiation_model` lives on the longwave model rather than on
`MicroConfig` because no other longwave model consumes it.
"""
@kwdef struct ViewFactorLongwave{ARM<:AbstractAtmosphericRadiationModel,GA<:AbstractGroundLongwaveAbsorption} <: AbstractLongwaveModel
    atmospheric_radiation_model::ARM = CampbellNormanAtmosphericRadiation()
    ground_absorption::GA = EmissivityAbsorption()
end

"""
    downward_longwave(atmospheric_radiation_model; air_temperature, relative_humidity, atmospheric_pressure,
        cloud_cover, cloud_emissivity, sky_view_fraction, supplied_longwave=nothing, vapour_pressure_equation=GoffGratch())

The longwave reaching open ground, as `(; sky, hillside, open_ground)`. `sky` is `supplied_longwave` when
given, else the clear sky's and the clouds' longwave. `hillside` is a black body at
`air_temperature`. `open_ground` is `sky` over the `sky_view_fraction` and `hillside` over the rest.
"""
function downward_longwave(atmospheric_radiation_model; air_temperature, relative_humidity, atmospheric_pressure,
    cloud_cover, cloud_emissivity, sky_view_fraction, supplied_longwave=nothing, vapour_pressure_equation=GoffGratch(),
)
    temperature = u"K"(air_temperature)
    sky = if isnothing(supplied_longwave)
        modelled_sky_longwave(atmospheric_radiation_model, temperature, relative_humidity, atmospheric_pressure,
            cloud_cover, cloud_emissivity, vapour_pressure_equation)
    else
        supplied_longwave
    end
    hillside = σ * temperature^4 # TODO this is assuming unit emissivity, but should have a vegetation emissivity term
    return (; sky, hillside, open_ground = sky * sky_view_fraction + hillside * (1.0 - sky_view_fraction))
end

function modelled_sky_longwave(atmospheric_radiation_model, temperature, relative_humidity, atmospheric_pressure,
    cloud_cover, cloud_emissivity, vapour_pressure_equation,
)
    vapour_pressure = wet_air_properties(temperature, relative_humidity, atmospheric_pressure; vapour_pressure_equation).vapour_pressure
    clear_sky = atmospheric_radiation(atmospheric_radiation_model, vapour_pressure, temperature)
    cloud = σ * cloud_emissivity * (temperature - 2.0u"K")^4
    return clear_sky * (1.0 - cloud_cover) + cloud * cloud_cover
end

precompute_longwave_sky(model::ViewFactorLongwave; kw...) =
    precompute_longwave_sky(model.atmospheric_radiation_model; kw...)

"""
    precompute_longwave_sky(atmospheric_radiation_model; site, environment_instant,
                             vapour_pressure_equation, surface_emissivity, shade)

Clear-sky + cloud + hillshade downwelling longwave, keyed directly on the
`AbstractAtmosphericRadiationModel` rather than a full [`ViewFactorLongwave`](@ref)
— lets other longwave schemes (e.g. the canopy's per-layer exchange) reuse
this exact sky calculation without needing a `ViewFactorLongwave` wrapper.

If `environment_instant.longwave_radiation` is set, uses it directly as the
sky term instead of deriving one from `cloud_cover`.
"""
function precompute_longwave_sky(atmospheric_radiation_model::AbstractAtmosphericRadiationModel;
    site,
    environment_instant,
    vapour_pressure_equation=GoffGratch(),
    # Optional explicit overrides so callers (e.g. snow surface) can pass
    # specific surface_emissivity / shade without rebuilding environment_instant.
    surface_emissivity=environment_instant.surface_emissivity,
    shade=environment_instant.shade,
)
    (; atmospheric_pressure, reference_humidity, reference_temperature, cloud_emissivity, cloud_cover) = environment_instant
    (; sky_view_fraction) = site
    downward = downward_longwave(atmospheric_radiation_model; air_temperature=reference_temperature,
        relative_humidity=reference_humidity, atmospheric_pressure, cloud_cover, cloud_emissivity, sky_view_fraction,
        supplied_longwave=get(environment_instant, :longwave_radiation, nothing), vapour_pressure_equation)
    hillshade_radiation = downward.hillside
    longwave_radiation_sky = downward.sky * (1.0 - shade)
    longwave_radiation_vegetation = shade * hillshade_radiation
    longwave_radiation_hillshade = hillshade_radiation

    incoming_longwave = (longwave_radiation_sky + longwave_radiation_vegetation) * sky_view_fraction +
                        longwave_radiation_hillshade * (1.0 - sky_view_fraction)
    outgoing_coeff = (1.0 - shade) * σ * surface_emissivity
    ground_shade_term = shade * hillshade_radiation
    sky_temperature = (incoming_longwave / σ)^(1//4)

    return (;
        incoming_longwave,
        open_downward_longwave = downward.open_ground,
        outgoing_coeff,
        ground_shade_term,
        sky_temperature=u"K"(sky_temperature),
        longwave_radiation_sky,
        longwave_radiation_vegetation,
        longwave_radiation_hillshade,
    )
end

function longwave_radiation(model::ViewFactorLongwave;
    site,
    environment_instant,
    surface_temperature,
    vapour_pressure_equation=GoffGratch(),
)
    sky = precompute_longwave_sky(model; site, environment_instant, vapour_pressure_equation)
    (; incoming_longwave, ground_shade_term) = sky

    surface_radiation = σ * environment_instant.surface_emissivity * (u"K"(surface_temperature))^4
    shade = environment_instant.shade
    longwave_radiation_ground = (1 - shade) * surface_radiation + ground_shade_term
    net_longwave_radiation = incoming_longwave - (1 - shade) * surface_radiation - ground_shade_term

    return (;
        sky.sky_temperature,
        sky.longwave_radiation_sky,
        sky.longwave_radiation_vegetation,
        sky.longwave_radiation_hillshade,
        longwave_radiation_ground,
        net_longwave_radiation,
    )
end
