"""
    MultilayerCanopy(; canopy_height, plant_area_index, woody_area_fraction=0.0,
                       shortwave_model=TwoStreamRadiation(),
                       longwave_model=LayeredRadiosityExchange(),
                       wind_model=ExponentialCanopyWindAttenuation(),
                       air_profile_model=RaupachLTheoryAirProfile(),
                       interception_model=NoInterception(),
                       condensation_model=MonteithLeafCondensation(),
                       leaf_parameters=LeafParameters(),
                       stomatal_model=PrescribedStomatalConductance(),
                       leaf_convection_model=ElaborateLeafConvection(),
                       convergence_model=NonlinearCanopyConvergence())

Multi-layer canopy structure resolved by height. Splits the canopy into as
many layers as there are `MicroModel.atmosphere.heights` entries at or below
`canopy_height`.

Composes swappable sub-models the same way `MicroModel` composes
`RadiationModel`/`boundary_layer_model`/`soil_energy_model`:

- `canopy_height`, `plant_area_index` — solver geometry shared by every
  sub-model, not owned by any one of them. `plant_area_index` (PAI: leaves
  plus woody material) is a scalar total (split evenly across layers) or an
  `AbstractVector` of per-layer values, top-to-bottom (see
  [`plant_area_index_from_density`](@ref) for a density-profile builder).
- `woody_area_fraction` — fraction of `plant_area_index` that is woody
  rather than leaf. Every sub-model here uses the full PAI (bark intercepts
  light/rain/momentum too); only [`canopy_leaf_area_index`](@ref)'s Campbell
  hookup narrows PAI down to LAI with it. Default `0.0` (PAI ≈ LAI).
- `shortwave_model::AbstractCanopyShortwaveModel` — layer-resolved shortwave
  radiative transfer (default [`TwoStreamRadiation`](@ref), Dickinson/Sellers)
- `longwave_model::AbstractCanopyLongwaveModel` — layer-resolved longwave
  exchange (default [`LayeredRadiosityExchange`](@ref); see also
  [`LayeredLongwaveExchange`](@ref))
- `wind_model::AbstractCanopyWindModel` — in-canopy wind attenuation (default
  [`ExponentialCanopyWindAttenuation`](@ref); see also [`MixingLengthCanopyWindAttenuation`](@ref))
- `air_profile_model::AbstractCanopyAirProfileModel` — in-canopy air-
  temperature profile (default [`RaupachLTheoryAirProfile`](@ref); see also
  [`KTheoryAirProfile`](@ref))
- `ground_coupling::AbstractGroundCoupling` — ground resistance (default
  [`FloorTurbulence`](@ref); see also [`SakaguchiZeng2009`](@ref), [`ProfileGroundResistance`](@ref))
- `roughness_sublayer::AbstractRoughnessSublayer` — mixing just above canopy top
  (default [`NoRoughnessSublayer`](@ref))
- `interception_model::AbstractCanopyInterceptionModel` — rain interception
  (default [`NoInterception`](@ref): off)
- `condensation_model::AbstractCondensationModel` — leaf dew/frost formation
  (default [`MonteithLeafCondensation`](@ref); [`NoCondensation`](@ref) to
  switch it off). Independent of `MicroModel.surface.condensation_model` (ground).
- `leaf_parameters::LeafParameters` — per-leaf structural/physiological data
- `stomatal_model::AbstractStomatalConductanceModel` — stomatal response
  (default [`PrescribedStomatalConductance`](@ref): day/night gating only,
  ignores `leaf_water_potential`). Pair with `CampbellSoilHydraulics` via
  [`MoistureResponsiveStomatalConductance`](@ref) instead if soil-moisture
  stress should actually reach the leaf's transpiration/temperature —
  `PrescribedStomatalConductance` leaves Campbell's `leaf_water_potential`
  output unused.
- `leaf_convection_model::AbstractLeafConvectionModel` — leaf boundary-layer
  Nusselt correlation (default [`ElaborateLeafConvection`](@ref); see also
  [`SimpleLeafConvection`](@ref))
- `convergence_model::AbstractCanopyConvergenceModel` — how the hourly
  leaf/air-temperature state is solved (default
  [`NonlinearCanopyConvergence`](@ref))

Ground reflectance is not stored here — supplied to
[`canopy_shortwave!`](@ref) by the caller (e.g. `Site.albedo`), matching how
other radiation models read albedo from their caller.
"""
@kwdef struct MultilayerCanopy{H,PAI,WAF,RM,LWM,WM,APM,GC,RS,IM,CDM,LP,SM,LCM,CVM} <: AbstractCanopyModel
    canopy_height::H
    plant_area_index::PAI
    woody_area_fraction::WAF = 0.0
    shortwave_model::RM = TwoStreamRadiation()
    longwave_model::LWM = LayeredRadiosityExchange()
    wind_model::WM = ExponentialCanopyWindAttenuation()
    air_profile_model::APM = RaupachLTheoryAirProfile()
    ground_coupling::GC = FloorTurbulence()
    roughness_sublayer::RS = NoRoughnessSublayer()
    interception_model::IM = NoInterception()
    condensation_model::CDM = MonteithLeafCondensation()
    leaf_parameters::LP = LeafParameters()
    stomatal_model::SM = PrescribedStomatalConductance()
    leaf_convection_model::LCM = ElaborateLeafConvection()
    convergence_model::CVM = NonlinearCanopyConvergence()
end

"""
    example_multilayer_canopy(; kwargs...)

Example [`MultilayerCanopy`](@ref), 1 m tall with a plant area index of 3. `kwargs`
override any field.
"""
function example_multilayer_canopy(;
    canopy_height = 1.0u"m",
    plant_area_index = 3.0,
    woody_area_fraction = 0.0,
    shortwave_model = TwoStreamRadiation(; leaf_reflectance = 0.1),
    longwave_model = LayeredRadiosityExchange(),
    wind_model = ExponentialCanopyWindAttenuation(),
    air_profile_model = RaupachLTheoryAirProfile(),
    ground_coupling = FloorTurbulence(),
    roughness_sublayer = NoRoughnessSublayer(),
    interception_model = NoInterception(),
    condensation_model = MonteithLeafCondensation(),
    leaf_parameters = LeafParameters(),
    stomatal_model = PrescribedStomatalConductance(),
    leaf_convection_model = ElaborateLeafConvection(),
    convergence_model = NonlinearCanopyConvergence(),
)
    MultilayerCanopy(;
        canopy_height, plant_area_index, woody_area_fraction, shortwave_model, longwave_model, wind_model,
        air_profile_model, ground_coupling, roughness_sublayer, interception_model, condensation_model, leaf_parameters, stomatal_model,
        leaf_convection_model, convergence_model,
    )
end

# Number of canopy layers: one per `heights` entry at or below canopy_height.
n_canopy_layers(model::MultilayerCanopy, heights) = count(h -> h <= model.canopy_height, heights)

"""
    allocate_canopy(model::MultilayerCanopy, heights, boundary_layer_model; site, vapour_pressure_equation, soil_hydraulic_model)

Buffers nest one sub-struct per sub-model — `buffers.shortwave`,
`buffers.longwave`, `buffers.wind`, `buffers.air_profile`,
`buffers.interception`, `buffers.leaf` — mirroring `MicroBuffers`'s own
composition, rather than flattening every sub-model's fields into one
`NamedTuple` (which would let same-named fields silently collide).
`buffers.leaf` bundles orchestrator-owned per-layer state: leaf body,
temperature, and the source terms other sub-models read/write.
"""
function allocate_canopy(model::MultilayerCanopy, heights, boundary_layer_model; site, vapour_pressure_equation, soil_hydraulic_model)
    n_layers = n_canopy_layers(model, heights)
    lowest_layer_height = minimum(h for h in heights if h <= model.canopy_height)
    lowest_layer_height > site.roughness_height || throw(ArgumentError(
        "lowest canopy layer height ($lowest_layer_height) must be above the ground roughness height ($(site.roughness_height))"))
    canopy_projection_ratio = model.leaf_parameters.canopy_projection_ratio
    wind = allocate_wind(model.wind_model, model.canopy_height, model.plant_area_index, boundary_layer_model, heights, n_layers)

    buffers = (;
        shortwave = allocate_shortwave(model.shortwave_model, model.canopy_height, model.plant_area_index, n_layers, canopy_projection_ratio),
        longwave = allocate_longwave(model.longwave_model, model.plant_area_index, n_layers, canopy_projection_ratio),
        wind,
        roughness_sublayer = (; enhancement=canopy_top_enhancement(model.roughness_sublayer, model.canopy_height,
            sum(model.plant_area_index), site.roughness_height, boundary_layer_model.karman_constant)),
        air_profile = allocate_air_profile(model.air_profile_model, model.canopy_height, wind.displacement_height,
            model.plant_area_index, heights, n_layers, boundary_layer_model),
        interception = allocate_interception(model.interception_model, model.canopy_height, model.plant_area_index, heights, n_layers, boundary_layer_model),
        leaf = (;
            leaf_body = leaf_body(model.leaf_parameters.leaf_length, model.leaf_parameters.leaf_width, canopy_projection_ratio),
            leaf_temperature = zeros(typeof(0.0u"K"), n_layers),
            sensible_heat_source = zeros(typeof(0.0u"W/m^2"), n_layers),
            evaporation_mass_flow = zeros(typeof(0.0u"g/s"), n_layers),
            potential_evaporation_mass_flow = zeros(typeof(0.0u"g/s"), n_layers),
            leaf_dew_formed = zeros(typeof(0.0u"kg/m^2"), n_layers),
            leaf_frost_formed = zeros(typeof(0.0u"kg/m^2"), n_layers),
            absorbed_radiation = zeros(typeof(0.0u"W/m^2"), n_layers),
            net_balance = zeros(typeof(0.0u"W"), n_layers),
        ),
    )
    convergence = allocate_canopy_convergence(model.convergence_model, buffers, model, boundary_layer_model;
        site, vapour_pressure_equation, soil_hydraulic_model)
    return (; buffers..., convergence)
end

canopy_shortwave!(buffers, model::MultilayerCanopy; kw...) =
    canopy_shortwave!(buffers.shortwave, model.shortwave_model, model.plant_area_index; kw...)

canopy_longwave!(buffers, model::MultilayerCanopy; kw...) =
    canopy_longwave!(buffers.longwave, model.longwave_model, model.leaf_parameters.leaf_emissivity; kw...)

canopy_incoming_longwave(model::MultilayerCanopy; kw...) = canopy_incoming_longwave(model.longwave_model; kw...)

canopy_wind_profile!(buffers, model::MultilayerCanopy, boundary_layer_model; kw...) =
    canopy_wind_profile!(buffers.wind, model.wind_model, boundary_layer_model; kw...)

canopy_interception!(buffers, model::MultilayerCanopy; kw...) =
    canopy_interception!(buffers.interception, model.interception_model, model.leaf_parameters.canopy_projection_ratio; kw...)

wet_canopy_fraction(model::MultilayerCanopy, buffers, layer) =
    wet_canopy_fraction(model.interception_model, buffers.interception, layer)

available_canopy_water(model::MultilayerCanopy, buffers, layer) =
    available_canopy_water(model.interception_model, buffers.interception, layer)

canopy_water_capacity(model::MultilayerCanopy, buffers, layer) =
    canopy_water_capacity(model.interception_model, buffers.interception, layer)

deplete_canopy_water!(model::MultilayerCanopy, buffers, layer, evaporated_mass) =
    deplete_canopy_water!(model.interception_model, buffers.interception, layer, evaporated_mass)

leaf_standing_dew(model::MultilayerCanopy, buffers, layer) =
    leaf_standing_dew(model.interception_model, buffers.interception, layer)

leaf_standing_frost(model::MultilayerCanopy, buffers, layer) =
    leaf_standing_frost(model.interception_model, buffers.interception, layer)

add_leaf_standing_dew!(model::MultilayerCanopy, buffers, layer, Δ) =
    add_leaf_standing_dew!(model.interception_model, buffers.interception, layer, Δ)

add_leaf_standing_frost!(model::MultilayerCanopy, buffers, layer, Δ) =
    add_leaf_standing_frost!(model.interception_model, buffers.interception, layer, Δ)

clamp_leaf_standing_dew!(model::MultilayerCanopy, buffers, layer, ceiling) =
    clamp_leaf_standing_dew!(model.interception_model, buffers.interception, layer, ceiling)

snapshot_canopy_water!(model::MultilayerCanopy, buffers) =
    snapshot_interception!(model.interception_model, buffers.interception)

restore_canopy_water!(model::MultilayerCanopy, buffers) =
    restore_interception!(model.interception_model, buffers.interception)

function reset_canopy_scratch!(model::MultilayerCanopy, buffers)
    reset_interception!(model.interception_model, buffers.interception)
    return nothing
end

initial_ground_overrides(::MultilayerCanopy) = (;
    ground_shortwave_transmission = 1.0, ground_incoming_longwave = 0.0u"W/m^2",
    ground_wind_speed = 0.0u"m/s", ground_air_temperature = 288.0u"K",
    ground_air_relative_humidity = 0.5, ground_reference_height = 1.0u"m",
    ground_heat_conductance = 0.0u"W/m^2/K", ground_vapor_conductance = 0.0u"m/s",
)

canopy_leaf_area_index(model::MultilayerCanopy) = sum(model.plant_area_index) * (1 - model.woody_area_fraction)

# Canopy-displaced log-profile origin for solve_air!'s output.profile.
canopy_profile_origin(::MultilayerCanopy, buffers, site, boundary_layer_model) =
    (; displacement_height = buffers.wind.displacement_height, roughness_length = buffers.wind.roughness_length,
        thermal_roughness_model = buffers.wind.thermal_roughness_model,
        temperature_anchor_height = buffers.wind.canopy_height)

# Top-of-canopy air temperature for solve_air!'s output.profile above canopy.
profile_surface_temperature(::MultilayerCanopy, output, i, ground_surface_temperature) =
    output.canopy.air_temperature[i, 1]

# Top-of-canopy relative humidity, RHabove's counterpart to profile_surface_temperature.
profile_surface_humidity(::MultilayerCanopy, output, i) = output.canopy.canopy_top_relative_humidity[i]
