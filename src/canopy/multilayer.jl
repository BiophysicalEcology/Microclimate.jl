"""
    MultilayerCanopy(; canopy_height, plant_area_index, woody_area_fraction=0.0,
                       shortwave_model=TwoStreamRadiation(),
                       longwave_model=LayeredLongwaveExchange(),
                       wind_model=CanopyWindAttenuation(),
                       air_profile_model=KTheoryAirProfile(),
                       interception_model=NoInterception(),
                       leaf_parameters=LeafParameters(),
                       stomatal_model=PrescribedStomatalConductance(),
                       leaf_temperature_solver=LinearizedLeafTemperature(),
                       convergence=FixedSoilTemperatureIterations(3))

Multi-layer canopy structure resolved by height. Splits the canopy into as
many layers as there are `MicroModel.heights` entries at or below
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
  exchange (default [`LayeredLongwaveExchange`](@ref))
- `wind_model::AbstractCanopyWindModel` — in-canopy wind attenuation (default
  [`CanopyWindAttenuation`](@ref))
- `air_profile_model::AbstractCanopyAirProfileModel` — in-canopy air-
  temperature profile (default [`KTheoryAirProfile`](@ref))
- `interception_model::AbstractCanopyInterceptionModel` — rain interception
  (default [`NoInterception`](@ref): off)
- `leaf_parameters::LeafParameters` — per-leaf structural/physiological data
- `stomatal_model::AbstractStomatalConductanceModel` — stomatal response
  (default [`PrescribedStomatalConductance`](@ref): day/night gating only).
  Recommended pairing with `CampbellSoilHydraulics`: keep the default here —
  Campbell's own demand/supply solve is already the water-stress mechanism
  once its `leaf_water_potential` output feeds this canopy each hour, so
  [`MoistureResponsiveStomatalConductance`](@ref) would apply the same
  closure curve twice.
- `leaf_temperature_solver::AbstractLeafTemperatureSolver` — how leaf
  temperature is solved each call (default [`LinearizedLeafTemperature`](@ref))
- `convergence::AbstractSoilTemperatureConvergence` — hourly Picard-loop
  convergence criterion (default `FixedSoilTemperatureIterations(3)`); reuses
  the soil spin-up loop's dispatch (generic over any temperature array).

Ground reflectance is not stored here — supplied to
[`canopy_shortwave!`](@ref) by the caller (e.g. `Site.albedo`), matching how
other radiation models read albedo from their caller.
"""
@kwdef struct MultilayerCanopy{H,PAI,WAF,RM,LWM,WM,APM,IM,LP,SM,LTS,CV} <: AbstractCanopyModel
    canopy_height::H
    plant_area_index::PAI
    woody_area_fraction::WAF = 0.0
    shortwave_model::RM = TwoStreamRadiation()
    longwave_model::LWM = LayeredLongwaveExchange()
    wind_model::WM = CanopyWindAttenuation()
    air_profile_model::APM = KTheoryAirProfile()
    interception_model::IM = NoInterception()
    leaf_parameters::LP = LeafParameters()
    stomatal_model::SM = PrescribedStomatalConductance()
    leaf_temperature_solver::LTS = LinearizedLeafTemperature()
    convergence::CV = FixedSoilTemperatureIterations(3)
end

function example_multilayer_canopy(;
    canopy_height = 1.0u"m",
    plant_area_index = 3.0,
    woody_area_fraction = 0.0,
    shortwave_model = TwoStreamRadiation(; leaf_reflectance = 0.1),
    longwave_model = LayeredLongwaveExchange(),
    wind_model = CanopyWindAttenuation(),
    air_profile_model = KTheoryAirProfile(),
    interception_model = NoInterception(),
    leaf_parameters = LeafParameters(),
    stomatal_model = PrescribedStomatalConductance(),
    leaf_temperature_solver = LinearizedLeafTemperature(),
    convergence = FixedSoilTemperatureIterations(3),
)
    MultilayerCanopy(;
        canopy_height, plant_area_index, woody_area_fraction, shortwave_model, longwave_model, wind_model,
        air_profile_model, interception_model, leaf_parameters, stomatal_model, leaf_temperature_solver, convergence,
    )
end

# Number of canopy layers: one per `heights` entry at or below canopy_height.
n_canopy_layers(model::MultilayerCanopy, heights) = count(h -> h <= model.canopy_height, heights)

"""
    allocate_canopy(model::MultilayerCanopy, heights, boundary_layer_model)

Buffers nest one sub-struct per sub-model — `buffers.shortwave`,
`buffers.longwave`, `buffers.wind`, `buffers.air_profile`,
`buffers.interception`, `buffers.leaf` — mirroring `MicroBuffers`'s own
composition, rather than flattening every sub-model's fields into one
`NamedTuple` (which would let same-named fields silently collide).
`buffers.leaf` bundles orchestrator-owned per-layer state: leaf body,
temperature, and the source terms other sub-models read/write.
"""
function allocate_canopy(model::MultilayerCanopy, heights, boundary_layer_model)
    n_layers = n_canopy_layers(model, heights)
    leaf_angle_distribution_parameter = model.leaf_parameters.leaf_angle_distribution_parameter

    return (;
        shortwave = allocate_shortwave(model.shortwave_model, model.canopy_height, model.plant_area_index, n_layers, leaf_angle_distribution_parameter),
        longwave = allocate_longwave(model.longwave_model, model.plant_area_index, n_layers),
        wind = allocate_wind(model.wind_model, model.canopy_height, model.plant_area_index, boundary_layer_model, heights, n_layers),
        air_profile = allocate_air_profile(model.air_profile_model, model.canopy_height, model.plant_area_index, heights, n_layers, boundary_layer_model),
        interception = allocate_interception(model.interception_model, model.canopy_height, model.plant_area_index, heights, n_layers, boundary_layer_model),
        leaf = (;
            leaf_body = leaf_body(model.leaf_parameters.leaf_length, model.leaf_parameters.leaf_width),
            leaf_temperature = zeros(typeof(0.0u"K"), n_layers),
            leaf_temperature_prev = zeros(typeof(0.0u"K"), n_layers),
            sensible_heat_source = zeros(typeof(0.0u"W/m^2"), n_layers),
            evaporation_mass_flow = zeros(typeof(0.0u"g/s"), n_layers),
            potential_evaporation_mass_flow = zeros(typeof(0.0u"g/s"), n_layers),
        ),
    )
end

canopy_shortwave!(buffers, model::MultilayerCanopy; kw...) =
    canopy_shortwave!(buffers.shortwave, model.shortwave_model, model.plant_area_index; kw...)

canopy_longwave!(buffers, model::MultilayerCanopy; kw...) =
    canopy_longwave!(buffers.longwave, model.longwave_model, model.leaf_parameters.leaf_emissivity; kw...)

canopy_wind_profile!(buffers, model::MultilayerCanopy, boundary_layer_model; kw...) =
    canopy_wind_profile!(buffers.wind, model.wind_model, boundary_layer_model; kw...)

canopy_air_profile!(buffers, model::MultilayerCanopy, boundary_layer_model; kw...) =
    canopy_air_profile!(buffers.air_profile, model.air_profile_model, boundary_layer_model; kw...)

canopy_interception!(buffers, model::MultilayerCanopy; kw...) =
    canopy_interception!(buffers.interception, model.interception_model, model.leaf_parameters.leaf_angle_distribution_parameter; kw...)

wet_canopy_fraction(model::MultilayerCanopy, buffers, layer) =
    wet_canopy_fraction(model.interception_model, buffers.interception, layer)

deplete_canopy_water!(model::MultilayerCanopy, buffers, layer, evaporated_mass) =
    deplete_canopy_water!(model.interception_model, buffers.interception, layer, evaporated_mass)

reset_canopy_scratch!(model::MultilayerCanopy, buffers) =
    reset_interception!(model.interception_model, buffers.interception)

# Bootstrap values only — refreshed before first real use by apply_canopy_overrides.
function allocate_canopy_inputs(model::MultilayerCanopy; site, environment_instant, boundary_layer_model)
    CanopyEnergyBalanceInputs(model;
        site, environment_instant, zenith_angle = environment_instant.zenith_angle,
        direct_horizontal_irradiance = 0.0u"W/m^2", diffuse_horizontal_irradiance = 0.0u"W/m^2",
        ground_reflectance = _baseline_albedo(site, environment_instant),
        ground_temperature = environment_instant.reference_temperature,
        ground_emissivity = environment_instant.surface_emissivity,
        ground_relative_humidity = environment_instant.reference_humidity,
        canopy_source_temperature = environment_instant.reference_temperature,
    )
end

initial_ground_overrides(::MultilayerCanopy) = (; ground_shortwave_transmission = 1.0, ground_incoming_longwave = 0.0u"W/m^2")

canopy_leaf_area_index(model::MultilayerCanopy) = sum(model.plant_area_index) * (1 - model.woody_area_fraction)
