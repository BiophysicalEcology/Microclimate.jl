"""
    KTheoryAirProfile()

Air temperature and vapor density through the canopy, from steady diffusion between
neighbouring layers: `d/dz(K dT/dz) = -S`, with each layer's leaves as the source `S`.

Layer 1 is the canopy-top value; each layer exchanges only with its neighbours. The
diffusivity `K`, shared by heat and vapor, follows the wind profile, scaled at the top to
von Kármán constant × friction velocity × (`canopy_height` − `displacement_height`).

# References
- Raupach, M. R. (1989). A practical Lagrangian method for relating scalar
  concentrations to source distributions in vegetation canopies. *Quarterly
  Journal of the Royal Meteorological Society*, 115(487), 609-632.
"""
struct KTheoryAirProfile <: AbstractCanopyAirProfileModel end

function allocate_air_profile(::KTheoryAirProfile, canopy_height, displacement_height, plant_area_index, heights, n_layers, boundary_layer_model)
    (; layer_thickness) = canopy_layer_heights(heights, canopy_height, n_layers)
    (; layer_plant_area_index) = canopy_layer_geometry(plant_area_index, n_layers)
    relative_eddy_diffusivity = wind_attenuation_profile(
        layer_plant_area_index, canopy_height, sum(plant_area_index), boundary_layer_model.karman_constant,
    )

    # Layer 1 is fixed, so only the layers below it are solved
    free_layer_count = n_layers - 1
    heat_conductance = typeof(0.0u"W/m^2/K")
    heat_flux = typeof(0.0u"W/m^2")
    vapor_conductance = typeof(0.0u"m/s")
    vapor_flux = typeof(0.0u"kg/m^2/s")

    return (;
        canopy_height, displacement_height,
        layer_spacing = layer_thickness, relative_eddy_diffusivity,
        heat_layer_total_conductance = zeros(heat_conductance, free_layer_count),
        heat_conductance_between_layers = zeros(heat_conductance, free_layer_count - 1),
        heat_known_flux = zeros(heat_flux, free_layer_count),
        vapor_layer_total_conductance = zeros(vapor_conductance, free_layer_count),
        vapor_conductance_between_layers = zeros(vapor_conductance, free_layer_count - 1),
        vapor_known_flux = zeros(vapor_flux, free_layer_count),
        air_temperature = zeros(typeof(0.0u"K"), n_layers),
        current_air_temperature = zeros(typeof(0.0u"K"), n_layers),
        vapor_density = zeros(typeof(0.0u"kg/m^3"), n_layers),
        current_vapor_density = zeros(typeof(0.0u"kg/m^3"), n_layers),
        relative_humidity = zeros(n_layers),
    )
end

# The floor stops the diffusivity reaching zero when displacement reaches the canopy top
_ktheory_canopy_top_eddy_diffusivity(buffers, boundary_layer_model, friction_velocity) =
    boundary_layer_model.karman_constant * friction_velocity * max(buffers.canopy_height - buffers.displacement_height, 1.0e-3u"m")

# The lowest layer's spacing over its eddy diffusivity
function ground_turbulent_resistance(buffers, ::KTheoryAirProfile, boundary_layer_model; friction_velocity, inverse_obukhov_length)
    eddy_diffusivity = _ktheory_canopy_top_eddy_diffusivity(buffers, boundary_layer_model, friction_velocity) *
        last(buffers.relative_eddy_diffusivity)
    return last(buffers.layer_spacing) / eddy_diffusivity
end

# Writes the conductance between each pair of free layers into `conductance_between_layers`, and
# returns the conductance from the fixed layer 1 to the first free layer. `scale` turns the
# diffusivity over the spacing into the flux's units.
function _ktheory_conductances!(conductance_between_layers, relative_eddy_diffusivity, canopy_top_eddy_diffusivity, layer_spacing, scale)
    free_layer_count = length(conductance_between_layers) + 1
    top_conductance = convert(eltype(conductance_between_layers), scale * canopy_top_eddy_diffusivity * 0.5 *
        (relative_eddy_diffusivity[1] + relative_eddy_diffusivity[2]) / layer_spacing[1])
    for i in 1:(free_layer_count - 1)
        boundary_eddy_diffusivity = canopy_top_eddy_diffusivity * 0.5 *
            (relative_eddy_diffusivity[i + 1] + relative_eddy_diffusivity[i + 2])
        conductance_between_layers[i] = scale * boundary_eddy_diffusivity / layer_spacing[i + 1]
    end
    return top_conductance
end

# Fills in and solves one flux balance on the free layers, writing the answer into `profile`:
#   g₋(C₋ - Cᵢ) - g₊(Cᵢ - C₊) + Sᵢ = 0
# `source(layer)` is Sᵢ, and layer 1 and the ground give the known end values. `known_flux` is
# what doesn't depend on the unknown values: each layer's source, and the flow from the canopy top
# and the ground into the layers next to them.
function _ktheory_solve!(profile, layer_total_conductance, conductance_between_layers, known_flux, top_conductance, ground_conductance, top_value, ground_value, source)
    free_layer_count = length(layer_total_conductance)
    for i in 1:free_layer_count
        conductance_above = i == 1 ? top_conductance : conductance_between_layers[i - 1]
        conductance_below = i == free_layer_count ? ground_conductance : conductance_between_layers[i]
        layer_total_conductance[i] = -(conductance_above + conductance_below)
        known_flux[i] = -source(i + 1)
    end
    known_flux[1] -= top_conductance * top_value
    known_flux[free_layer_count] -= ground_conductance * ground_value
    # Eliminate from the top layer down, overwriting the totals and known fluxes, then substitute back up
    for i in 1:(free_layer_count - 1)
        ratio = conductance_between_layers[i] / layer_total_conductance[i]
        layer_total_conductance[i + 1] -= ratio^2 * layer_total_conductance[i]
        known_flux[i + 1] -= ratio * known_flux[i]
    end
    profile[free_layer_count + 1] = known_flux[free_layer_count] / layer_total_conductance[free_layer_count]
    for i in (free_layer_count - 1):-1:1
        profile[i + 1] = known_flux[i] / layer_total_conductance[i] - conductance_between_layers[i] / layer_total_conductance[i] * profile[i + 2]
    end
    return profile
end

function canopy_heat_vapor_transport!(buffers, ::KTheoryAirProfile, boundary_layer_model,
    forcing::CanopyAirProfileForcing,
)
    (; friction_velocity, canopy_top_air_temperature, canopy_top_vapor_density, ground_temperature, ground_vapor_density,
       ground_resistance, air_volumetric_heat_capacity, sensible_heat_flux, evaporation_mass_flux) = forcing
    (; layer_spacing, relative_eddy_diffusivity,
       heat_layer_total_conductance, heat_conductance_between_layers, heat_known_flux,
       vapor_layer_total_conductance, vapor_conductance_between_layers, vapor_known_flux,
       air_temperature, vapor_density) = buffers

    layer_count = length(air_temperature)
    length(sensible_heat_flux) == layer_count || throw(ArgumentError("sensible_heat_flux must have one entry per canopy layer"))
    length(evaporation_mass_flux) == layer_count || throw(ArgumentError("evaporation_mass_flux must have one entry per canopy layer"))

    air_temperature[1] = canopy_top_air_temperature
    vapor_density[1] = canopy_top_vapor_density

    canopy_top_eddy_diffusivity = _ktheory_canopy_top_eddy_diffusivity(buffers, boundary_layer_model, friction_velocity)

    heat_top_conductance = _ktheory_conductances!(heat_conductance_between_layers, relative_eddy_diffusivity,
        canopy_top_eddy_diffusivity, layer_spacing, air_volumetric_heat_capacity)
    heat_ground_conductance = convert(eltype(heat_conductance_between_layers), air_volumetric_heat_capacity / ground_resistance.heat)
    _ktheory_solve!(air_temperature, heat_layer_total_conductance, heat_conductance_between_layers, heat_known_flux, heat_top_conductance, heat_ground_conductance,
        canopy_top_air_temperature, ground_temperature, layer -> sensible_heat_flux[layer])

    # A diffusivity over a spacing is already a velocity, so vapor needs no scaling
    vapor_top_conductance = _ktheory_conductances!(vapor_conductance_between_layers, relative_eddy_diffusivity,
        canopy_top_eddy_diffusivity, layer_spacing, 1.0)
    vapor_ground_conductance = convert(eltype(vapor_conductance_between_layers), 1 / ground_resistance.vapor)
    _ktheory_solve!(vapor_density, vapor_layer_total_conductance, vapor_conductance_between_layers, vapor_known_flux, vapor_top_conductance, vapor_ground_conductance,
        canopy_top_vapor_density, ground_vapor_density, layer -> evaporation_mass_flux[layer])

    return nothing
end
