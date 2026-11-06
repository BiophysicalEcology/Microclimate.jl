"""
    KTheoryAirProfile()

Air temperature and vapor density through the canopy, from steady diffusion between
neighbouring layers: `d/dz(K dT/dz) = -S`, with each layer's leaves as the source `S`.

Layer 1 is fixed at the canopy-top value, and heat and vapor diffuse from there down to the
ground. Unlike [`RaupachLTheoryAirProfile`](@ref), each layer only exchanges directly with
its neighbours.

The eddy diffusivity `K` follows the shape of the wind profile
(`wind_attenuation_profile`), scaled to `κ u* (canopy_height - displacement_height)`
at the canopy top. Heat and vapor share that diffusivity, and each is a tridiagonal system
factorized in place in preallocated buffers.

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
    heat_diagonal, heat_off_diagonal = zeros(free_layer_count), zeros(free_layer_count - 1)
    vapor_diagonal, vapor_off_diagonal = zeros(free_layer_count), zeros(free_layer_count - 1)

    return (;
        canopy_height, displacement_height,
        layer_spacing = layer_thickness, relative_eddy_diffusivity,
        heat_diagonal, heat_off_diagonal, heat_right_hand_side = zeros(free_layer_count),
        vapor_diagonal, vapor_off_diagonal, vapor_right_hand_side = zeros(free_layer_count),
        heat_matrix = SymTridiagonal(heat_diagonal, heat_off_diagonal),
        vapor_matrix = SymTridiagonal(vapor_diagonal, vapor_off_diagonal),
        air_temperature = zeros(typeof(0.0u"K"), n_layers),
        current_air_temperature = zeros(typeof(0.0u"K"), n_layers),
        vapor_density = zeros(typeof(0.0u"kg/m^3"), n_layers),
        current_vapor_density = zeros(typeof(0.0u"kg/m^3"), n_layers),
        relative_humidity = zeros(n_layers),
        ground_heat_conductance = Ref(0.0u"W/m^2/K"),
        ground_vapor_conductance = Ref(0.0u"m/s"),
    )
end

# Writes the conductance between each pair of free layers into `off_diagonal`, and returns the
# conductances from the fixed layer 1 to the first free layer and from the last free layer to
# the ground. `scale` turns the diffusivity over the spacing into the flux's units.
function _ktheory_conductances!(off_diagonal, relative_eddy_diffusivity, canopy_top_eddy_diffusivity, layer_spacing, scale, unit)
    free_layer_count = length(off_diagonal) + 1
    top_conductance = ustrip(unit, scale * canopy_top_eddy_diffusivity * 0.5 *
        (relative_eddy_diffusivity[1] + relative_eddy_diffusivity[2]) / layer_spacing[1])
    ground_conductance = ustrip(unit, scale * canopy_top_eddy_diffusivity * relative_eddy_diffusivity[free_layer_count + 1] /
        layer_spacing[free_layer_count + 1])
    for i in 1:(free_layer_count - 1)
        boundary_eddy_diffusivity = canopy_top_eddy_diffusivity * 0.5 *
            (relative_eddy_diffusivity[i + 1] + relative_eddy_diffusivity[i + 2])
        off_diagonal[i] = ustrip(unit, scale * boundary_eddy_diffusivity / layer_spacing[i + 1])
    end
    return top_conductance, ground_conductance
end

# Fills in and solves one flux balance on the free layers, leaving the answer in `right_hand_side`:
#   g₋(C₋ - Cᵢ) - g₊(Cᵢ - C₊) + Sᵢ = 0
# `source(layer)` is Sᵢ, and layer 1 and the ground give the known end values.
function _ktheory_solve!(matrix, diagonal, off_diagonal, right_hand_side, top_conductance, ground_conductance, top_value, ground_value, source)
    free_layer_count = length(diagonal)
    for i in 1:free_layer_count
        conductance_above = i == 1 ? top_conductance : off_diagonal[i - 1]
        conductance_below = i == free_layer_count ? ground_conductance : off_diagonal[i]
        diagonal[i] = -(conductance_above + conductance_below)
        right_hand_side[i] = -source(i + 1)
    end
    right_hand_side[1] -= top_conductance * top_value
    right_hand_side[free_layer_count] -= ground_conductance * ground_value
    ldiv!(ldlt!(matrix), right_hand_side)
    return right_hand_side
end

function canopy_heat_vapor_transport!(buffers, ::KTheoryAirProfile, boundary_layer_model,
    forcing::CanopyAirProfileForcing; site, vapour_pressure_equation=GoffGratch(),
)
    (; friction_velocity, canopy_top_air_temperature, canopy_top_vapor_density, ground_temperature,
       ground_relative_humidity, sensible_heat_source, evaporation_mass_flow, atmospheric_pressure) = forcing
    (; canopy_height, displacement_height, layer_spacing, relative_eddy_diffusivity,
       heat_diagonal, heat_off_diagonal, heat_right_hand_side,
       vapor_diagonal, vapor_off_diagonal, vapor_right_hand_side,
       heat_matrix, vapor_matrix, air_temperature, vapor_density,
       ground_heat_conductance, ground_vapor_conductance) = buffers

    layer_count = length(air_temperature)
    length(sensible_heat_source) == layer_count || throw(ArgumentError("sensible_heat_source must have one entry per canopy layer"))
    length(evaporation_mass_flow) == layer_count || throw(ArgumentError("evaporation_mass_flow must have one entry per canopy layer"))

    air_temperature[1] = canopy_top_air_temperature
    vapor_density[1] = canopy_top_vapor_density

    # The floor stops the diffusivity reaching zero when displacement reaches the canopy top
    canopy_top_eddy_diffusivity = boundary_layer_model.karman_constant * friction_velocity * max(canopy_height - displacement_height, 1.0e-3u"m")

    # Heat conductances are in W/m²/K
    heat_top_conductance, heat_ground_conductance = _ktheory_conductances!(heat_off_diagonal, relative_eddy_diffusivity,
        canopy_top_eddy_diffusivity, layer_spacing, calc_ρ_cp(canopy_top_air_temperature), u"W/m^2/K")
    ground_heat_conductance[] = heat_ground_conductance * u"W/m^2/K"
    _ktheory_solve!(heat_matrix, heat_diagonal, heat_off_diagonal, heat_right_hand_side, heat_top_conductance, heat_ground_conductance,
        ustrip(u"K", canopy_top_air_temperature), ustrip(u"K", ground_temperature), layer -> ustrip(u"W/m^2", sensible_heat_source[layer]))
    air_temperature[2:end] .= heat_right_hand_side .* u"K"

    # A diffusivity over a spacing is already a velocity, so vapor conductances are in m/s
    vapor_top_conductance, vapor_ground_conductance = _ktheory_conductances!(vapor_off_diagonal, relative_eddy_diffusivity,
        canopy_top_eddy_diffusivity, layer_spacing, 1.0, u"m/s")
    ground_vapor_conductance[] = vapor_ground_conductance * u"m/s"
    ground_vapor_density = wet_air_properties(ground_temperature, ground_relative_humidity, atmospheric_pressure; vapour_pressure_equation).vapour_density
    _ktheory_solve!(vapor_matrix, vapor_diagonal, vapor_off_diagonal, vapor_right_hand_side, vapor_top_conductance, vapor_ground_conductance,
        ustrip(u"kg/m^3", canopy_top_vapor_density), ustrip(u"kg/m^3", ground_vapor_density),
        layer -> ustrip(u"kg/m^2/s", evaporation_mass_flow[layer] / 1.0u"m^2"))
    vapor_density[2:end] .= vapor_right_hand_side .* u"kg/m^3"

    return nothing
end
