"""
    AbstractRaupachFormulation

How [`RaupachLTheoryAirProfile`](@ref) calculates its far and near fields:
[`ExactRaupachFormulation`](@ref) or [`BulkRaupachFormulation`](@ref).
"""
abstract type AbstractRaupachFormulation end

"""
    ExactRaupachFormulation()

Follows Raupach (1989): eq. 19a for the far field and eq. 37 for the near field. A layer's near
field from its own source averages the kernel over the layer, since the kernel is infinite at
zero distance.
"""
struct ExactRaupachFormulation <: AbstractRaupachFormulation end

"""
    BulkRaupachFormulation()

Cheaper: the far field is each layer's flux times the resistance to the top; the near field
evaluates the kernel once, skips a layer's own source, and is scaled by the layer count.
"""
struct BulkRaupachFormulation <: AbstractRaupachFormulation end

"""
    RaupachLTheoryAirProfile(; ground_velocity_deviation_factor=0.25, canopy_top_velocity_deviation_factor=1.25,
                              formulation=ExactRaupachFormulation(), max_lagrangian_timescale=200.0u"s")

Air temperature and vapor density through the canopy, from Raupach's (1989) Lagrangian
transport model. Each layer's value is the canopy-top value plus two parts:

- a far field, from all the heat and vapor sources mixed together
- a near field, a correction close to each source, where mixing hasn't happened yet

Every layer depends on every layer's source.

`ground_velocity_deviation_factor` (`a0`) and `canopy_top_velocity_deviation_factor` (`a1`) set the
standard deviation of vertical velocity `σ_w`, from `a0·u*` at the ground to `a1·u*` at the top
(`u*` friction velocity, `h` canopy height): `σ_w(z) = (a1+a0)/2·u* + (a1-a0)/2·u*·cos(π(1 - z/h))`.

The ground's resistance is the lowest layer's turbulent resistance plus a thin sublayer at
the ground ([`sublayer_stanton`](@ref), with a friction velocity from the wind in that layer).

`max_lagrangian_timescale` caps the Lagrangian timescale, which grows in calm and unstable air.

`formulation` is an [`AbstractRaupachFormulation`](@ref), which sets how both fields are calculated.

# References
- Raupach, M. R. (1989). A practical Lagrangian method for relating scalar
  concentrations to source distributions in vegetation canopies. *Quarterly
  Journal of the Royal Meteorological Society*, 115(487), 609-632.
"""
@kwdef struct RaupachLTheoryAirProfile{GVSF,CVSF,F<:AbstractRaupachFormulation,MTL} <: AbstractCanopyAirProfileModel
    ground_velocity_deviation_factor::GVSF = 0.25
    canopy_top_velocity_deviation_factor::CVSF = 1.25
    formulation::F = ExactRaupachFormulation()
    max_lagrangian_timescale::MTL = 200.0u"s"
end

function allocate_air_profile(model::RaupachLTheoryAirProfile, canopy_height, displacement_height, plant_area_index, heights, n_layers, boundary_layer_model)
    (; layer_heights, layer_thickness) = canopy_layer_heights(heights, canopy_height, n_layers)
    air_temperature = zeros(typeof(0.0u"K"), n_layers)
    vapor_density = zeros(typeof(0.0u"kg/m^3"), n_layers)
    return (;
        canopy_height,
        displacement_height,
        layer_heights,
        layer_thickness,
        inverse_near_field_length = zeros(typeof(1.0 / (1.0u"m/s" * 1.0u"s")), n_layers),
        layer_resistance = zeros(typeof(0.0u"s/m"), n_layers),
        far_field_heat = zeros(typeof(0.0u"J/m^3"), n_layers),
        far_field_vapor = zeros(typeof(0.0u"kg/m^3"), n_layers),
        resistance_from_top = zeros(typeof(0.0u"s/m"), n_layers),
        resistance_to_ground = zeros(typeof(0.0u"s/m"), n_layers),
        cumulative_sensible_below = zeros(typeof(0.0u"W/m^2"), n_layers),
        sensible_near_field_weight = zeros(typeof(0.0u"W/m^2" / 1.0u"m/s"), n_layers),
        cumulative_latent_below = zeros(typeof(0.0u"kg/m^2/s"), n_layers),
        latent_near_field_weight = zeros(typeof(1.0u"kg/m^2/s" / 1.0u"m/s"), n_layers),
        air_temperature,
        current_air_temperature = similar(air_temperature),
        vapor_density,
        current_vapor_density = similar(vapor_density),
        relative_humidity = zeros(n_layers),
    )
end

# Raupach (1989), eq. 34, with ζ the source distance over σw·T_L:
#   kₙ(ζ) ≈ c₁ log(1 - exp(-|ζ|)) + c₂ exp(-|ζ|)
# The constants satisfy ∫₀^∞ kₙ(ζ)dζ = 1/2.
const _RAUPACH_KERNEL_LOG_COEFFICIENT = -1 / sqrt(2π)
const _RAUPACH_KERNEL_EXPONENTIAL_COEFFICIENT = 0.5 - π^2 / (6sqrt(2π))

function _raupach_kernel(scaled_distance)
    ζ, c₁, c₂ = scaled_distance, _RAUPACH_KERNEL_LOG_COEFFICIENT, _RAUPACH_KERNEL_EXPONENTIAL_COEFFICIENT
    decay = exp(-abs(ζ))
    return c₁ * log1p(-decay) + c₂ * decay
end

# Raupach (1989), eq. 37:
#   Cₙ(z) = ∫ S(z₀)/σw(z₀) [kₙ((z-z₀)/λ) + kₙ((z+z₀)/λ)] dz₀
# The two terms are the direct path and its reflection from the ground. Both use products
# of the same two exponentials, so the pair takes one exp and one log each.
function _raupach_kernel_weight(evaluation_height, source_layer, layer_heights, inverse_near_field_length)
    z, z₀, inverse_λ = evaluation_height, layer_heights[source_layer], inverse_near_field_length[source_layer]
    c₁, c₂ = _RAUPACH_KERNEL_LOG_COEFFICIENT, _RAUPACH_KERNEL_EXPONENTIAL_COEFFICIENT
    evaluation_decay = exp(-z * inverse_λ)
    source_decay = exp(-z₀ * inverse_λ)
    direct = z > z₀ ? evaluation_decay / source_decay : source_decay / evaluation_decay
    reflected = evaluation_decay * source_decay
    return c₁ * log((1 - direct) * (1 - reflected)) + c₂ * (direct + reflected)
end

"""
    _raupach_kernel_integral(scaled_distance)

Exact integral of the Raupach kernel from 0 to `scaled_distance` (≥ 0), via the dilogarithm `reli2`.
"""
function _raupach_kernel_integral(scaled_distance)
    ζ, c₁, c₂ = scaled_distance, _RAUPACH_KERNEL_LOG_COEFFICIENT, _RAUPACH_KERNEL_EXPONENTIAL_COEFFICIENT
    decay = exp(-ζ)
    return c₁ * (reli2(decay) - π^2 / 6) + c₂ * (1 - decay)
end

# The kernel is infinite at zero distance, so a source's effect on a point within its own
# layer is the kernel's average over the layer. The kernel is even, so a point inside the
# layer adds the integrals to each edge, and a point outside takes their difference.
function _raupach_self_kernel_weight(evaluation_height, source_layer, layer_heights, layer_thickness, inverse_near_field_length)
    z, inverse_λ = evaluation_height, inverse_near_field_length[source_layer]
    bottom = layer_heights[source_layer] - layer_thickness[source_layer] / 2
    top = layer_heights[source_layer] + layer_thickness[source_layer] / 2

    to_bottom = _raupach_kernel_integral(abs((z - bottom) * inverse_λ))
    to_top = _raupach_kernel_integral(abs((top - z) * inverse_λ))
    direct = bottom <= z <= top ? to_bottom + to_top : abs(to_bottom - to_top)
    reflected = _raupach_kernel_integral((z + top) * inverse_λ) - _raupach_kernel_integral((z + bottom) * inverse_λ)
    return (direct + reflected) / (layer_thickness[source_layer] * inverse_λ)
end

# The bulk formulation evaluates the kernel once, at the direct distance, and multiplies by the
# sum of the direct and reflected distances
function _raupach_bulk_kernel_weight(evaluation_height, source_layer, layer_heights, inverse_near_field_length)
    z, z₀, inverse_λ = evaluation_height, layer_heights[source_layer], inverse_near_field_length[source_layer]
    direct_distance = (z - z₀) * inverse_λ
    reflected_distance = (z + z₀) * inverse_λ
    return _raupach_kernel(direct_distance) * (direct_distance + reflected_distance)
end

# The near-field weight at the canopy top of source `layer`. Layer 1 is always the self
# term, as in _raupach_pair_kernel_weight, because its node need not sit exactly at the
# canopy top. Using height here would treat it differently in the two places.
function _raupach_top_kernel_weight(::ExactRaupachFormulation, layer, canopy_height, layer_heights, layer_thickness, inverse_near_field_length)
    return layer == 1 ?
        _raupach_self_kernel_weight(canopy_height, layer, layer_heights, layer_thickness, inverse_near_field_length) :
        _raupach_kernel_weight(canopy_height, layer, layer_heights, inverse_near_field_length)
end
function _raupach_top_kernel_weight(::BulkRaupachFormulation, layer, canopy_height, layer_heights, layer_thickness, inverse_near_field_length)
    return layer == 1 ? 0.0 : _raupach_bulk_kernel_weight(canopy_height, layer, layer_heights, inverse_near_field_length)
end

# The near-field weight at `layer` of a source in `source_layer`
function _raupach_pair_kernel_weight(::ExactRaupachFormulation, layer, source_layer, layer_heights, layer_thickness, inverse_near_field_length)
    return source_layer == layer ?
        _raupach_self_kernel_weight(layer_heights[layer], source_layer, layer_heights, layer_thickness, inverse_near_field_length) :
        _raupach_kernel_weight(layer_heights[layer], source_layer, layer_heights, inverse_near_field_length)
end
function _raupach_pair_kernel_weight(::BulkRaupachFormulation, layer, source_layer, layer_heights, layer_thickness, inverse_near_field_length)
    return source_layer == layer ? 0.0 : _raupach_bulk_kernel_weight(layer_heights[layer], source_layer, layer_heights, inverse_near_field_length)
end

# Near-field heat and vapor at `layer` from every layer's source
function _raupach_near_field(formulation, layer, layer_heights, layer_thickness, inverse_near_field_length,
    sensible_near_field_weight, latent_near_field_weight, near_field_scale,
)
    heat = zero(eltype(sensible_near_field_weight))
    vapor = zero(eltype(latent_near_field_weight))
    for source_layer in eachindex(sensible_near_field_weight)
        kernel_weight = _raupach_pair_kernel_weight(formulation, layer, source_layer, layer_heights, layer_thickness, inverse_near_field_length)
        heat += sensible_near_field_weight[source_layer] * kernel_weight
        vapor += latent_near_field_weight[source_layer] * kernel_weight
    end
    return (; heat=near_field_scale * heat, vapor=near_field_scale * vapor)
end

# Eq. 19a: Cf(z) - Cf(zR) = ∫[z,zR] F(z')/Kf(z') dz'. Each layer's flux over its own
# resistance, summed from the canopy top down.
#
# The flux F(z') is the sources below z' plus one ground flux. The bottom layer's value is
# its offset (the canopy top plus the near field) plus the fluxes across the in-canopy
# resistance, so the ground flux solves directly through the ground resistance and the
# in-canopy resistance in series, without iterating on the bottom layer's own value.
function _raupach_far_field!(::ExactRaupachFormulation, buffers, ground_temperature, ground_vapor_density,
    ground_heat_resistance, ground_vapor_resistance, volumetric_heat_capacity, layer_count;
    bottom_offset_heat, bottom_offset_vapor,
)
    (; far_field_heat, far_field_vapor, cumulative_sensible_below, cumulative_latent_below, layer_resistance) = buffers

    total_resistance = sum(layer_resistance)
    offset_heat = bottom_offset_heat
    offset_vapor = bottom_offset_vapor
    for layer in 1:layer_count
        offset_heat += cumulative_sensible_below[layer] * layer_resistance[layer]
        offset_vapor += cumulative_latent_below[layer] * layer_resistance[layer]
    end
    ground_heat_flux = (volumetric_heat_capacity * ground_temperature - offset_heat) / (ground_heat_resistance + total_resistance)
    ground_vapor_flux = (ground_vapor_density - offset_vapor) / (ground_vapor_resistance + total_resistance)

    far_field_heat[1] = (cumulative_sensible_below[1] + ground_heat_flux) * layer_resistance[1]
    far_field_vapor[1] = (cumulative_latent_below[1] + ground_vapor_flux) * layer_resistance[1]
    for layer in 2:layer_count
        far_field_heat[layer] = far_field_heat[layer - 1] + (cumulative_sensible_below[layer] + ground_heat_flux) * layer_resistance[layer]
        far_field_vapor[layer] = far_field_vapor[layer - 1] + (cumulative_latent_below[layer] + ground_vapor_flux) * layer_resistance[layer]
    end
    return nothing
end

# The bulk formulation's far field: the flux at each layer times the whole resistance from there to the
# top, with the ground flux from that layer's current value. Heat and vapor share the vapor
# resistance to the ground. Heat capacity comes from each layer's current temperature, as
# in micropoint. Fixing it instead made the canopy-top solve unstable.
function _raupach_far_field!(::BulkRaupachFormulation, buffers, ground_temperature, ground_vapor_density,
    _ground_heat_resistance, ground_vapor_resistance, _volumetric_heat_capacity, layer_count; _...,
)
    (; far_field_heat, far_field_vapor, resistance_from_top, resistance_to_ground, cumulative_sensible_below,
       cumulative_latent_below, layer_resistance, current_air_temperature, current_vapor_density) = buffers

    resistance_from_top[1] = layer_resistance[1]
    for layer in 2:layer_count
        resistance_from_top[layer] = resistance_from_top[layer - 1] + layer_resistance[layer]
    end
    resistance_to_ground[layer_count] = ground_vapor_resistance
    for layer in (layer_count - 1):-1:1
        resistance_to_ground[layer] = resistance_to_ground[layer + 1] + layer_resistance[layer]
    end
    for layer in 1:layer_count
        volumetric_heat_capacity = calc_air_volumetric_heat_capacity(current_air_temperature[layer])
        ground_heat_flux = volumetric_heat_capacity * (ground_temperature - current_air_temperature[layer]) / resistance_to_ground[layer]
        ground_vapor_flux = (ground_vapor_density - current_vapor_density[layer]) / resistance_to_ground[layer]
        far_field_heat[layer] = (cumulative_sensible_below[layer] + ground_heat_flux) * resistance_from_top[layer]
        far_field_vapor[layer] = (cumulative_latent_below[layer] + ground_vapor_flux) * resistance_from_top[layer]
    end
    return nothing
end

# The bulk formulation's correction to the near field, which depends only on the layer count
_raupach_near_field_scale(::ExactRaupachFormulation, layer_count) = 1.0
_raupach_near_field_scale(::BulkRaupachFormulation, layer_count) = 1.0 + 0.894 * exp(-0.01386 * layer_count) + 9.82 * exp(-0.15 * layer_count)

# The heat capacity that converts a layer's heat into a temperature: fixed for the exact
# formulation, and from the layer's current temperature for the bulk one (see _raupach_far_field!)
_raupach_layer_heat_capacity(::ExactRaupachFormulation, volumetric_heat_capacity, _current_air_temperature, _layer) = volumetric_heat_capacity
_raupach_layer_heat_capacity(::BulkRaupachFormulation, _volumetric_heat_capacity, current_air_temperature, layer) = calc_air_volumetric_heat_capacity(current_air_temperature[layer])

# The Lagrangian timescale T_L is the same at every height. It makes the diffusivity at the
# canopy top, K(h) = σw(h)²T_L, match the Monin-Obukhov diffusivity above it, κu*(h-d)/Φh:
#   T_L = min(a₂h/u*, max_lagrangian_timescale), with a₂ = κ(1-d/h)/(a₁²Φh)
function _raupach_lagrangian_timescale(model, boundary_layer_model, canopy_height, displacement_height, friction_velocity,
    inverse_obukhov_length,
)
    h, d, u★, κ = canopy_height, displacement_height, friction_velocity, boundary_layer_model.karman_constant
    a₁ = model.canopy_top_velocity_deviation_factor
    canopy_top_depth = max(h - d, 1.0e-3u"m")
    Φh = heat_gradient_stability_factor(boundary_layer_model.stability_functions,
        stability_parameter_at_height(InverseObukhovLength(), canopy_top_depth, inverse_obukhov_length))
    a₂ = κ * (1.0 - d / h) / (a₁^2 * Φh)
    return min(a₂ * h / u★, model.max_lagrangian_timescale)
end

# The standard deviation of vertical velocity at `height`
function _raupach_vertical_velocity_deviation(model, height, canopy_height, friction_velocity)
    a₀, a₁, u★ = model.ground_velocity_deviation_factor, model.canopy_top_velocity_deviation_factor, friction_velocity
    return (a₁ + a₀) * u★ / 2 + (a₁ - a₀) * u★ / 2 * cos(π * (1.0 - height / canopy_height))
end

# The lowest layer's turbulent resistance plus a thin sublayer at the ground, with a friction
# velocity from the log law in that layer. allocate_canopy checks the layer is above the roughness.
function canopy_ground_resistance(buffers, model::RaupachLTheoryAirProfile, boundary_layer_model;
    friction_velocity, inverse_obukhov_length, lowest_wind_speed, lowest_air_temperature, ground_temperature,
    atmospheric_pressure, site,
)
    (; canopy_height, displacement_height, layer_heights, layer_thickness) = buffers
    T_L = _raupach_lagrangian_timescale(model, boundary_layer_model, canopy_height, displacement_height, friction_velocity,
        inverse_obukhov_length)
    σw = _raupach_vertical_velocity_deviation(model, last(layer_heights), canopy_height, friction_velocity)
    lowest_layer_resistance = last(layer_thickness) / (σw^2 * T_L)

    κ = boundary_layer_model.karman_constant
    z₀ = site.roughness_height
    ground_friction_velocity = κ * lowest_wind_speed / log(last(layer_heights) / z₀)
    kinematic_viscosity = dry_air_properties((lowest_air_temperature + ground_temperature) / 2, atmospheric_pressure).kinematic_viscosity
    sublayer_heat_resistance = 1 / (ground_friction_velocity * sublayer_stanton(z₀, ground_friction_velocity, kinematic_viscosity))
    return (; heat=lowest_layer_resistance + sublayer_heat_resistance,
        vapor=lowest_layer_resistance + sublayer_heat_resistance / LEWIS_HEAT_TO_MASS_RATIO)
end

function canopy_heat_vapor_transport!(buffers, model::RaupachLTheoryAirProfile, boundary_layer_model,
    forcing::CanopyAirProfileForcing,
)
    (; friction_velocity, inverse_obukhov_length, canopy_top_air_temperature, canopy_top_vapor_density,
       ground_temperature, ground_vapor_density, ground_resistance, air_volumetric_heat_capacity,
       sensible_heat_flux, evaporation_mass_flux) = forcing
    (; canopy_height, displacement_height, layer_heights, layer_thickness, inverse_near_field_length, layer_resistance,
       far_field_heat, far_field_vapor, cumulative_sensible_below, sensible_near_field_weight, cumulative_latent_below,
       latent_near_field_weight, air_temperature, current_air_temperature, vapor_density) = buffers
    (; formulation) = model

    layer_count = length(air_temperature)
    length(sensible_heat_flux) == layer_count || throw(ArgumentError("sensible_heat_flux must have one entry per canopy layer"))
    length(evaporation_mass_flux) == layer_count || throw(ArgumentError("evaporation_mass_flux must have one entry per canopy layer"))

    T_L = _raupach_lagrangian_timescale(model, boundary_layer_model, canopy_height, displacement_height, friction_velocity,
        inverse_obukhov_length)
    for layer in 1:layer_count
        σw = _raupach_vertical_velocity_deviation(model, layer_heights[layer], canopy_height, friction_velocity)
        inverse_near_field_length[layer] = 1 / (σw * T_L)
        layer_resistance[layer] = layer_thickness[layer] / (σw^2 * T_L)
        sensible_near_field_weight[layer] = sensible_heat_flux[layer] / σw
        latent_near_field_weight[layer] = evaporation_mass_flux[layer] / σw
    end

    cumulative_sensible_below[layer_count] = sensible_heat_flux[layer_count]
    cumulative_latent_below[layer_count] = evaporation_mass_flux[layer_count]
    for layer in (layer_count - 1):-1:1
        cumulative_sensible_below[layer] = cumulative_sensible_below[layer + 1] + sensible_heat_flux[layer]
        cumulative_latent_below[layer] = cumulative_latent_below[layer + 1] + evaporation_mass_flux[layer]
    end

    near_field_scale = _raupach_near_field_scale(formulation, layer_count)

    # The near field at the canopy top, subtracted from every layer's own near field below
    top_near_field_heat = zero(eltype(far_field_heat))
    top_near_field_vapor = zero(canopy_top_vapor_density)
    for source_layer in 1:layer_count
        kernel_weight = _raupach_top_kernel_weight(formulation, source_layer, canopy_height, layer_heights, layer_thickness, inverse_near_field_length)
        top_near_field_heat += sensible_near_field_weight[source_layer] * kernel_weight
        top_near_field_vapor += latent_near_field_weight[source_layer] * kernel_weight
    end
    top_near_field_heat *= near_field_scale
    top_near_field_vapor *= near_field_scale

    bottom_near_field = _raupach_near_field(formulation, layer_count, layer_heights, layer_thickness, inverse_near_field_length,
        sensible_near_field_weight, latent_near_field_weight, near_field_scale)
    _raupach_far_field!(formulation, buffers, ground_temperature, ground_vapor_density, ground_resistance.heat, ground_resistance.vapor,
        air_volumetric_heat_capacity, layer_count;
        bottom_offset_heat=canopy_top_air_temperature * air_volumetric_heat_capacity + bottom_near_field.heat - top_near_field_heat,
        bottom_offset_vapor=canopy_top_vapor_density + bottom_near_field.vapor - top_near_field_vapor)

    for layer in 1:layer_count
        near_field = _raupach_near_field(formulation, layer, layer_heights, layer_thickness, inverse_near_field_length,
            sensible_near_field_weight, latent_near_field_weight, near_field_scale)
        layer_heat_capacity = _raupach_layer_heat_capacity(formulation, air_volumetric_heat_capacity, current_air_temperature, layer)
        air_temperature[layer] = canopy_top_air_temperature + (far_field_heat[layer] + near_field.heat - top_near_field_heat) / layer_heat_capacity
        vapor_density[layer] = canopy_top_vapor_density + far_field_vapor[layer] + near_field.vapor - top_near_field_vapor
    end

    return nothing
end
