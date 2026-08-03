# Shared leaf-level heat/mass transfer physics, used by every
# AbstractLeafTemperatureSolver variant. Built on HeatExchange.jl's
# convection()/evaporation(): shape-specific Nusselt correlations and
# water-potential-driven surface humidity (Kelvin equation).

# The oblate spheroid's equatorial radius is the boundary-layer length scale
# directly -- no extra empirical scale factor on top of it.
const LEAF_CHARACTERISTIC_DIMENSION = ScaledDimension(1.0, :b_semi_minor_skin)

# Flattening ratio (polar axis / equatorial radius) for the leaf-as-oblate-
# spheroid approximation. Free/tunable -- only the equatorial radius (set
# below to sqrt(leaf_length*leaf_width)) affects convection; this only
# affects the geometry's own (unused downstream) thinness.
const LEAF_FLATTENING = 0.01

# Campbell (1990): leaf-angle-dependent reduction of the aerodynamic width,
# same `1.774*(y+1.182)^(-0.733)` constants as `ellipsoidal_extinction_coefficient`
# but evaluated at y = 1/leaf_angle_distribution_parameter. x==1 (spherical)
# gives leaf_angle_width_factor(1) = 1/(1+1.774*2.182^-0.733) ≈ 0.5.
function leaf_angle_width_factor(leaf_angle_distribution_parameter)
    inverse_x = 1.0 / leaf_angle_distribution_parameter
    return 1.0 / (inverse_x + 1.774 * (inverse_x + 1.182)^(-0.733))
end

"""
    leaf_body(leaf_length, leaf_width, leaf_angle_distribution_parameter)

A `BiophysicalGeometry.Body` representing a leaf as a flattened (oblate)
ellipsoid, sized so its equatorial radius (`b_semi_minor_skin` ==
`c_semi_minor_skin`) equals `sqrt(leaf_length * leaf_width)` scaled by
[`leaf_angle_width_factor`](@ref).
"""
function leaf_body(leaf_length, leaf_width, leaf_angle_distribution_parameter)
    equatorial_radius = sqrt(leaf_length * leaf_width) * leaf_angle_width_factor(leaf_angle_distribution_parameter)
    reference_density = 1.0u"kg/m^3"
    volume = equatorial_radius^3 * (4π / 3) * LEAF_FLATTENING
    shape = Ellipsoid(volume * reference_density, reference_density, LEAF_FLATTENING, LEAF_FLATTENING)
    return Body(shape, Naked())
end

"""
    AbstractLeafConvectionModel

Selects the Nusselt correlation `leaf_convection` uses.
"""
abstract type AbstractLeafConvectionModel end

"""
    ElaborateLeafConvection()

Combined free+forced convection, shape-dispatched via the leaf's
oblate-ellipsoid geometry (`HeatExchange.ShapeCorrelation`). Default.
"""
struct ElaborateLeafConvection <: AbstractLeafConvectionModel end

"""
    SimpleLeafConvection()

Forced convection only, `Nu = 0.6·√Re`, no shape dependence
(`HeatExchange.SimpleForcedCorrelation`). Cheaper and less complete than
[`ElaborateLeafConvection`](@ref), not more accurate.
"""
struct SimpleLeafConvection <: AbstractLeafConvectionModel end

_convection_correlation(::ElaborateLeafConvection) = HeatExchange.ShapeCorrelation()
_convection_correlation(::SimpleLeafConvection) = HeatExchange.SimpleForcedCorrelation()

# Boundary-layer heat/mass transfer coefficients for one leaf layer.
leaf_convection(convection_model, body, leaf_area, air_temperature, leaf_temperature, wind_speed, atmospheric_pressure) =
    HeatExchange.convection(; body, area=leaf_area, air_temperature, surface_temperature=leaf_temperature,
        wind_speed, atmospheric_pressure, fluid=Air(), characteristic_dimension_formula=LEAF_CHARACTERISTIC_DIMENSION,
        correlation=_convection_correlation(convection_model))

"""
    leaf_heat_balance(leaf_temperature, absorbed_radiation, air_temperature, relative_humidity,
                       wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance,
                       leaf_water_potential, body, leaf_area, convection_model)

Net leaf energy balance (W, positive = leaf gaining energy): absorbed
radiation minus emitted longwave, convection, and transpiration.
"""
function leaf_heat_balance(leaf_temperature, absorbed_radiation, air_temperature, relative_humidity,
    wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential,
    body, leaf_area, convection_model=ElaborateLeafConvection(),
)
    conv = leaf_convection(convection_model, body, leaf_area, air_temperature, leaf_temperature, wind_speed, atmospheric_pressure)
    atmos = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
    evap = HeatExchange.evaporation(stomatal_conductance, conv.mass_transfer_coefficient, atmos, leaf_area,
        leaf_temperature, air_temperature; water_potential=leaf_water_potential)
    emitted_longwave = leaf_emissivity * σ * leaf_temperature^4 * leaf_area
    net = absorbed_radiation * leaf_area - emitted_longwave - conv.convection_flow - evap.evaporation_heat_flow
    return (; net, conv, evap, emitted_longwave)
end
