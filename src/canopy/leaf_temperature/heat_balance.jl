# Shared leaf-level heat/mass transfer physics, used by every
# AbstractLeafTemperatureSolver variant. Built on HeatExchange.jl's
# convection()/evaporation(): shape-specific Nusselt correlations (Gates 1980;
# Bird, Stewart & Lightfoot 1960) and water-potential-driven surface humidity
# (Kelvin equation).

# Campbell & Norman (1998): leaf boundary-layer characteristic dimension is
# 0.7 x leaf width.
const LEAF_CHARACTERISTIC_DIMENSION = ScaledDimension(0.7, :width_skin)

"""
    leaf_body(leaf_length, leaf_width)

A `BiophysicalGeometry.Body` representing a leaf as a thin plate, sized so
its `width_skin` equals `sqrt(leaf_length * leaf_width)`. Mass, density, and
thinness only affect `width_skin` through this one constraint (fixed here to
convenient values); nothing else about them is used downstream.
"""
function leaf_body(leaf_length, leaf_width)
    characteristic_width = sqrt(leaf_length * leaf_width)
    reference_density = 1.0u"kg/m^3"
    shape = Plate(characteristic_width^3 * reference_density, reference_density, 1.0, 1.0)
    return Body(shape, Naked())
end

# Boundary-layer heat/mass transfer coefficients for one leaf layer.
leaf_convection(body, leaf_area, air_temperature, leaf_temperature, wind_speed, atmospheric_pressure) =
    HeatExchange.convection(; body, area=leaf_area, air_temperature, surface_temperature=leaf_temperature,
        wind_speed, atmospheric_pressure, fluid=Air(), characteristic_dimension_formula=LEAF_CHARACTERISTIC_DIMENSION)

"""
    leaf_heat_balance(leaf_temperature, absorbed_radiation, air_temperature, relative_humidity,
                       wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance,
                       leaf_water_potential, body, leaf_area)

Net leaf energy balance (W, positive = leaf gaining energy): absorbed
radiation minus emitted longwave, convection, and transpiration.
"""
function leaf_heat_balance(leaf_temperature, absorbed_radiation, air_temperature, relative_humidity,
    wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential,
    body, leaf_area,
)
    conv = leaf_convection(body, leaf_area, air_temperature, leaf_temperature, wind_speed, atmospheric_pressure)
    atmos = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
    evap = HeatExchange.evaporation(stomatal_conductance, conv.mass_transfer_coefficient, atmos, leaf_area,
        leaf_temperature, air_temperature; water_potential=leaf_water_potential)
    emitted_longwave = leaf_emissivity * σ * leaf_temperature^4 * leaf_area
    net = absorbed_radiation * leaf_area - emitted_longwave - conv.convection_flow - evap.evaporation_heat_flow
    return (; net, conv, evap, emitted_longwave)
end
