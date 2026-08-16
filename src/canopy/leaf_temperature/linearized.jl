"""
    LinearizedLeafTemperature()

One frozen-coefficient quasi-Newton step of `leaf_heat_balance` around
`leaf_temperature_guess`: longwave uses the exact `4σT³` slope, but
convection/evaporation use their coefficients as evaluated at the guess,
not re-differentiated (so their own T-dependence, e.g. free-convection h(T),
is omitted from the slope). Cheap, no bracket — intended to be called from
an outer loop supplying an improving guess each time; convergence to the
true root doesn't depend on this slope being exact, only its sign.
"""
struct LinearizedLeafTemperature <: AbstractLeafTemperatureSolver end

function leaf_temperature(::LinearizedLeafTemperature, absorbed_radiation, air_temperature, relative_humidity,
    wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body;
    leaf_area=1.0u"m^2", leaf_temperature_guess=air_temperature,
    convection_model=ElaborateLeafConvection(),
)
    out = leaf_heat_balance(leaf_temperature_guess, absorbed_radiation, air_temperature, relative_humidity,
        wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body, leaf_area,
        convection_model)

    radiative_heat_transfer_coefficient = 4.0 * out.emitted_longwave / leaf_temperature_guess
    convective_heat_transfer_coefficient = out.conv.heat_transfer_coefficient.combined * leaf_area

    δ = 0.1u"K"  # finite-difference step, free/tunable, uncited
    atmos = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
    evap_perturbed = _clamped_leaf_evaporation(stomatal_conductance, out.conv.mass_transfer_coefficient, atmos,
        leaf_area, leaf_temperature_guess + δ, air_temperature, leaf_water_potential)
    evaporative_heat_transfer_coefficient = (evap_perturbed.evaporation_heat_flow - out.evap.evaporation_heat_flow) / δ

    return leaf_temperature_guess + out.net / (
        radiative_heat_transfer_coefficient + convective_heat_transfer_coefficient + evaporative_heat_transfer_coefficient
    )
end
