"""
    BulkTransferCondensation()

Ground dew/frost from the model's own bulk-transfer surface evaporation flux
(assumed-saturated surface): condenses whenever the surface's saturated
vapour density falls below the actual air vapour density. Self-consistent
with the model's own boundary-layer conductance, but has no
energy-availability gating (see [`GarrattSegalCondensation`](@ref)), so it
does not distinguish nocturnal radiative dew from any other
surface-cooler-than-air event and can produce implausibly large fluxes when
the boundary-layer conductance is large.
"""
struct BulkTransferCondensation <: AbstractCondensationModel end

function condensation_energy_flux(::BulkTransferCondensation;
    evaporation_model, surface_temperature, air_temperature, relative_humidity,
    atmospheric_pressure, mass_transfer_coefficient, vapour_pressure_equation=GoffGratch(), kw...,
)
    Q_evaporation, _ = evaporation(evaporation_model;
        surface_temperature, air_temperature, relative_humidity,
        surface_relative_humidity=1.0, mass_transfer_coefficient, atmospheric_pressure,
        soil_wetness=1.0, saturated=true, vapour_pressure_equation,
    )
    return Q_evaporation
end
