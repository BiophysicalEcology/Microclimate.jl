abstract type AbstractSoilHydraulicsModel end

"""
    allocate_soil_water_balance(soil_hydraulic_model, num_layers)

Allocate per-layer scratch buffers for the variant's water-balance solver.
"""
function allocate_soil_water_balance end

"""
    soil_water_balance(soil_hydraulic_model; num_layers, kwargs...)

Convenience wrapper that allocates scratch via `allocate_soil_water_balance`
and forwards to `soil_water_balance!`.
"""
function soil_water_balance end

"""
    soil_water_balance!(buffers, soil_hydraulic_model; soil_profile, depths, pool, evaporation_potential, local_relative_humidity, niter_moist, soil_moisture, moisture_timestep, moisture_tolerance, moisture_max_iterations, max_surface_pool, T0, vapour_pressure_equation, canopy_transpiration_potential, canopy_leaf_area_index, environment_instant)

One-hour water-balance step: `niter_moist` sub-iterations of the variant's
per-timestep infiltration solver, driven by `evaporation_potential` (the
ground-surface evaporative demand after [`ground_condensation_step!`](@ref)
has already shielded it with any standing dew/frost) and
`local_relative_humidity` (also from that call).
"""
function soil_water_balance! end
