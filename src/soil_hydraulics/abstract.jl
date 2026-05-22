abstract type AbstractSoilHydraulicsModel end

"""
    allocate_soil_water_balance(soil_hydraulics, num_layers)

Allocate per-layer scratch buffers for the variant's water-balance solver.
"""
function allocate_soil_water_balance end

"""
    soil_water_balance(soil_hydraulics; num_layers, kwargs...)

Convenience wrapper that allocates scratch via `allocate_soil_water_balance`
and forwards to `soil_water_balance!`.
"""
function soil_water_balance end

"""
    soil_water_balance!(buffers, soil_hydraulics; depths, site, boundary_layer_model, environment_instant, T0, pool, niter_moist, soil_wetness, soil_moisture, moisture_timestep, moisture_tolerance, moisture_max_iterations, max_surface_pool, vapour_pressure_equation, snow_present)

One-hour water-balance step: atmospheric profile, evaporation, surface pool
dynamics, and `niter_moist` sub-iterations of the variant's per-timestep solver.
"""
function soil_water_balance! end
