abstract type AbstractSoilHydraulicsModel end

"""
    allocate_soil_water_balance(soil_hydraulic_model, num_layers)

Allocate per-layer scratch buffers for the variant's water-balance solver.
"""
function allocate_soil_water_balance end

"""
    soil_water_balance(soil_hydraulic_model, moisture_mode; num_layers, kwargs...)

Convenience wrapper that allocates scratch via `allocate_soil_water_balance`
and forwards to `soil_water_balance!`.
"""
function soil_water_balance end

"""
    soil_water_balance!(buffers, soil_hydraulic_model, moisture_mode::DynamicSoilMoisture;
        soil_profile, depths, site, boundary_layer_model, environment_instant, T0, pool, soil_moisture,
        evaporation_model, vapour_pressure_equation, snow_present)

One-hour water-balance step: atmospheric profile, evaporation, surface pool
dynamics, and sub-steps of the variant's solver. Solver settings and the
current soil wetness come from `moisture_mode`.
"""
function soil_water_balance! end
