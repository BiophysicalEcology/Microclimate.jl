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
    soil_water_balance!(buffers, soil_hydraulic_model; soil_profile, depths, site, boundary_layer_model, environment_instant, T0, pool, niter_moist, soil_wetness, soil_moisture, moisture_timestep, moisture_tolerance, moisture_max_iterations, max_surface_pool, evaporation_model, frozen_water_content, vapour_pressure_equation, snow_present)

One-hour water-balance step: atmospheric profile, evaporation, surface pool
dynamics, and `niter_moist` sub-iterations of the variant's per-timestep solver.
"""
function soil_water_balance! end

"""
    NoIce()

Default `frozen_water_content`: every layer reports 0.
"""
struct NoIce end
Base.getindex(::NoIce, ::Int) = 0.0

"""
    ice_impeded_conductivity(hydraulic_conductivity, ice_content, porosity)

Bloomsburg & Wang (1969) ice-blocking of conductivity, floored at
[`ICE_CONDUCTIVITY_FLOOR_FACTOR`](@ref) (numerical regularization only).
"""
@inline function ice_impeded_conductivity(hydraulic_conductivity, ice_content, porosity)
    ice_content <= zero(ice_content) && return hydraulic_conductivity
    available = porosity - ice_content
    threshold = oftype(available, ICE_IMPEDANCE_MIN_POROSITY)
    floor_fraction = oftype(available, ICE_CONDUCTIVITY_FLOOR_FACTOR)
    available <= threshold && return hydraulic_conductivity * floor_fraction
    return hydraulic_conductivity * max(floor_fraction, (available - threshold) / (porosity - threshold))
end

"""
    ice_free_capacity(porosity, ice_content)

Ice-free porosity remaining; not a storage ceiling (`soil_moisture` already
includes ice) — used only to test whether conductivity has collapsed.
"""
@inline ice_free_capacity(porosity, ice_content) = max(zero(porosity), porosity - ice_content)
