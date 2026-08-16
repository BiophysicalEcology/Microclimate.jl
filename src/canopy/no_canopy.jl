"""
    NoCanopy()

Default canopy model: no canopy-layer structure. Vegetation influence is a
scalar `environment_daily.shade` fraction consumed directly by the longwave
(`ViewFactorLongwave`) and soil-energy (`SoilHeatTransport1D`) budgets.
Reproduces existing behaviour bit-for-bit.
"""
struct NoCanopy <: AbstractCanopyModel end

allocate_canopy(::NoCanopy, heights, boundary_layer_model) = nothing

reset_canopy_scratch!(::NoCanopy, buffers) = nothing

snapshot_canopy_water!(::NoCanopy, buffers) = nothing

restore_canopy_water!(::NoCanopy, buffers) = nothing

allocate_canopy_inputs(::NoCanopy; kw...) = nothing

initial_ground_overrides(::NoCanopy) = (;
    ground_shortwave_transmission = nothing, ground_incoming_longwave = nothing,
    ground_wind_speed = nothing, ground_air_temperature = nothing,
    ground_air_relative_humidity = nothing, ground_reference_height = nothing,
    ground_heat_conductance = nothing, ground_vapor_conductance = nothing,
)

n_canopy_layers(::NoCanopy, heights) = 0

canopy_leaf_area_index(::NoCanopy) = nothing

# Bare-ground origin -- temperature_anchor_height=nothing keeps the default z0-anchored behavior.
canopy_profile_origin(::NoCanopy, buffers, site, boundary_layer_model) =
    (; displacement_height = 0.0u"m", roughness_length = site.roughness_height,
        thermal_roughness_model = boundary_layer_model.thermal_roughness_model,
        temperature_anchor_height = nothing)

# Ground/snow surface temperature for solve_air!'s output.profile.
profile_surface_temperature(::NoCanopy, output, i, ground_surface_temperature) = ground_surface_temperature

# Unused (temperature_anchor_height is nothing too) -- default constant-vapour-pressure behavior applies.
profile_surface_humidity(::NoCanopy, output, i) = nothing
