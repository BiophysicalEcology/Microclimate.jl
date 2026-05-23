"""
    NoSnow()

Snow model that disables snow physics entirely. All dispatches are no-ops or
return sentinel values so callers can run the rest of the microclimate solver
without paying any snow cost (and without runtime branches on `nothing`).
"""
struct NoSnow <: AbstractSnowModel end

n_snow_nodes(::NoSnow) = 0

init_snow_temperatures(::NoSnow, _) = nothing

initial_snow_state(::NoSnow) = nothing
initial_snow_state(::NoSnow, ::Any, ::Any) = nothing

allocate_snow_scratch(::NoSnow, args...) = nothing
reset_snow_scratch!(::NoSnow, _) = nothing

activate_snow_nodes!(::NoSnow, args...) = nothing

compute_effective_depths!(::NoSnow, _, _, _) = nothing

snow_surface_overrides(::NoSnow, _, _, _) = (;
    albedo=nothing, emissivity=nothing, soil_wetness=nothing, shade=nothing
)

snow_phase_transition(::NoSnow, args...) = (nothing, nothing, nothing)

update_snow!(::NoSnow, args...; kw...) = (nothing, 0.0u"cm", 0.0u"cm")

adjust_snow_near_nodes!(::NoSnow, args...) = nothing

@inline _n_snow(::NoSnow) = 0

@inline _recompute_soil_snow_properties!(::NoSnow, p, st, ap, vpe) =
    _recompute_soil_snow_properties_no_snow!(p, st, ap, vpe)
