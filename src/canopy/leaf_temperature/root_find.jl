"""
    RootFindLeafTemperature()

Root-find (`HeatExchange.zbrent`) of `leaf_heat_balance` each call. More
expensive than [`LinearizedLeafTemperature`](@ref), but doesn't depend on
the quality of a supplied guess. "Root" means a root of the (evaporation-
clamped) balance within `bracket`, not necessarily unique or the true
physical equilibrium if that lies outside `bracket`.
"""
struct RootFindLeafTemperature <: AbstractLeafTemperatureSolver end

function leaf_temperature(::RootFindLeafTemperature, absorbed_radiation, air_temperature, relative_humidity,
    wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body;
    leaf_area=1.0u"m^2", aerodynamic_area=leaf_area, leaf_temperature_guess=air_temperature,
    convection_model=ElaborateLeafConvection(),
    # Free/tunable, uncited; same asymmetric bracket energy_balance.jl's clamp uses.
    bracket=(air_temperature - 30.0u"K", air_temperature + 40.0u"K"),
)
    residual(T) = ustrip(u"W", leaf_heat_balance(T * u"K", absorbed_radiation, air_temperature, relative_humidity,
        wind_speed, atmospheric_pressure, leaf_emissivity, stomatal_conductance, leaf_water_potential, body, leaf_area,
        convection_model; aerodynamic_area).net)
    lo, hi = ustrip(u"K", bracket[1]), ustrip(u"K", bracket[2])
    q_lo, q_hi = residual(lo), residual(hi)
    # zbrent doesn't validate its own bracket (nor throw for a bad one) and
    # NaN comparisons never throw, so a mis-bracketed or NaN residual would
    # otherwise pass through silently instead of falling back to the guess.
    bracketed = isfinite(q_lo) && isfinite(q_hi) && (q_lo == 0.0 || q_hi == 0.0 || signbit(q_lo) != signbit(q_hi))
    solved = if bracketed
        try
            zbrent(residual, lo, hi, 1e-3)
        catch e
            e isa DomainError || rethrow()
            @warn "RootFindLeafTemperature: leaf_heat_balance raised $(typeof(e)) during root-find, falling back to leaf_temperature_guess" maxlog=40
            ustrip(u"K", leaf_temperature_guess)
        end
    else
        @warn "RootFindLeafTemperature: [$lo, $hi] K doesn't bracket a root (q_lo=$q_lo W, q_hi=$q_hi W), falling back to leaf_temperature_guess" maxlog=40
        ustrip(u"K", leaf_temperature_guess)
    end
    return solved * u"K"
end
