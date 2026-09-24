# Soil hydraulics algorithms

[Configuring the model](../manual/configuring_the_model.md) emphasises that composing the
simulation from swappable models makes it cheap to ask "what changes if I use a
different algorithm?" This tutorial illustrates how this works for
[`CampbellSoilHydraulics`](@ref)'s two independent algorithm choices.

## Infiltration algorithm

`infiltration_algorithm` selects how [Campbell's soil water mass balance](../manual/soil_moisture.md#Soil-water-mass-balance)
is solved: [`MatricPotentialAlgorithm`](@ref) (Campbell's Program 8.1, using matric
potential as the dependent variable directly — most efficient in dry soil) or
[`MatricFluxPotentialAlgorithm`](@ref) (Program 8.2, using matric flux potential —
fewer iterations, more stable in wet soil). At moderate-to-wet soil moisture the two
should agree closely:

```@example algorithms
using Microclimate, Unitful

config = MicroConfig(; soil_moisture_strategy = DynamicSoilMoisture())
n = length(Microclimate.DEFAULT_DEPTHS)
wet_initial_moisture = fill(0.30, n)

function run_with(infiltration_algorithm)
    soil_hydraulic_model = example_soil_hydraulic_model(; infiltration_algorithm)
    solve(example_microclimate_problem(; soil_hydraulic_model, config,
        initial_soil_moisture = wet_initial_moisture))
end

psi_out = run_with(MatricPotentialAlgorithm())
phi_out = run_with(MatricFluxPotentialAlgorithm())
maximum(abs.(psi_out.soil_moisture .- phi_out.soil_moisture))
```

In dry soil, `MatricFluxPotentialAlgorithm` is documented (Campbell 1985) to be less
accurate — a case where the choice actually matters rather than being interchangeable.

## Rainfall entry mode

`rainfall_entry_mode` selects how rain reaches the soil column each hour:
[`PoolCapacityRainfall`](@ref) (the default — a surface pool fills the top node's own
storage capacity before infiltration runs), [`ImplicitFluxRainfall`](@ref) (rain
enters as a flux boundary condition inside the same implicit solve infiltration
itself uses, alongside evaporation), or [`RateLimitedFrontRainfall`](@ref) (a
stateless, rate-limited multi-layer wetting-front march). Comparing them under a
rainy period:

```@example algorithms
environment_daily = example_daily_environment(; rainfall = fill(30.0, 12)u"kg/m^2")

function run_with_rain(rainfall_entry_mode)
    soil_hydraulic_model = example_soil_hydraulic_model(; rainfall_entry_mode)
    inputs_kw = (; environment_daily, initial_soil_moisture = fill(0.15, n))
    solve(example_microclimate_problem(; soil_hydraulic_model, config, inputs_kw...))
end

pool_out = run_with_rain(PoolCapacityRainfall())
flux_out = run_with_rain(ImplicitFluxRainfall())
(pool_out.soil_moisture[end, 1], flux_out.soil_moisture[end, 1])
```

## Frozen soil

Once a layer is frozen, water is treated as unable to move through it — both through
the soil matrix and through roots — so infiltration and root uptake stop at whichever
depth freezing reaches, rather than continuing through ice-filled pores. See
[Phase transition](../manual/phase_transition.md) for the full mechanism, and how it
interacts with [Snow](../manual/snow.md), since a snow-covered, freezing soil column
is exactly where this matters.

## Next steps

- [Soil moisture](../manual/soil_moisture.md) — the full Campbell (1985) theory these
  algorithms implement.
- [Configuring the model](../manual/configuring_the_model.md) — the general pattern
  behind comparing model choices this way.
