# Introduction

Microclimates are the thermal, hydric, and radiative conditions in the first metre or
so above and below the earth's surface — "the climate near the ground". The topic
encompasses the effects of terrain and vegetation on radiation, air temperature, wind
speed and humidity, as well as the dynamics of soil temperature, soil moisture and snow, 
and phase transitions of water including freeze/thaw of water in the soil and surface 
condensation (dew and frost). 

Microclimates are the physical conditions experienced by organisms, and
constraining their energy and mass budgets and ultimately their
behaviour, distribution and abundance. Weather station measurements are made 1–2 m
above the ground specifically to avoid local terrain and vegetation effects, so that
they give a regionally representative measure — but it is the local, near-ground
conditions that matter to an organism.

Microclimate.jl computes those near-ground conditions — short- and long-wavelength
radiation, air temperature, wind speed, humidity, substrate temperature, soil
moisture, snow, dew and frost — from weather, terrain, soil and vegetation inputs. 
It is designed to compute the vertical distribution of microclimatic conditions near 
and below the ground at a point, given the properties of the habitat and information 
on the weather conditions 1–2 m above the ground. It does not compute meso-scale phenomena 
that depend on surrounding conditions, such as cold-air drainage, though such effects can be 
added directly via the driving weather input data. It  assumes habitat properties are uniform 
across an infinite plane, so it does not capture spatial dynamics from lateral heat or moisture 
flow. Some pixel-aware effects, like cold-air drainage and runoff routing between neighbouring 
cells, are captured by [MicroclimateMapper.jl](https://github.com/BiophysicalEcology/MicroclimateMapper.jl), which integrates Microclimate.jl with existing weather, climate and terrain data sets via the
[RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl), [Rasters.jl](https://github.com/rafaqz/Rasters.jl) and [Geomorphometry.jl](https://github.com/Deltares/Geomorphometry.jl) packages.

## Land-surface models

A number of other models compute microclimatic or land-surface conditions from
atmospheric and terrain data, differing in scope and audience. Full land-surface
models — CABLE, JULES, Noah-MP, and, in the Julia ecosystem,
[ClimaLand.jl](https://clima.github.io/ClimaLand.jl/stable/) — are built to provide
biosphere feedbacks to general circulation models, and are therefore not directly
tailored to producing microclimatic outputs required for solving the energy and mass 
budgets of organisms. They are generally more difficult to set up as stand-alone tools 
without specialist skills. TrenchR's
["Estimating microclimates"](https://cran.r-project.org/web/packages/TrenchR/vignettes/MicroclimateTutorial.html)
vignette offers a general-purpose set of R functions for scaling weather
station data to microclimate variables, built on Campbell and Norman (1998). Microclimate.jl
and its ancestor [NicheMapR's](https://github.com/mrke/NicheMapR)
microclimate model, occupy the middle ground of general-purpose and ecologically focused models 
(hourly time step, fine-scale topographic adjustments, detailed soil profiles), accessible
without land-surface-model expertise.

## Relationship to NicheMapR, micropoint/microclimc and SHAW

The soil and above-ground core of Microclimate.jl is a Julia port of the NicheMapR
microclimate model (Kearney and Porter 2017), and the initial port was checked for fidelity 
against the original Fortran by comparing against NicheMapR runs (see `test/R/` in the package repository).

The more modular structure of Microclimate.jl makes it easier to build new functionality,
including interactions with other parts of the BiophysicalEcology software ecosystem. 
An example is the inclusion of a multilayer canopy ([Canopy](canopy.md)) scheme inspired
by Maclean's [micropoint](https://github.com/ilyamaclean/microclimc)/
[microclimlearn](https://rpubs.com/ilyamaclean/microclimlearn), and the SHAW model
(Flerchinger), combined with leaf temperature calculations that use
[HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl). See [Canopy](canopy.md).

## Package ecosystem

Microclimate.jl is one package in the growing BiophysicalEcology.jl ecosystem,
where features of the NicheMapR package have been separated, modularised and
integrated into a more composable and extensible system:

```
FluidProperties.jl   SolarRadiation.jl   BiophysicalGeometry.jl
        │                    │                     │
        ├────────────────────┴──────────┐          │
        ▼                                ▼          ▼
  Microclimate.jl  ◀─────────────────▶  HeatExchange.jl
        │                                │
        ▼                                ▼
  MicroclimateMapper.jl  ──────▶  AnimalMapper.jl / PlantMapper.jl / MicrobeMapper.jl
                                                │
                                                ▼
                                          NicheMapper.jl
```

(arrows show the direction of dependency or data flow; Microclimate.jl and
HeatExchange.jl depend on each other in different ways, see below.)
Microclimate.jl itself depends on
[FluidProperties.jl](https://biophysicalecology.github.io/FluidProperties.jl/stable/)
(air and water properties) and
[SolarRadiation.jl](https://biophysicalecology.github.io/SolarRadiation.jl/stable/)
(clear-sky solar geometry), and, for the canopy's leaf-temperature solve, on
[BiophysicalGeometry.jl](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl)
(leaf shape) and [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl)
(leaf convection/evaporation) — see [Canopy](canopy.md).

Two companion packages sit either side of Microclimate.jl in this ecosystem.
[MicroclimateMapper.jl](https://github.com/BiophysicalEcology/MicroclimateMapper.jl) 
(in development) connects this model to gridded spatial and
temporal driving datasets, supplying its inputs — the Julia analogue of NicheMapR's
`micro_global`/`micro_usa`/`build.global.climate` machinery. Downstream,
[HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) computes
organism heat, water and metabolic budgets and requires microclimate conditions —
air temperature, wind, humidity, radiation — as its environmental driving data; a
`MicroResult` is the natural source for that input, the same relationship
NicheMapR's ectotherm and endotherm models have to its own microclimate output. Beyond
that, MicroclimateMapper.jl and HeatExchange.jl together feed a planned layer of
organism-specific packages (AnimalMapper.jl, PlantMapper.jl, MicrobeMapper.jl), which
in turn feed a NicheMapper.jl umbrella package.

## How the package is organised

Physical processes are individually swappable models composed together into a
[`MicroModel`](@ref) — see [Configuring the model](configuring_the_model.md) for the
design and why it's structured this way. The manual covers each subsystem in turn:

- [Configuring the model](configuring_the_model.md) — assembling a `MicroModel`
- [Diel curves and time handling](diel_time_handling.md) — turning daily/monthly data
  into hourly cycles
- [Boundary layer](boundary_layer.md) — wind, temperature and humidity profiles
- [Radiation](radiation.md) — shortwave and longwave budgets at the surface
- [Soil thermal properties](soil_thermal_properties.md) — the soil heat balance
- [Phase transition](phase_transition.md) — freezing/thawing and its effect on water
  transport
- [Soil moisture](soil_moisture.md) — infiltration and the soil water balance
- [Snow](snow.md) — the snow node scheme
- [Canopy](canopy.md) — the multilayer canopy
- [Evaporation and condensation](evaporation_condensation.md) — surface latent heat
  flux and dew/frost formation

## References

See the [References](references.md) page for the full bibliography.
