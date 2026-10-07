# Microclimate

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://BiophysicalEcology.github.io/Microclimate.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://BiophysicalEcology.github.io/Microclimate.jl/dev)
[![CI](https://github.com/BiophysicalEcology/Microclimate.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/BiophysicalEcology/Microclimate.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/BiophysicalEcology/Microclimate.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/BiophysicalEcology/Microclimate.jl/tree/main)

Microclimate modelling in the Julia language: hourly above- and below-ground
microclimates — radiation, wind, soil heat and moisture, snow and canopy — from site,
weather and soil inputs. A Julia implementation of the microclimate model of
[NicheMapR](https://github.com/mrke/NicheMapR).

```julia
using Pkg
Pkg.add(url = "https://github.com/BiophysicalEcology/Microclimate.jl")
```

```julia
using Microclimate

out = solve(example_microclimate_problem())
```

See the [docs](https://BiophysicalEcology.github.io/Microclimate.jl/stable) for more,
starting with [Get started](https://BiophysicalEcology.github.io/Microclimate.jl/stable/get_started/).
