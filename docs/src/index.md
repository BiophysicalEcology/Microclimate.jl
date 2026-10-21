```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "Microclimate.jl"
  text: "The climate near the ground"
  tagline: "hourly above- and below-ground microclimates — radiation, wind, soil heat and moisture, snow and canopy — for biophysical ecology, with units."
  actions:
    - theme: brand
      text: Get Started
      link: /get_started
    - theme: alt
      text: View on Github
      link: https://github.com/BiophysicalEcology/Microclimate.jl
    - theme: alt
      text: API Reference
      link: /api

features:
  - title: 🌡️ Soil heat balance
    details: Hourly soil temperature at any depth from a 1-D node-based heat budget, with <a class="highlight-link">Campbell and de Vries</a> thermal properties that vary with moisture and texture.
    link: /manual/soil_thermal_properties
  - title: 💧 Soil moisture
    details: Infiltration, redistribution and root uptake from <a class="highlight-link">Campbell's (1985) Soil-Plant-Atmosphere-Continuum model</a>, tested at continental scale against SCAN, CosmOz and satellite soil moisture.
    link: /manual/soil_moisture
  - title: ❄️ Snow
    details: A snow node scheme that grows and shrinks with accumulation and melt, buffering soil temperature from the air above it.
    link: /manual/snow
  - title: 🌳 Canopy
    details: A layer-resolved multilayer canopy — <a class="highlight-link">two-stream shortwave, longwave exchange, wind attenuation and leaf temperature</a> — for vegetated sites.
    link: /manual/canopy
  - title: ☀️ Radiation
    details: Clear-sky solar geometry from <a class="highlight-link">SolarRadiation.jl</a>, cloud-adjusted shortwave and a full longwave/sky-emissivity budget.
    link: /manual/radiation
  - title: 🧩 Composable models
    details: Every physical process — soil, snow, canopy, evaporation, convergence — is a swappable component, so alternative model choices can be compared by changing one argument.
    link: /manual/configuring_the_model
  - title: 📏 Units
    details: Every input and output uses <a class="highlight-link">Unitful.jl</a>, and calculations build on <a class="highlight-link">FluidProperties.jl</a> for air and water properties.
    link: /get_started
---
```

## How to install Microclimate.jl?

Microclimate.jl can be installed from the Julia REPL:

```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/BiophysicalEcology/Microclimate.jl")
```

## Model

Microclimate.jl is a Julia implementation of the microclimate model of
[NicheMapR](https://github.com/mrke/NicheMapR) (Kearney and Porter 2017), extended
with a new multilayer canopy. See the [Introduction](manual/introduction.md) for what
a microclimate is, how this package relates to NicheMapR and other land-surface
models, and how the package is organised.
