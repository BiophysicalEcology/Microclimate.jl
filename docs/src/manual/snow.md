# Snow

From a microclimate perspective, snow buffers soil temperature from the air above it, resulting in a 
highly stable 'subnivean' zone that remains near zero degrees celcius as phase transition occurs. Snow 
models range from phenomenological schemes that predict presence and depth, to detailed, physics-explicit 
models that capture 'ripening' and complex melting and refreezing processes. The snow algorithm in 
Microclimate.jl sits towards the latter end of the spectrum, with a fully physical model of heat transfer
between the atmosphere and the soil via the snow layer, including the melting and freezing processes,
with empricial functions of the ripening processes, specifically the evolution of snow density and snow
albedo as it ages.

## Names and symbols

| Symbol | Julia name | Meaning | Units |
| :----- | :--------- | :------ | :---- |
| ``\rho_{snow}`` | `snow_density` | snow density | g cm⁻³ |
| ``h`` | `current_depth` | snow height | cm |
| ``D`` | `snow_age` | age of the snow | days |
| ``\alpha`` | `albedo` | fractional surface albedo | 0–1 |
| ``d_{snow}`` | `days_since_snow` | days since the last snowfall | days |
| ``k_{snow}`` | `snow_conductivity` | snow thermal conductivity | W m⁻¹ K⁻¹ |

## Node scheme

[`SnowModel`](@ref) extends the soil depth nodes (`MicroModel.depths`) with up to 8
additional snow nodes at fixed depths (`snow_node_thresholds`), thermal properties set
to that of ice while snow is present. The nodes default to 2, 5, 10, 20, 50, 100, 200
and 300 cm (`DEFAULT_SNOW_NODE_THRESHOLDS`), with the original 0 cm surface
node continuing to act as the snow (or soil) surface node, and the link between the
deepest snow node and the true soil surface allowed to grow arbitrarily as the pack
deepens. As the snowpack grows and shrinks, `activate_snow_nodes!` checks
whether to bring more nodes into use or drop them; snow is set to zero below
`min_snow_depth` (2 cm by default, for numerical stability), with small rainfall
events that would otherwise produce a thinner layer accumulated hour to hour until the
threshold is met, so mass balance is conserved.

## Growth, density and albedo

Snow growth is driven by the rainfall input for each day (or hour, if hourly rainfall
data are supplied) and by `snow_temperature_threshold`, the air temperature at which
rain falls as snow instead. All of a day's snow is assumed to fall at midnight when
daily rainfall is used, and the `undercatch` parameter scales precipitation input to 
account for wind-driven gauge under-catch (Rasmussen et al. 2012).

By default, snow density is held constant at `snow_density` (0.375 g/cm³ by default).
Setting `density_function = (a, b, c, d)` switches to one of two time-varying forms:
with `c == 0`, density grows linearly with snow age, `min(0.9167, a·age + b)`
g/cm³; with `c > 0`, density follows the nonlinear asymptotic form of Sturm et al.
(2010, eq. 5):

```math
\rho_{snow} = (\rho_{max} - \rho_0)\left[1 - \exp(-k_1 h - k_2 D)\right] + \rho_0
```

with `density_function = (ρmax, ρ0, k1, k2)`: ``\rho_{max}`` the maximum allowable
density, ``\rho_0`` the initial density, ``h`` the snow height (cm) and ``D`` the age
of the snow (days — note this differs from Sturm et al. (2010), who use day of year).
Typical parameters: ``\rho_{max} = 0.9167``, ``\rho_0 = 0.27``, ``k_1 = 0.004``,
``k_2 = 0.009`` g/cm³.

Fractional surface albedo changes with the age of the snowpack,

```math
\alpha = \frac{-9.874\ln(d_{snow}) + 78.3434}{100}
```

where ``d_{snow}`` is the number of days since the last snowfall (floored at 0.3 to
keep the logarithm finite), based on regressions fitted to Anderson (2006, fig. A4),
reverting to the site's own albedo when snow is absent. Snow thermal conductivity, by
default, follows a cubic polynomial in density (Aggarwal 2009):

```math
k_{snow} = 0.00395 + 0.00084\rho - 1.7756\times10^{-6}\rho^2 + 3.80635\times10^{-9}\rho^3
```

with ``\rho`` = `snow_density` in kg/m³; `snow_conductivity` overrides this with a
fixed value instead, if set to a nonzero conductivity.

## Melt

Snow melt follows the heat load computed by the substrate (including snow) heat
budget, plus melting attributable to rain. When the temperature at a node rises, the
temperature difference between timesteps is multiplied by the snow heat capacity
(including the heat of fusion — see
[Phase transition](phase_transition.md#Apparent-heat-capacity,-inside-the-ODE-(snow))),
snow density, and the distance
between nodes, to get the heat input; dividing by the heat of fusion gives the mass of
snow melted, with the node's temperature reset to 0°C. Rain-melt follows Anderson
(2006): the product of rainfall and mean daily air temperature is multiplied by
`rain_melt_factor` (default 0.0125) whenever rain falls at air temperatures above
`snow_temperature_threshold` — the only purely empirical (as opposed to physically
derived) part of the scheme.

## Background

Kearney (2020) integrated this snow scheme (as originally implemented in NicheMapR's
Fortran) with the gridMET daily historical climatology of the continental USA and
tested it against hourly snow depth, snow water equivalent, and soil
temperature/moisture observations from 590 SCAN/SNOTEL sites (1979–2017), finding
predicted snow depth and soil temperature within about 10–15% of the observed range
(correlation around 0.85–0.96 at 5–50 cm depth) without site-specific tuning, and a
buffering effect on minimum soil temperature of roughly 10°C on average, as much as
19°C at some sites — shallow soil temperature warmed at roughly one-third the rate of
air temperature at snow-affected sites over the study period, a decoupling that
persisted even under a simulated 3°C warming of the historical record. This
Microclimate.jl port has since had further corrections made to it, so these figures
describe the scheme's original validation — ongoing comparisons against real station 
can be found in [MicroclimateTests.jl](https://github.com/BiophysicalEcology/MicroclimateTests.jl).

## References

See the [References](references.md) page for the full bibliography.
