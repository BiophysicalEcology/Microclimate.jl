# Radiation

Solar radiation has a dominating effect on the microclimates available at a site. The
clear-sky solar geometry, spectral irradiance and terrain effects (slope, aspect,
hillshade) are computed by [SolarRadiation.jl](https://biophysicalecology.github.io/SolarRadiation.jl/stable/)
— see that package's [Introduction](https://biophysicalecology.github.io/SolarRadiation.jl/stable/manual/introduction/)
for the solar geometry and atmospheric attenuation equations, following McCullough
and Porter (1971). Microclimate.jl adds two things SolarRadiation.jl deliberately
leaves out: adjustment of the clear-sky irradiance for cloud cover, and the longwave
radiation budget.

## Names and symbols

| Symbol | Julia name | Meaning | Units |
| :----- | :--------- | :------ | :---- |
| ``G_{clear}``, ``G_{cld}`` | `global_horizontal`, `global_radiation` | clear-sky / cloud-adjusted global irradiance | W m⁻² |
| ``c`` | `cloud_cover` | fractional cloud cover | 0–1 |
| ``k_t`` | (clearness index) | ratio of global to extraterrestrial irradiance | – |
| ``D`` | `diffuse_horizontal_irradiance` | diffuse horizontal irradiance | W m⁻² |
| ``k_d = D/G`` | `diffuse_fraction` | diffuse fraction of global irradiance | 0–1 |
| ``\sigma`` | `σ` | Stefan-Boltzmann constant | W m⁻² K⁻⁴ |
| ``\epsilon`` | `surface_emissivity`/`cloud_emissivity` | longwave emissivity | 0–1 |
| ``T_a`` | `reference_temperature` | reference-height air temperature | K |
| ``e_a`` | (from `wet_air_properties`) | vapour pressure of the air | kPa |
| ``\epsilon_{sky}`` | (from `atmospheric_radiation`) | clear-sky emissivity | 0–1 |
| ``s`` | `shade` | fractional vegetation shade | 0–1 |
| ``f`` | `sky_view_fraction` | fraction of sky unobscured by terrain | 0–1 |
| ``T_{sky}`` | `sky_temperature` | effective sky radiant temperature | K |

## Shortwave

[`AngstromMaxwellShortwave`](@ref) adjusts SolarRadiation.jl's clear-sky global
irradiance for cloud cover via the Ångström-Prescott scaling of [`Angstrom`](@ref)
(defaults matching Linacre 1992, eq. 5.33):

```math
G_{cld} = G_{clear} \left(a + b(1 - c)^{\gamma}\right)
```

where ``G_{clear}`` is the clear-sky global irradiance, ``c`` is `cloud_cover` as a
fraction, and ``a = 0.36``, ``b = 0.64``, ``\gamma = 1`` are
`Angstrom`'s own fields. The cloud-adjusted global irradiance is then split back into
direct and diffuse components via a Maxwell (1987) clearness index ``k_t`` passed to
[`ErbsDiffuseFraction`](@ref), a piecewise fit (Erbs et al. 1982) of diffuse fraction
against clearness index:

```math
k_d = \begin{cases}
1 - 0.09 k_t & k_t \le 0.22 \\
0.9511 - 0.1604 k_t + 4.388 k_t^2 - 16.638 k_t^3 + 12.336 k_t^4 & 0.22 < k_t \le 0.80 \\
0.165 & k_t > 0.80
\end{cases}
```

(see `shortwave_radiation!`).

## Longwave

All objects in a habitat emit longwave radiation at a rate proportional to the 4th
power of their temperature, ``\sigma \epsilon T^4``, with ``\sigma`` the
Stefan-Boltzmann constant, ``\epsilon`` the emissivity and ``T`` the temperature in
Kelvin. [`ViewFactorLongwave`](@ref) combines the
clear-sky downwelling longwave from an [`AbstractAtmosphericRadiationModel`](@ref)
(default [`CampbellNormanAtmosphericRadiation`](@ref)) with cloud emissivity, sky view
fraction, vegetation shade and hillshade contributions to give the net longwave at the
surface, computed by `precompute_longwave_sky` and `longwave_radiation`.

For clear skies, sky emissivity follows Campbell and Norman (1998, eq. 10.10):

```math
\epsilon_{sky} = 1.72 \left(\frac{e_a}{T_a}\right)^{1/7}
```

where ``T_a`` is `reference_temperature` in Kelvin and ``e_a`` its vapour pressure in
kPa ([`CampbellNormanAtmosphericRadiation`](@ref) implements this as
`atmospheric_radiation`, applying the ``\sigma T_a^4`` term for the final output). Cloud
radiation is approximated from ``T_a - 2\,\mathrm{K}``:

```math
\text{cloud\_radiation} = \sigma\, \epsilon_{cloud} (T_a - 2\,\mathrm{K})^4
```

with `cloud_emissivity` ``\epsilon_{cloud}`` the daily cloud emissivity forcing. `atmospheric_longwave`
(clear-sky) and `cloud_radiation` are combined by `clear_sky_fraction = 1 -
cloud_cover` and `cloud_cover`, then reduced by `shade` (fractional vegetation shade,
0–1) to give the sky's own contribution:

```math
\text{longwave\_radiation\_sky} = \bigl[\text{atmospheric\_longwave}\,(1-c) + \text{cloud\_radiation}\, c\bigr] (1 - s)
```

with ``c`` = `cloud_cover`, ``s`` = `shade`. Vegetation casting that shade is assumed
to radiate at unit emissivity, ``\sigma T_a^4`` (`hillshade_radiation`), contributing
`longwave_radiation_vegetation = shade * hillshade_radiation`; hillshade contributes
`longwave_radiation_hillshade = hillshade_radiation` directly.

`sky_view_fraction` (from `Site`) is the fraction of the sky not obscured by
surrounding terrain, typically derived from horizon angles as ``1 - \sum_i
\sin(\text{horizon\_angle}_i)/n`` for ``n`` horizon directions. Sky and vegetation
terms are weighted by `sky_view_fraction`, hillshade by the remainder:

```math
\text{incoming\_longwave} = \bigl(\text{longwave\_radiation\_sky} + \text{longwave\_radiation\_vegetation}\bigr)\,f + \text{longwave\_radiation\_hillshade}\,(1-f)
```

with ``f`` = `sky_view_fraction`. `longwave_radiation` then adds the ground's
own outgoing emission, `surface_radiation = σ * surface_emissivity * surface_temperature^4`,
weighted by `(1 - shade)` plus a shaded-ground term at `hillshade_radiation`, to give
the net longwave gain for the surface heat budget,
`net_longwave_radiation = incoming_longwave - longwave_radiation_ground`. The
effective *sky temperature* — an output of the model, `sky_temperature` in
`MicroResult` — is

```math
T_{sky} = \left(\frac{\text{incoming\_longwave}}{\sigma}\right)^{1/4}.
```

If `HourlyTimeseries.longwave_radiation` is supplied directly (measured sky longwave),
`precompute_longwave_sky` uses it in place of the cloud-cover-derived estimate for
`longwave_radiation_sky`.

## References

See the [References](references.md) page for the full bibliography.
