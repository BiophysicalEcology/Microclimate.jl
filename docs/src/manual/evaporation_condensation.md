# Evaporation and condensation

Surface latent heat exchange — evaporation from wet ground and dew/frost formation —
draws on the hydric and thermal properties of air computed by
[FluidProperties.jl](https://biophysicalecology.github.io/FluidProperties.jl/stable/)
(dry/wet air density, viscosity, vapour density, saturation vapour pressure); see that
package's docs for the underlying equations.

## Names and symbols

| Symbol | Julia name | Meaning | Units |
| :----- | :--------- | :------ | :---- |
| ``m_{evap}`` | `evaporation_mass_flux` | evaporated mass flux | kg s⁻¹ m⁻² |
| ``Q_{evap}`` | `Q_evaporation` | evaporative heat loss | W m⁻² |
| ``\Delta H_{vap}`` | `enthalpy_of_vaporisation(T)` (FluidProperties.jl) | latent heat of vaporisation | J kg⁻¹ |
| ``E`` | `dew_frost_energy_flux` | dew/frost energy flux (Monteith convention: +ve loss, −ve condensation) | W m⁻² |
| ``q_g^\star`` | `q_g★` | surface saturated specific humidity | kg kg⁻¹ |
| ``q_1`` | `q_1` | ambient (actual) specific humidity | kg kg⁻¹ |
| ``q_1^\star`` | `q_1★` | ambient saturated specific humidity | kg kg⁻¹ |
| ``s`` | `(q_g★ - q_1★) / ΔT` | slope of saturated specific humidity vs. temperature, finite-differenced between two `wet_air_properties` calls | kg kg⁻¹ K⁻¹ |
| ``\gamma`` | `specific_heat/λ` | psychrometric constant | kg J⁻¹ kg⁻¹ K⁻¹ |
| ``R`` | `R_N` (net radiation) | net radiation | W m⁻² |
| ``G`` | `G` (ground heat flux term) | heat flux to/from the surface | W m⁻² |
| ``\delta q`` | `δq` | air's own saturation specific-humidity deficit | kg kg⁻¹ |
| ``C_e`` | `C_E` | bulk transfer coefficient for heat/vapour | – |
| ``\bar u`` | `wind_speed` | wind speed at reference height | m s⁻¹ |

## Evaporation

[`BulkTransferEvaporation`](@ref) computes surface evaporation as a Penman-style
mass-transfer flux from the surface and air vapour densities, scaled by surface
wetness:

```math
m_{evap} = w\, h_d (v_{d,surf} - v_{d,air})
```

where ``w`` = `soil_wetness` (0–1) is the fraction of the unit surface area acting as a free water
surface. The mass transfer coefficient ``h_d`` follows from the convective heat
transfer coefficient via the Chilton-Colburn analogy (Bird et al. 2002), which assumes
the temperature and humidity profiles share the same shape:

```math
h_d = \frac{h_c}{c_p \rho_{air}}\left(\frac{0.71}{0.60}\right)^{0.666}, \qquad
h_c = \max\!\left(\left|\frac{Q_{conv}}{T_{loc}-T_{ref}}\right|, 0.5\right)
```

(see `calc_mass_transfer_coefficient`,
`calc_heat_transfer_coefficient`, `LEWIS_HEAT_TO_MASS_RATIO`).
Evaporated mass converts to a heat loss ``Q_{evap} = m_{evap}\,\Delta H_{vap}``, with
the latent heat of vaporisation itself temperature-dependent (see
`enthalpy_of_vaporisation` in FluidProperties.jl).

## Dew and frost

Dew forms on a surface (e.g., of the ground, a leaf or an insect) when the surface
temperature drops sufficiently below the air temperature such that the surrounding air
cools below the dew point. The dew point depends on the specific humidity of the air —
the mass of water per mass of moist air. Surface temperatures vary dramatically in
both space and time in natural microclimates, so dew formation is very sensitive to
the microclimatic setting.

Surface temperature is the outcome of a balance between incoming and outgoing
radiation, sensible heat exchange by convection, and the heat gained or lost to the
latent heat of vaporisation as water condenses or evaporates. The rate of evaporation
``E`` (kg s⁻¹, or mm assuming a water density of 1 g cm⁻³) can be expressed as the sum
of two terms (Garratt and Segal 1988; Monteith 1963):

```math
E = \frac{s}{s+\gamma}\,\frac{R-G}{\lambda} + \frac{\gamma}{s+\gamma}\,\rho\,\delta q\,C_e\,\bar u
```

In the first term, ``R`` is net radiation (W m⁻²), ``G`` the heat flux to the soil
surface (W m⁻²), ``\lambda`` the latent heat of vaporisation (J kg⁻¹), ``\gamma`` the
psychrometric constant ``C_p/\lambda`` (kg J⁻¹ kg⁻¹ K⁻¹), and ``s`` the change in
saturated specific humidity with temperature, from the reference height to the surface
(kg kg⁻¹ K⁻¹). In the second term, ``\rho`` is air density (kg m⁻³), ``\delta q`` the
difference between the saturated and actual specific humidity of the air (kg kg⁻¹),
``C_e`` a dimensionless bulk transfer coefficient for heat or water vapour that
depends on the roughness height, and ``\bar u`` the wind speed at the reference
height (m s⁻¹). [`GarrattSegalCondensation`](@ref) (the `MicroModel` default)
implements this directly (eq. 6b of Garratt and Segal 1988), as a Penman-style
combination of an energy-availability term and an aerodynamic term — the energy term
suppresses condensation whenever the surface is net gaining energy (the usual daytime
case). [`BulkTransferCondensation`](@ref) is a simpler alternative that reuses the
model's own [`BulkTransferEvaporation`](@ref) flux directly (condensing whenever the
saturated surface vapour density falls below the actual air vapour density): it's
self-consistent with the model's own boundary-layer conductance, but has no
energy-availability gating, so it can't distinguish nocturnal radiative dew from any
other surface-cooler-than-air event and can produce implausibly large fluxes when
boundary-layer conductance is high. [`NoCondensation`](@ref) switches ground dew/frost
off. The equivalent leaf-level process, for a [`MultilayerCanopy`](@ref), is
[`MonteithLeafCondensation`](@ref).

Dew forms whenever ``E`` is negative — i.e. when condensation is occurring rather than
evaporation. The source of dew may be moisture in the soil (distillation) or
condensation from the atmosphere (dewfall). When the air is saturated (``\delta q =
0 kg kg⁻¹``) the second term drops out; for an isolated object such as a leaf, ``G`` drops
out. The ``s`` term increases with temperature, because warm air holds more water
vapour, but ``R`` also increases with temperature when the atmosphere is saturated,
because sky emissivity rises with absolute humidity — these balancing factors mean the
maximum potential rate of dew formation is largely independent of temperature, around
0.07–0.08 mm/h (Garratt and Segal 1988; Monteith 1963). The actual rate is usually
lower because of atmospheric relative humidity, wind and heat from the ground: wind
speed can augment or diminish dew formation, since it both replenishes moist air and
enhances convective heat exchange. The processes involved in dew formation interact in
complex ways that cannot be easily intuited, which is exactly why a mechanistic model
of the full surface energy balance — rather than an empirical dew-point rule — is
useful here.

## References

See the [References](references.md) page for the full bibliography.
