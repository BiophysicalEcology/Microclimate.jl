# Boundary layer

Air temperature and wind speed vary with height above the ground as a function of the
conditions at the reference height, the surface temperature, and the roughness of the
surface. Microclimate.jl computes vertical profiles of wind speed, air temperature and
relative humidity between the surface and the reference height with
[`atmospheric_surface_profile`](@ref)/[`atmospheric_surface_profile!`](@ref), via
[`MoninObukhov`](@ref) similarity theory (MOST).

## Names and symbols

| Symbol | Julia name | Meaning | Units |
| :----- | :--------- | :------ | :---- |
| ``u^*`` | `friction_velocity` | friction velocity | m s⁻¹ |
| ``L`` | `obukhov_length` | Monin-Obukhov length | m |
| ``\kappa`` | `karman_constant` | von Kármán constant | – |
| ``z`` | (from `heights`) | height above the surface | m |
| ``z_0`` | `roughness_height`/`roughness_length` | momentum roughness length | m |
| ``d_0`` | `displacement_height` | zero-plane displacement | m |
| ``z_H`` | (via `thermal_roughness_model`) | roughness length for heat transfer | m |

## Monin-Obukhov similarity theory

[`MoninObukhov`](@ref) solves for the friction velocity ``u^*`` and Obukhov length
``L`` iteratively (`calc_Obukhov_length`), so stable and unstable atmospheric
regimes are both handled directly by the same solve. Wind and temperature at a height
``z`` above a roughness length ``z_0`` follow the log-law with a stability correction:

```math
u(z) = \frac{u^*}{\kappa}\left[\ln\!\left(\frac{z}{z_0}\right) - \psi_m\right]
```

with ``\kappa`` the von Kármán constant and ``\psi_m`` the Businger-Dyer/Paulson
momentum stability correction (`calc_ψ_m`), following Businger et al. (1971)
and Dyer (1974) for the unstable branch and a Dyer (1974) linear/saturating form for
the stable branch — the two meet continuously at neutral (``L \to \pm\infty``). The
analogous heat correction ``\psi_h`` (`calc_ψ_h`) and the ``\Phi_h`` bulk
diffusivity multiplier (`calc_Φ_h`) place the surface temperature and
convective heat flux on the same log-law basis.

Convective heat transfer at the surface combines a bulk (log-law) and a sublayer
Stanton number (`bulk_stanton`, `sublayer_stanton`,
`convective_flux`) — the default [`SublayerStantonRoughness`](@ref)
roughness-to-heat-transfer correction, appropriate for bare surfaces. A
fixed-ratio alternative, [`ScalarRoughnessRatio`](@ref), is used for the canopy
surface instead (see [Canopy](canopy.md)), where a Stanton-number formula tuned to
bare ground is less appropriate.

```@example boundary_layer
using Microclimate, Unitful

boundary_layer_model = MoninObukhov()
site = example_site()
environment_instant = (;
    reference_temperature = 25u"°C",
    reference_wind_speed  = 2.0u"m/s",
    reference_humidity    = 0.6,
    atmospheric_pressure  = 101.325u"kPa",
    zenith_angle          = 45u"°",
)
profile = atmospheric_surface_profile(boundary_layer_model;
    heights = [0.01, 0.5, 2.0]u"m",
    reference_height = 2.0u"m",
    site,
    environment_instant,
    surface_temperature = 35u"°C",
)
profile.air_temperature
```

## References

See the [References](references.md) page for the full bibliography.
