# Soil moisture

Campbell (1985) developed an algorithm to compute soil water balance as a function of
infiltration and evapotranspiration, including a full soil-plant-atmosphere continuum
(SPAC) as a function of soil properties. This page summarises the theory and equations
of Campbell's Program 11.1 (C11.1) as implemented by [`CampbellSoilHydraulics`](@ref).

Like [soil thermal properties](soil_thermal_properties.md), the soil hydraulic properties 
can be made to vary with depth, specified per node on [`CampbellHydraulicProfile`](@ref) 
(`air_entry_water_potential`, `saturated_hydraulic_conductivity`, `campbell_b_parameter`, 
`root_density`). These can be taken from Campbell and Norman's (1998) Table 9.1 soil-texture 
values, or derived from soil texture data by a pedotransfer function (see Kearney and Maino 
2018, and below).

## Names and symbols

Campbell's (1985) own notation is used in the equations below; these are the Julia
names/fields that carry the same quantities (`i` subscripts a depth node):

| Symbol | Julia name | Meaning | Units |
| :----- | :--------- | :------ | :---- |
| ``\psi_i`` | `water_potential` | soil water potential at node `i` | J kg⁻¹ |
| ``\theta_i`` | `soil_moisture` | volumetric soil moisture at node `i` | m³ m⁻³ |
| ``\theta_{s,i}`` | (from `bulk_density`/`mineral_density`) | saturation water content, ``1-\rho_b/\rho_s`` | m³ m⁻³ |
| ``k_i`` | `hydraulic_conductivity` | hydraulic conductivity at node `i` | kg s m⁻³ |
| ``\psi_{e,i}`` | `air_entry_water_potential` | air entry potential | J kg⁻¹ |
| ``k_{s,i}`` | `saturated_hydraulic_conductivity` | saturated hydraulic conductivity | kg s m⁻³ |
| ``b_i`` | `campbell_b_parameter` | Campbell's dimensionless exponent | – |
| ``L_i`` | `root_density` | root length density at node `i` | m m⁻³ |
| ``R_w`` | `root_resistance` | resistance per unit root length | m³ s⁻¹ kg⁻¹ |
| ``r_1`` | `root_radius` | root radius | m |
| ``R_L`` | `leaf_resistance` | leaf resistance | m⁴ s⁻¹ kg⁻¹ |
| ``\psi_c`` | `stomatal_closure_potential` | critical leaf water potential | J kg⁻¹ |
| ``n`` (eq. 15–16 only) | `stomatal_stability_parameter` | stomatal-closure stability exponent | – |
| ``LAI`` | `leaf_area_index` | leaf area index | m² m⁻² |

Note `n` is reused for two different things in Campbell's own notation: the exponent
``n=2+3/b`` derived from `campbell_b_parameter` (eqs. 3, 22), and the unrelated
stomatal-closure stability exponent `stomatal_stability_parameter` (eqs. 15–16) — kept
as in the original rather than introduced here.

The overall driving equation for transpiration ``E`` (kg m⁻² s⁻¹) is

```math
E = \frac{\psi_{xL} - \psi_{L}}{R_L} = \frac{\psi_{xr} - \psi_{xL}}{R_x} = \frac{\psi_{xr} - \psi_{r}}{R_r} = \frac{\psi_{s} - \psi_{r}}{R_s}
```

where ``R`` is resistance (m⁴ s⁻¹ kg⁻¹), ``\psi`` is water potential (J kg⁻¹), and the
subscripts are xylem (``x``), leaf (``L``), root (``r``) and soil (``s``) — a chain of
resistances in series from the soil, through the root and xylem, to the leaf, driving
evaporation ``E``.

Soil moisture is solved on the same depth nodes as [soil temperature](soil_thermal_properties.md)
— `MicroModel.depths` (19 nodes by default, `DEFAULT_DEPTHS`, but any node
count or spacing the user configures) — plus one extra node beyond the deepest
supplied depth, whose hydraulic properties are copied from the last real layer and
which is held saturated throughout as the lower boundary condition. Each hour, after
soil temperature has been solved, `infiltration_step!` runs:

1. initialise water potential and hydraulic conductivity;
2. initialise the root water uptake variables;
3. partition potential evapotranspiration into potential evaporation and transpiration;
4. compute plant water uptake;
5. calculate actual transpiration rate and water extraction by roots from each layer;
6. solve for the mass balance of water in the soil (`soil_water_balance!`).

## Initialisation

Water potential ``\psi`` (J kg⁻¹) and hydraulic conductivity ``k`` (kg s m⁻³) at each
node ``i`` are computed from the previous hour's volumetric soil moisture ``\theta``
(m³ m⁻³):

```math
\psi_i = \psi_{e,i} (\theta_i/\theta_{s,i})^{-b_i}, \qquad
k_i = k_{s,i} (\psi_{e,i}/\psi_i)^n
```

where ``\psi_e`` is the air entry potential, ``\theta_s = 1 - \rho_b/\rho_s`` is the
saturation water content (``\rho_b`` bulk density, ``\rho_s`` mineral density),
``n = 2 + 3/b``, and ``b`` is Campbell's dimensionless exponent
(`campbell_b_parameter`).

Root resistance at each node is

```math
R_{r,i} = R_w / (L_i \Delta z)
```

where ``R_w`` (`root_resistance`) is the resistance per unit root length and ``L``
(`root_density`) the root-length density; layers without roots get an arbitrarily
large resistance. Water uptake from each layer is

```math
E_i = \frac{k_{r,i}\psi_{r,i} - k_{s,i}\psi_{s,i}}{B_{z,i}}, \qquad
B_{z,i} = \frac{(1-n)\ln(\pi r_1^2 L_i)}{4\pi L_i \Delta z}
```

(cylindrical-root geometry, `root_radius` = ``r_1``).

## Partitioning evapotranspiration

The soil heat budget's evaporative heat loss, converted to a mass flux of total
evapotranspiration ``E_T``, is partitioned into potential evaporation on the
assumption that ``E_P/E_T = e^{-0.82\,LAI}`` (`leaf_area_index`) — i.e. the ratio of
transpiration to evapotranspiration equals the ratio of radiation intercepted by
leaves to total incident radiation:

```math
E_P = e^{-0.82\,LAI}\,E_T, \qquad T_P = E_T - E_P
```

## Plant water uptake

Transpiration is the sum of the per-layer extraction rates, ``E = \sum E_i``. If axial
resistances are small compared with the others, this can be solved for the root xylem
potential ``\psi_{xr}``:

```math
\psi_{xr} = \frac{-E + \sum[\psi_{s,i}/(R_{s,i}+R_{r,i})]}{\sum[1/(R_{s,i}+R_{r,i})]}
```

with soil resistance at each layer ``R_{s,i} = B_{z,i}/k_i``, and a weighted mean soil
water potential

```math
\bar\psi_s = \frac{\sum[\psi_{s,i}/(R_{s,i}+R_{r,i})]}{\sum[1/(R_{s,i}+R_{r,i})]}
```

(the denominator is the weighted mean root-soil resistance ``\bar R_s``). Stem
resistance is neglected. Leaf water potential is then

```math
\psi_L = \bar\psi_s - E(\bar R_s - R_L)
```

— at ``E = 0``, ``\psi_L`` equals ``\bar\psi_s``: leaf water potential equilibrates to
the weighted mean soil water potential.

## Stomatal closure and actual transpiration

Transpiration is inversely related to stomatal resistance, which in turn varies with
leaf water potential:

```math
r_{vs} = r_{vs}^0 \left[1 + (\psi_L/\psi_c)^n\right]
```

where ``\psi_c`` is `stomatal_closure_potential` and ``n`` is
`stomatal_stability_parameter` (typically 3–20). If leaf boundary-layer resistance is
negligible, ``E = E_P/(\psi_L/\psi_c)^n``. A Newton-Raphson iteration finds the
``\psi_L`` that balances water supply and demand (see `infiltration_step!`); once
found, actual transpiration is ``T_R = T_P/(1+X_\psi)`` and water extracted per layer
is ``E_i = (\psi_{s,i} - \psi_L - R_L T_R)/(R_{r,i}+R_{s,i})``.

## Soil water mass balance

`soil_water_balance!` solves the simultaneous per-node balances with the
Newton-Raphson method and the tridiagonal (Thomas) algorithm, accounting for liquid
flow, gravitational flow, vapour flow and root water extraction. Vapour flux at the
surface is ``J_{v,1} = E_P (h_2-h_a)/(1-h_a)`` (``h_a`` the fractional atmospheric
relative humidity); below the surface, ``J_{v,i} = k_v(h_{i+1}-h_i)`` with vapour
conductivity ``k_v = D_v c_v' b \phi_g^m \Delta z``. The soil hydraulic capacity at
each node is ``C_i = -v_i\theta_i/(b_i\psi_i\Delta t)``. `MatricPotentialAlgorithm`
implements this using matric potential as the dependent variable directly (Campbell's
Program 8.1 — most efficient in dry soil); `MatricFluxPotentialAlgorithm` uses matric
flux potential instead (Program 8.2 — fewer iterations, more stable in wet soil).

The emergent behaviour: when soil water potential is uniform, water is extracted in
proportion to rooting density. As soil in high-root-density layers dries, water is
increasingly taken from layers with fewer roots, and the mean soil water potential
``\bar\psi_s`` declines even while some roots still sit in moist soil. Plant water
potential drops, driving stomatal closure and reduced transpiration.

## Prescribed or dynamic

Running the full Campbell solve every hour is optional — `MicroConfig.soil_moisture_strategy`
chooses between it and a constant/prescribed alternative:

- [`PrescribedSoilMoisture`](@ref) (the default) skips the solver entirely. Surface
  wetness for the evaporation calculation comes straight from
  `environment_daily.soil_wetness`; soil moisture itself is held at
  `initial_soil_moisture` unless a `precomputed_soil_moisture` matrix (depths × days)
  is supplied, in which case that day's column overrides it — a way to drive the
  model from externally computed or observed soil moisture without running Campbell's
  algorithm at all.
- [`DynamicSoilMoisture`](@ref) runs `infiltration_step!`/
  `soil_water_balance!` every hour as described above, tracking surface
  wetness as an evolving state rather than reading it from the environment.
  `moisture_tolerance`, `moisture_max_iterations` and `moisture_timestep` (≤ 1 hour)
  tune the Newton-Raphson sub-stepping.

## Background

Kearney and Maino (2018) tested this scheme (with soil hydraulic parameters derived
by pedotransfer from gridded soil texture/bulk-density products — SLGA and SoilGrids —
rather than by hand) against three years of root-zone soil moisture at 35 OzNet
monitoring sites, cosmic-ray neutron soil moisture (CosmOz), and satellite soil
moisture (ASCAT), all driven by continent-scale gridded weather data, finding
predictions within experimental error of the observations at every scale tested: root
mean square error around 0.07–0.08 m³ m⁻³ and a correlation coefficient around
0.65–0.76 against the in-situ and cosmic-ray datasets, and next-generation
depth-specific soil products (SLGA, SoilGrids) improving predictions by roughly 17%
over an older chloropleth soil map, with the largest gains at depth. That study
predates this Julia implementation and corrections made since, so these figures
describe the algorithm's original validation, not a from-scratch re-validation of this
port.

## References

See the [References](references.md) page for the full bibliography.
