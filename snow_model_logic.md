# NicheMapR Snow Model Logic

## Overview

The snow model adds up to 8 snow nodes on top of the 10 soil nodes, giving 18 total nodes when snow is active. Snow accumulation, melt, phase change, and heat exchange with the soil are all computed hourly within the ODE framework.

---

## 1. Node Structure

### Configuration
- **Without snow**: 10 soil nodes, `T(1:10)`
- **With snow**: 18 nodes total — `T(1:8)` snow nodes + `T(9:18)` soil nodes

### Snow Node Depths (fixed thresholds)
| Node | Depth (cm) |
|------|------------|
| 1    | 2          |
| 2    | 5          |
| 3    | 10         |
| 4    | 20         |
| 5    | 50         |
| 6    | 100        |
| 7    | 200        |
| 8    | 300        |

Nodes are finer near the surface where temperature gradients are largest. Only nodes whose threshold is exceeded by the current snow depth are active (`maxsnode` tracks the count). Inactive snow nodes are set equal to the surface soil node temperature.

### Unified Depth Array
A single depth array `DEPP(1:18)` links snow and soil seamlessly:
- `DEPP(1:8)` = snow node depths (`snode` array)
- `DEPP(9:18)` = soil depths + current snow depth offset

This offset means the soil grid automatically adjusts downward as snow accumulates.

---

## 2. Snow Accumulation

### Precipitation Partitioning
If air temperature ≤ `snowtemp` (user parameter), precipitation falls as snow:

```
snowfall = (rainfall × rainmult × undercatch × 0.1) / snowdens
```

- `undercatch`: correction for gauge undercatch (< 1)
- Factor 0.1 converts mm rain to cm snow depth
- Vegetation intercepts a fraction (`intercept`) before it reaches the ground

### Snow Density Evolution
Density evolves either exponentially (if `densfun(3) > 0`):
```
snowdens = (densfun(1) - densfun(2)) × (1 - exp(-densfun(3)×cursnow - densfun(4)×snowage)) + densfun(2)
```
or linearly:
```
snowdens = min(0.9167, densfun(1)×snowage + densfun(2))
```
`snowage` increments by 1/25 each hour. When density changes between hours, the previous hour's snow depth is rescaled to conserve mass.

### Node Reassignment (SNOWLAYER)
Each hour, `maxsnode` is recalculated from current snow depth, and the `snode` depth array is reassigned bottom-up. Nodes that fall above the current snow surface are zeroed and their temperatures reset to the surface soil temperature.

See Section 7 for the full detail of SNOWLAYER logic.

---

## 3. Snow Melt

### Energy-Balance Melt
For each active snow layer with mean temperature > `melthresh` (0.4°C):

```
meltheat = (meanT - melthresh) × cpsnow × layer_depth × 10000 × snowdens
         + HTOFN × layer_depth × 10000 × snowdens
melt += meltheat / HTOFN
```

where:
- `cpsnow = 2.100×snowdens + (1.005 + 1.82×RW/(1+RW))×(1−snowdens)` — weighted specific heat of ice and humid air
- `HTOFN = 333.55 J/g` — latent heat of fusion
- Depth factor 10000 converts cm² to m²

Melt depth is then scaled by user parameter `snowmelt` (0–1). Unmelted energy enters `QFREZE`, a latent heat term that feeds back into the surface energy balance.

### Rainfall-on-Snow Melt
When air temperature ≥ `snowtemp` and rainfall is occurring:
```
rainmelt = rainmeltf × rainfall × T_air × 24 / 10 / snowdens
```
This adds to the melt total beyond the energy-balance melt.

### Net Snowpack Update
```
netsnow = snowfall - sublimation×(0.0001/snowdens) + rainmelt
cursnow = cursnow + netsnow - melted
```

---

## 4. Phase Change

### Concept
When a snow or soil layer crosses 0°C, latent heat is consumed or released. Temperatures are held at 0°C until the cumulative phase-change energy budget is exhausted.

### Detection
- **Freezing**: `meanT_past > 0` and `meanT ≤ 0`
- **Thawing**: `meanT_past ≤ 0` and `meanT > 0`

### Energy Accounting (snow layers)
```
qphase(layer) += (meanT_past - meanT) × layermass × cpsnow
```
Temperatures of freezing/thawing nodes are pinned to 0°C:
```
T(node) = 0.0,  TT(node) = 0.0
```

### Capacity Limit
`sumphase` (cumulative phase energy) is compared against `HTOFN × total_water_mass`. If the limit is exceeded, node temperatures are reset to −0.5°C and the accumulator is cleared — preventing energy conservation violations in deep snow.

### Soil Phase Change
Identical logic applies to soil nodes (tracked separately in `qphase2`, `sumphase2`), with the maximum latent heat per layer capped at `HTOFN × layermass`.

### Snowfall Temperature Enforcement
When new snow falls, all corresponding nodes are forced to ≤ 0°C to prevent freshly-arrived snow from inheriting warm pre-existing temperatures.

---

## 5. Snow–Soil Heat Exchange

### Surface Energy Balance
The surface ODE is:
```
dT/dt[surface] = (Q_solar + Q_IR + Q_conduction + Q_convection + Q_freeze - Q_evap) / WC[j]
```
`j` is the first active node below the surface. When snow is present, `j` is the topmost active snow node; when no snow is present, `j` is soil node 1.

`QFREZE` (latent heat from snow melt/freeze) enters directly into this balance, coupling the snowpack thermodynamics to the surface temperature evolution.

### Conduction at the Snow–Soil Interface
Conduction across each interface uses:
```
C(i) = thermal_conductivity(i) / (DEPP(i+1) - DEPP(i))
DTDT(i) = (C(i-1)×(T(i-1)-T(i)) + C(i)×(T(i+1)-T(i))) / WC(i)
```
The soil node immediately below the snowpack (index 9 when all 8 snow nodes active) uses this same formula, so heat flows continuously from snow to soil and vice versa — there is no special jump condition at the interface.

### Thermal Properties
Thermal capacitance `WC(i) = density(i) × specific_heat(i) × layer_thickness(i)` and conductance `C(i)` are computed for every node at each timestep. Inactive snow nodes receive zero capacitance and borrow conductance from the first active node below, ensuring they do not create fictitious thermal resistance.

### Snow Albedo Feedback
`daysincesnow` tracks how long the current snowpack has been present. Fresh snow has albedo ~0.9; this decays with age, increasing solar absorption and accelerating melt. A new snowfall event resets `daysincesnow` to zero, temporarily boosting albedo.

---

## 6. Time-Stepping

| Level | Period | Description |
|-------|--------|-------------|
| ODE integration (SFODE) | Adaptive, ~5–60 min | Gear predictor-corrector advances all 18 temperatures |
| Output/snow update (OSUB) | 60 min | Snow mass, melt, phase change, node reassignment |
| Day loop | 24 h | Outer loop; repeats each day `ND` times for spin-up |

Snow mass accounting (snowfall, melt, evaporation) is done once per output interval (hourly), while temperature integration runs at the finer adaptive timestep. This means snow depth changes are stepwise (hourly) while temperatures evolve continuously.

---

## 7. SNOWLAYER Subroutine — Detailed Logic

`SNOWLAYER` is called once per hour (from `OSUB`) after the net snow depth for that hour has been committed to `snowhr(methour)`. Its job is to translate the scalar snow depth into a set of active node depths (`snode`) and a count (`maxsnode1`), and to keep inactive node temperatures consistent with the soil surface.

### Call Context
```
methour = int(SIOUT(1)/60) + 1 + 24×(DOY-1)
```
`methour` is the absolute hour index across the whole simulation. `snowhr(methour)` holds the net snow depth for the current hour, already updated by OSUB before SNOWLAYER is called.

### Step 1 — Hard Cap at 300 cm
```fortran
if (cursnow > 300) maxsnode1 = 0
```
If the snowpack somehow exceeds 300 cm (the deepest node threshold), the node count is reset. This acts as a safety valve; in practice the top node threshold is 300 cm so this condition should never be reached in normal operation.

### Step 2 — Snow Below Minimum Threshold
```fortran
if (cursnow < minsnow) then
    maxsnode  = 0
    maxsnode1 = 0
    snode(1:8) = 0
    snowhr(methour) = 0
    t(1:8)  = t(1)     ! soil surface temperature
    tt(1:8) = tt(1)
    cursnow = 0
    return
end if
```
When snow depth falls below `minsnow` (2 cm), the snowpack is wiped entirely:
- All `snode` depths set to zero (nodes become inactive).
- All snow node temperatures (`T(1:8)` and their previous-step counterparts `TT(1:8)`) are overwritten with the current soil-surface temperature `T(1)` (which, in the 18-node layout, is the topmost soil node). This prevents stale cold or warm temperatures persisting in unused nodes.
- `cursnow` is zeroed and the function returns immediately.

### Step 3 — Determine Number of Active Nodes
When `snowhr(methour) ≥ minsnow`, SNOWLAYER counts how many of the 8 threshold depths are exceeded:
```fortran
do i = 1, 8
    if (snowhr(methour) > snownode(i)) maxsnode = i
end do
if (maxsnode > 7) maxsnode = 7
```
- The loop increments `maxsnode` for each threshold exceeded, so it ends with the index of the deepest threshold still below the current snow depth.
- The cap at 7 ensures at least one "headroom" node above the surface is always available.

### Step 4 — Bottom-Up Node Assignment
```fortran
do i = 1, maxsnode
    snode(i + (8 - maxsnode)) = snownode(i)
end do
```
Nodes are packed into the **bottom** of the 8-element `snode` array. For example, with `maxsnode = 3` (snow between 10 and 20 cm):

| snode index | value |
|-------------|-------|
| 1–5         | 0 (inactive) |
| 6           | 2 cm  |
| 7           | 5 cm  |
| 8           | 10 cm |

This bottom-up convention means the **deepest** node always sits at `snode(8)` and the **shallowest active** node is at `snode(9 - maxsnode)`. Nodes above the snow surface remain at zero.

**Why bottom-up?** The unified depth array `DEPP` places snow nodes at indices 1–8 and soil nodes at 9–18. Packing snow nodes toward index 8 keeps the active nodes adjacent to the soil block (starting at index 9), minimising the number of zero-depth inactive nodes between the topmost snow node and the soil surface.

### Step 5 — Shrinking Snowpack (Node Count Decreases)
```fortran
if (maxsnode < maxsnode1) then
    do i = 1, 8 - maxsnode
        snode(i) = 0
        t(i) = t(1)    ! inherit soil surface temperature
    end do
end if
```
When snow depth drops enough to cross a node threshold downward (one fewer active node than last hour), the newly vacated upper slots are zeroed and their temperatures reset to `T(1)`. This is the same temperature-pinning used in the below-minimum case, applied node by node as the snowpack thins.

### Step 6 — Commit and Return
```fortran
cursnow   = snowhr(methour)
maxsnode1 = real(maxsnode, 8)
```
`cursnow` and `maxsnode1` (stored as double precision in the COMMON block) are updated to reflect the current hour's values. All downstream code (DSUB for thermal properties and ODE derivatives, OSUB for melt and phase change) reads from `maxsnode1`, `snode`, and `cursnow` after SNOWLAYER has returned.

### Interaction with the Rest of the Scheme

| What changes | Effect |
|---|---|
| `snode(1:8)` updated | DSUB rebuilds `DEPP`, `WC`, and `C` using the new node depths at the next integration step |
| `maxsnode1` updated | OSUB uses this to loop over active layers in melt and phase-change calculations |
| `cursnow` updated | OSUB uses this for density evolution, melt depth accounting, and to decide whether to compute `QFREZE` |
| Inactive node temps set to `T(1)` | Prevents zero-depth nodes from introducing spurious temperature gradients in the conduction scheme |

---

## 8. Key Parameters

| Parameter | Description |
|-----------|-------------|
| `snowtemp` | Air temperature threshold (°C) for rain/snow partitioning |
| `snowdens` | Fresh snow density (g/cm³) |
| `snowmelt` | Melt coefficient (0–1); scales energy-balance melt |
| `undercatch` | Gauge undercatch correction (fraction) |
| `intercept` | Vegetation interception fraction |
| `rainmeltf` | Rain-on-snow enhancement factor |
| `snowcond` | Snow thermal conductivity (cal/cm·min·°C) |
| `densfun(4)` | Density evolution parameters (a, b, c, d) |
| `minsnow` | Minimum depth to activate snow model (2 cm, hardcoded) |
| `HTOFN` | Latent heat of fusion (333.55 J/g, hardcoded) |
