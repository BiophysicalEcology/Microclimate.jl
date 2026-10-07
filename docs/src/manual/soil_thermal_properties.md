# Soil thermal properties

Change in soil temperature through depth ``z`` and time ``t`` is described by the
one-dimensional heat equation

```math
\frac{dT}{dt} = \frac{k}{\rho c} \frac{\delta^2 T}{\delta z^2}
```

(Carslaw and Jaeger 1959), which needs one initial condition (a soil temperature profile) and two 
boundary conditions. The initial profile could be based on measurements, or assumed as the mean daily
reference air temperature at all nodes, or — for a continuous multi-day run — the
profile at the end of the previous day; see [Time modes](diel_time_handling.md#Time-modes).
The deep soil boundary condition is the user-supplied deep soil temperature. This is typically
the mean annual air temperature for a maximum profile depth of around 2 m, where soil temperature 
tends towards the running mean annual shaded 2 m air temperature. The surface boundary condition comes from 
equating the energy conducted to the soil surface to the net heat transfer to the surface by
solar radiation, infrared radiation, convection and evaporation:

```math
Q_{cond} = -k_{so} \left.\frac{\delta T}{\delta z}\right|_{z_0} = Q_{solar} + Q_{IR} + Q_{conv} - Q_{evap}
```

## Solving the heat equation

[`SoilHeatTransport1D`](@ref) solves this per-node, as a system of ODEs — one per
depth node, the deepest pinned at the deep soil temperature as a Dirichlet boundary
condition. The system is handed to any [SciML](https://sciml.ai/) ODE solver
(`Tsit5()` by default, configurable via `ode_solver`), with `ode_kwargs` controlling
solver tolerances. For a given node ``i``, the discretized right-hand side is

```math
\frac{dT_i}{dt} = \frac{K_{i-1}(T_{i-1}-T_i) + K_i(T_{i+1}-T_i)}{C_i},
```

with node heat capacity ``C_i = \rho_i c_i \frac{z_{i+1}-z_{i-1}}{2}`` and
``K_i = k_i (z_{i+1}-z_i)``, ``k_i`` the thermal conductivity. The depth nodes
themselves are `MicroModel.depths` — any count or spacing the user chooses, 19 nodes
by default (`DEFAULT_DEPTHS`), spaced closer near the surface where
temperature changes fastest — see the soil-depth-grid note in
[Soil hydraulics algorithms](../tutorials/soil_hydraulics_algorithms.md) for why
near-surface spacing matters for both accuracy and solver stability. Be sure not to
make the near-surface spacing too short, however, to avoid numerical problems; 
typically ~1-2 cm is sufficient between the first two nodes.

## Soil thermal properties from bulk composition

Bulk density ``\rho``, specific heat capacity ``c_s`` and thermal conductivity ``k_s``
all vary with soil moisture and temperature, and with depth as a function of soil
composition. Each layer can be independently configured with differing soil properties. 
[`CampbelldeVriesSoilProperties`](@ref) computes them from the mineral
composition on [`SoilProfile`](@ref) (bulk density, mineral density, mineral
conductivity, mineral heat capacity) plus a de Vries shape factor, following Campbell
et al. (1994, eqs. 8, 9) and Campbell and Norman (1998, eqs. 8.13, 8.17, 8.20, and the
soil-texture properties in their table 9.1 — see [`example_soil_properties_model`](@ref)
for the mineral-fraction defaults). This computation has no phase-transition term of
its own — soil freeze/thaw is handled separately, as a correction applied to the ODE's
result once per hour rather than folded into the specific heat, described in
[Phase transition](phase_transition.md).

## References

See the [References](references.md) page for the full bibliography.
