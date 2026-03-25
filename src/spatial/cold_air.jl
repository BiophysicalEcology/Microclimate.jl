"""
Cold air drainage models for microclimate simulation.

Two approaches are available:

1. `KLAM21` - Dynamic model from Sievers & Kossmann (2016)
   Simulates time-evolution of heat deficit, cold air depth, and velocity.

2. `MacleanColdAirDrainage` - Static potential model from Maclean et al. (2019)
   Diagnostic approach using flow accumulation and lapse rate.

# References
- Sievers & Kossmann (2016) "The cold air drainage model KLAM_21", Weather and Climate, 36, 2-24.
- Maclean et al. (2019) "Microclima: An R package for modelling meso- and microclimate", Methods Ecol Evol, 10, 280-290.
"""

using Unitful: Unitful, @u_str, Quantity, ustrip
using Unitful: m, s, K, J, kg, W
using Stencils: Stencils, VonNeumann, StencilArray, mapstencil!, neighbors, center, Remove
using KernelAbstractions: KernelAbstractions, @kernel, @index, get_backend

# =============================================================================
# KernelAbstractions Kernels for Pointwise Operations
# =============================================================================

"""
Kernel to clamp velocity components to [-max_vel, max_vel].
"""
@kernel function _clamp_velocity_kernel!(u, v, max_vel)
    i, j = @index(Global, NTuple)
    @inbounds u[i, j] = clamp(u[i, j], -max_vel, max_vel)
    @inbounds v[i, j] = clamp(v[i, j], -max_vel, max_vel)
end

"""
Kernel to diagnose cold air depth and temperature disturbance from heat deficit.

Uses ParabolicProfile equations:
- H = H₀ × (E / (ρ₀ cₚ f H₀ ΔT₀))^(2/3)
- ΔT = ΔT₀ × √(H / H₀)
"""
@kernel function _diagnose_depth_kernel!(
    depth, temperature_disturbance, heat_deficit,
    H₀, ΔT₀, ρ₀, cₚ, f, max_depth, zero_depth, zero_temp
)
    i, j = @index(Global, NTuple)
    @inbounds begin
        E = heat_deficit[i, j]
        if E <= zero(E)
            depth[i, j] = zero_depth
            temperature_disturbance[i, j] = zero_temp
        else
            # H = H₀ × (E / (ρ₀ cₚ f H₀ ΔT₀))^(2/3)
            H = H₀ * (E / (ρ₀ * cₚ * f * H₀ * ΔT₀))^(2/3)
            H = min(H, max_depth)
            depth[i, j] = H

            # ΔT = ΔT₀ × √(H / H₀)
            temperature_disturbance[i, j] = ΔT₀ * sqrt(H / H₀)
        end
    end
end

# =============================================================================
# Abstract Type
# =============================================================================

"""
    ColdAirModel

Abstract type for cold air drainage models.

Subtypes:
- `KLAM21`: Dynamic model with heat deficit evolution
- `MacleanColdAirDrainage`: Static potential-based model
"""
abstract type ColdAirModel end

# =============================================================================
# Lapse Rate - temperature change with elevation
# =============================================================================

"""
    LapseRateMethod

Abstract type for computing environmental lapse rate.

Subtypes:
- `FixedLapseRate`: Constant lapse rate (default 6.5 K/km)
- `HessLapseRate`: Dynamic moist adiabatic lapse rate from Hess (1959)
"""
abstract type LapseRateMethod end

"""
    FixedLapseRate <: LapseRateMethod

Fixed environmental lapse rate.

Default is 6.5 K/km, the standard atmospheric lapse rate.
"""
Base.@kwdef struct FixedLapseRate{L} <: LapseRateMethod
    rate::L = 0.0065u"K/m"
end

"""
    HessLapseRate <: LapseRateMethod

Dynamic lapse rate from Hess (1959).

Computes moist adiabatic lapse rate from temperature and humidity.
The moist adiabatic rate is lower than dry adiabatic (~10 K/km) because
latent heat release during condensation partially offsets cooling.

Typical values:
- Saturated warm air: ~3-4 K/km
- Drier conditions: ~6-8 K/km
"""
struct HessLapseRate <: LapseRateMethod end

"""
    lapse_rate(method::FixedLapseRate; kwargs...)

Return fixed lapse rate.
"""
lapse_rate(method::FixedLapseRate; kwargs...) = method.rate

"""
    lapse_rate(method::HessLapseRate; temperature, mixing_ratio)

Compute environmental lapse rate using Hess (1959) formula.

    lapse_rate = g × (1 + Lv×rv/(R×T)) / (cpd + Lv²×rv/(R×T²))

where:
- g = gravitational acceleration (9.81 m/s²)
- Lv = latent heat of vaporization (2.5 MJ/kg)
- R = gas constant for dry air (287 J/kg/K)
- cpd = specific heat of dry air (1004 J/kg/K)
- T = temperature (K)
- rv = mixing ratio (kg water / kg dry air)

# Arguments
- `temperature`: Air temperature (K)
- `mixing_ratio`: Water vapor mixing ratio (dimensionless, kg/kg)

# Returns
Lapse rate (K/m)
"""
function lapse_rate(method::HessLapseRate; temperature, mixing_ratio)
    g = 9.8076u"m/s^2"
    Lv = 2501000.0u"J/kg"      # Latent heat of vaporization
    R = 287.0u"J/kg/K"         # Gas constant for dry air
    cpd = 1003.5u"J/kg/K"      # Specific heat of dry air

    T = temperature
    rv = mixing_ratio

    numerator = 1 + (Lv * rv) / (R * T)
    denominator = cpd + (Lv^2 * rv) / (R * T^2)

    return g * numerator / denominator
end

# =============================================================================
# Temperature Profile - how temperature varies with depth in cold air
# =============================================================================

"""
    TemperatureProfile

Abstract type for computing temperature disturbance from cold air depth.

Both cold air models compute temperature disturbance from depth, but use
different formulas:
- KLAM21 uses a parabolic profile (square root relationship)
- Maclean uses a linear profile (constant lapse rate)

The "depth" concept is:
- For KLAM21: computed cold air layer thickness from heat deficit
- For Maclean: elevation difference below basin maximum
"""
abstract type TemperatureProfile end

"""
    ParabolicProfile <: TemperatureProfile

KLAM21's parabolic temperature profile from Sievers & Kossmann (2016).

Temperature disturbance increases with square root of depth:
    temperature_disturbance = reference_temperature_disturbance × √(depth / reference_depth)

This comes from integrating the parabolic temperature profile T'(z) = ΔT × ((H-z)/H)²
"""
Base.@kwdef struct ParabolicProfile{D,T} <: TemperatureProfile
    reference_depth::D = 10.0u"m"
    reference_temperature_disturbance::T = 3.0u"K"
end

"""
    LinearProfile <: TemperatureProfile

Linear temperature profile using environmental lapse rate.

Temperature disturbance is proportional to depth:
    temperature_disturbance = lapse_rate × depth

The lapse rate can be fixed or computed dynamically (e.g., HessLapseRate).
"""
Base.@kwdef struct LinearProfile{L<:LapseRateMethod} <: TemperatureProfile
    lapse_rate_method::L = FixedLapseRate()
end


"""
    temperature_disturbance(model::ColdAirModel, depth; kwargs...)

Calculate surface temperature disturbance from cold air depth.

Delegates to `temperature_disturbance(model.temperature_profile, depth; kwargs...)`.
Works for any ColdAirModel subtype (KLAM21, MacleanColdAirDrainage, etc.)

# Returns
Temperature disturbance (K), positive = colder than ambient.
"""
function temperature_disturbance(model::ColdAirModel, depth; kwargs...)
    return temperature_disturbance(model.temperature_profile, depth; kwargs...)
end
"""
    temperature_disturbance(profile::ParabolicProfile, depth; kwargs...)

Compute temperature disturbance using KLAM21's parabolic profile.

Returns temperature disturbance (K), positive = colder than ambient.
"""
function temperature_disturbance(profile::ParabolicProfile, depth; kwargs...)
    depth <= zero(depth) && return zero(profile.reference_temperature_disturbance)
    return profile.reference_temperature_disturbance * sqrt(depth / profile.reference_depth)
end
"""
    temperature_disturbance(profile::LinearProfile, depth; temperature=288.0u"K",
                            mixing_ratio=0.01, kwargs...)

Compute temperature disturbance using linear lapse rate profile.

For HessLapseRate, requires temperature and mixing_ratio to compute the lapse rate.
For FixedLapseRate, these are ignored.

Returns temperature disturbance (K), positive = colder than ambient.
"""
function temperature_disturbance(
    profile::LinearProfile,
    depth;
    temperature = 288.0u"K",
    mixing_ratio = 0.01,
    kwargs...
)
    depth <= zero(depth) && return 0.0u"K"
    Γ = lapse_rate(profile.lapse_rate_method; temperature, mixing_ratio)
    return Γ * depth
end

# =============================================================================
# Drainage Conditions - when does drainage occur?
# =============================================================================

"""
    DrainageConditions

Abstract type for determining when cold air drainage occurs.

Subtypes:
- `AlwaysDraining`: Drainage always occurs (for dynamic models like KLAM21)
- `KlokOerlemansConditions`: Based on wind, sky emissivity, and time
- `SimpleWindConditions`: Based on wind speed only
"""
abstract type DrainageConditions end

"""
    AlwaysDraining <: DrainageConditions

Drainage is always active. Used by dynamic models like KLAM21 that
run continuously and let physics determine drainage intensity.
"""
struct AlwaysDraining <: DrainageConditions end

"""
    KlokOerlemansConditions <: DrainageConditions

Drainage conditions based on Klok & Oerlemans (2002) sky emissivity.

Used by Maclean et al. (2019). Drainage occurs when:
- Wind speed below threshold
- Sky emissivity below threshold (clear sky)
- Nighttime or within post-dawn hours
- Consecutive favorable hours met
"""
Base.@kwdef struct KlokOerlemansConditions{W,E} <: DrainageConditions
    wind_speed_threshold::W = 4.5u"m/s"
    emissivity_threshold::E = 0.5
    post_dawn_hours::Float64 = 3.0
    consecutive_hours_required::Int = 3
end

"""
    SimpleWindConditions <: DrainageConditions

Simple wind-only drainage conditions.

Drainage occurs when wind speed is below threshold at night.
"""
Base.@kwdef struct SimpleWindConditions{W} <: DrainageConditions
    wind_speed_threshold::W = 4.5u"m/s"
    post_dawn_hours::Float64 = 3.0
end

# =============================================================================
# Flow Weighting - how flow accumulation affects temperature
# =============================================================================

"""
    FlowWeighting

Abstract type for weighting temperature by flow accumulation.
"""
abstract type FlowWeighting end

"""
    LogFlowWeighting <: FlowWeighting

Logarithmic flow weighting from Maclean et al. (2019).

Temperature disturbance is scaled by log(F) where F is normalized flow accumulation.
"""
struct LogFlowWeighting <: FlowWeighting end

"""
    LinearFlowWeighting <: FlowWeighting

Linear flow weighting.

Temperature disturbance is scaled directly by normalized flow accumulation.
"""
struct LinearFlowWeighting <: FlowWeighting end

"""
    NoFlowWeighting <: FlowWeighting

No flow weighting - temperature depends only on elevation.
"""
struct NoFlowWeighting <: FlowWeighting end

# =============================================================================
# KLAM21 Model (Dynamic)
# =============================================================================

"""
    KLAM21 <: ColdAirModel

KLAM_21 cold air drainage model parameters.

Based on Sievers & Kossmann (2016). The model calculates depth, velocity, and
direction of a stably stratified boundary layer (cold air layer) that evolves
during nighttime over complex terrain.

# Composable Components
- `conditions`: When drainage occurs (default: AlwaysDraining)
- `temperature_profile`: How temperature varies with depth (default: ParabolicProfile)

# Model Parameters
- `effective_fraction`: Fraction β of cold air layer for gravitational forcing (dimensionless)
- `von_karman_constant`: von Karman constant κ (dimensionless)
- `canopy_drag_coefficient`: Drag coefficient for canopy friction (dimensionless)
- `mixing_length`: Mixing length for horizontal diffusion (m)

# Physical Constants
- `gravitational_acceleration`: g (m/s²)
- `heat_capacity`: Heat capacity of air at constant pressure cₚ (J/kg/K)
- `reference_density`: Reference air density ρ₀ (kg/m³)
- `reference_temperature`: Reference temperature T₀ (K)

# References
- Sievers, U. (2005): Das Kaltluft-Abfluss-Modell KLAM_21. Berichte des DWD 227.
- Sievers & Kossmann (2016): Weather and Climate, 36, 2-24.
"""
Base.@kwdef struct KLAM21{C<:DrainageConditions,P<:TemperatureProfile,L,T,A,H,D} <: ColdAirModel
    # Composable components
    conditions::C = AlwaysDraining()
    temperature_profile::P = ParabolicProfile()
    # Model parameters
    effective_fraction::Float64 = 5/12
    von_karman_constant::Float64 = 0.4
    canopy_drag_coefficient::Float64 = 0.2
    mixing_length::L = 1.0u"m"
    # Physical constants
    gravitational_acceleration::A = 9.81u"m/s^2"
    heat_capacity::H = 1005.0u"J/kg/K"
    reference_density::D = 1.2u"kg/m^3"
    reference_temperature::T = 288.0u"K"
end

# Derived constant: f = 1/3 for parabolic temperature profile
const KLAM21_PROFILE_FACTOR = 1/3

# =============================================================================
# KA Kernel Wrapper Functions
# =============================================================================

"""
    clamp_velocity!(state; max_velocity)

Clamp velocity components using KernelAbstractions for GPU compatibility.
"""
function clamp_velocity!(state; max_velocity)
    u = state.velocity_x
    v = state.velocity_y
    backend = get_backend(u)
    kernel = _clamp_velocity_kernel!(backend)
    kernel(u, v, max_velocity; ndrange=size(u))
    return state
end

"""
    diagnose_depth!(model::KLAM21, state; max_depth)

Diagnose cold air depth and temperature disturbance from heat deficit.
Uses KernelAbstractions for GPU compatibility.

Requires model.temperature_profile to be a ParabolicProfile.
"""
function diagnose_depth!(model::KLAM21, state; max_depth)
    profile = model.temperature_profile
    profile isa ParabolicProfile || error("diagnose_depth! requires ParabolicProfile")

    H₀ = profile.reference_depth
    ΔT₀ = profile.reference_temperature_disturbance
    ρ₀ = model.reference_density
    cₚ = model.heat_capacity
    f = KLAM21_PROFILE_FACTOR

    zero_depth = zero(H₀)
    zero_temp = zero(ΔT₀)

    backend = get_backend(state.heat_deficit)
    kernel = _diagnose_depth_kernel!(backend)
    kernel(
        state.depth, state.temperature_disturbance, state.heat_deficit,
        H₀, ΔT₀, ρ₀, cₚ, f, max_depth, zero_depth, zero_temp;
        ndrange=size(state.heat_deficit)
    )
    return state
end

# =============================================================================
# Diagnostic Functions
# =============================================================================

"""
    cold_air_depth(params::KLAM21, heat_deficit; volume_reduction_factor=1.0)

Calculate cold air layer depth from heat deficit.

Implements Equation 4 from Sievers & Kossmann (2016):
    H = H₀ × (E / (rᵥ ρ₀ cₚ f H₀ ΔT₀))^(2/3)

Requires `params.temperature_profile` to be a `ParabolicProfile` (the default).

# Arguments
- `params`: KLAM21 parameters
- `heat_deficit`: Heat deficit E (J/m²)

# Keywords
- `volume_reduction_factor`: rᵥ, reduction factor (< 1 for urban areas)

# Returns
Cold air layer depth (m)
"""
function cold_air_depth(params::KLAM21, heat_deficit; volume_reduction_factor=1.0)
    E = heat_deficit
    E <= 0.0u"J/m^2" && return 0.0u"m"

    profile = params.temperature_profile
    profile isa ParabolicProfile || error("cold_air_depth requires ParabolicProfile")

    # Short names for math
    H₀ = profile.reference_depth
    ΔT₀ = profile.reference_temperature_disturbance
    ρ₀ = params.reference_density
    cₚ = params.heat_capacity
    rᵥ = volume_reduction_factor
    f = KLAM21_PROFILE_FACTOR

    # Eq 4: H = H₀ × (E / (rᵥ ρ₀ cₚ f H₀ ΔT₀))^(2/3)
    return H₀ * (E / (rᵥ * ρ₀ * cₚ * f * H₀ * ΔT₀))^(2/3)
end

"""
    effective_cold_air_depth(params::KLAM21, depth)

Calculate effective cold air layer depth: H_eff = β × H.

Only the lower fraction β of the cold air layer contributes to drainage flow.
"""
function effective_cold_air_depth(params::KLAM21, depth)
    return params.effective_fraction * depth
end

# =============================================================================
# Friction Parameterization (Appendix A)
# =============================================================================

"""
    friction_coefficient(params::KLAM21, effective_depth, roughness_length)

Calculate friction coefficient for surface friction.

Implements Equation A1 from Sievers & Kossmann (2016):
    c* = (2κ / ln(0.25 H_eff / z₀))²

# Arguments
- `params`: KLAM21 parameters
- `effective_depth`: Effective cold air layer depth (m)
- `roughness_length`: Aerodynamic roughness length z₀ (m)

# Returns
Friction coefficient c* (dimensionless)
"""
function friction_coefficient(params::KLAM21, effective_depth, roughness_length)
    effective_depth <= 0.0u"m" && return 0.0

    κ = params.von_karman_constant
    height_ratio = 0.25 * effective_depth / roughness_length  # m/m = dimensionless
    height_ratio <= 1 && return 1.0  # Maximum friction when layer is very shallow

    return (2 * κ / log(height_ratio))^2
end

"""
    friction_coefficient_with_canopy(params::KLAM21, effective_depth, roughness_length,
                                      canopy_height, canopy_density)

Calculate friction coefficient including canopy effects.

Implements Equation A2 when wind maximum is within canopy:
    c* = (2κ / ln(0.25 H_eff / z₀))² + c_d × σ × min(H, h_canopy)

# Arguments
- `params`: KLAM21 parameters
- `effective_depth`: Effective cold air layer depth (m)
- `roughness_length`: Aerodynamic roughness length z₀ (m)
- `canopy_height`: Canopy height (m)
- `canopy_density`: Canopy surface density σ (1/m, from LAI or WAI)
"""
function friction_coefficient_with_canopy(
    params::KLAM21,
    effective_depth,
    roughness_length,
    canopy_height,
    canopy_density,
)
    base_friction = friction_coefficient(params, effective_depth, roughness_length)

    β = params.effective_fraction
    c_d = params.canopy_drag_coefficient

    # Add canopy friction if wind maximum is within canopy
    if 0.25 * effective_depth <= canopy_height
        # c_d is dimensionless, σ has units 1/m, height has units m → dimensionless
        canopy_friction = c_d * canopy_density * min(effective_depth / β, canopy_height)
        return base_friction + canopy_friction
    end

    return base_friction
end

# =============================================================================
# Evolution Equations
# =============================================================================

"""
    update_heat_deficit!(params::KLAM21, state, heat_loss_rate, elevation;
                         timestep, cellsize, effective_density=nothing)

Update heat deficit for one timestep.

Implements Equation 1:
    ∂E/∂t = P - ∇·(ρ_eff cₚ H ⟨vₕ T'⟩)

The advection term uses a donor-cell (upwind) scheme.

# Arguments
- `params`: KLAM21 parameters
- `state`: KLAM21State with heat_deficit, depth, temperature_disturbance, velocity_x, velocity_y
- `heat_loss_rate`: Local heat loss rate P (W/m²), from energy balance
- `elevation`: Digital elevation model (m)

# Keywords
- `timestep`: Timestep (s or any time unit)
- `cellsize`: Grid cell size (dx, dy) with length units
- `effective_density`: Effective density field (defaults to reference_density)
"""
function update_heat_deficit!(
    params::KLAM21,
    state,
    heat_loss_rate::AbstractMatrix,
    elevation::AbstractMatrix;
    timestep,
    cellsize,
    effective_density::Union{Nothing,AbstractMatrix} = nothing,
)
    E = state.heat_deficit
    H = state.depth
    ΔT = state.temperature_disturbance
    u = state.velocity_x
    v = state.velocity_y

    cₚ = params.heat_capacity
    ρ₀ = params.reference_density
    dx = abs(cellsize[1])
    dy = abs(cellsize[2])

    zero_energy = 0.0u"J/m^2"
    zero_vel = 0.0u"m/s"

    # Set up stencil arrays
    hood = VonNeumann{1}()
    E_sa = StencilArray(E, hood; boundary = Remove(zero_energy))
    H_sa = StencilArray(H, hood; boundary = Remove(zero(eltype(H))))
    ΔT_sa = StencilArray(ΔT, hood; boundary = Remove(zero(eltype(ΔT))))
    u_sa = StencilArray(u, hood; boundary = Remove(zero_vel))
    v_sa = StencilArray(v, hood; boundary = Remove(zero_vel))
    P_sa = StencilArray(heat_loss_rate, hood; boundary = Remove(zero(eltype(heat_loss_rate))))

    # Double buffer
    E_next = similar(E)

    # Compute new heat deficit using upwind stencil
    mapstencil!(E_next, E_sa, H_sa, ΔT_sa, u_sa, v_sa, P_sa) do E_hood, H_hood, ΔT_hood, u_hood, v_hood, P_hood
        E_c = center(E_hood)
        H_c = center(H_hood)
        ΔT_c = center(ΔT_hood)
        u_c = center(u_hood)
        v_c = center(v_hood)
        P_c = center(P_hood)

        # VonNeumann{1} neighbors: 1=N(0,1), 2=S(0,-1), 3=E(1,0), 4=W(-1,0)
        H_n = neighbors(H_hood)
        ΔT_n = neighbors(ΔT_hood)
        u_n = neighbors(u_hood)
        v_n = neighbors(v_hood)

        # Heat content at center: ρ cₚ H ΔT
        heat_content = ρ₀ * cₚ * H_c * ΔT_c

        # x-direction flux divergence (upwind scheme)
        # If u > 0: flow from W (index 4), out to E
        # If u < 0: flow from E (index 3), out to W
        if u_c > zero_vel
            flux_x_in = ρ₀ * cₚ * H_n[4] * ΔT_n[4] * u_n[4]  # W neighbor
            flux_x_out = heat_content * u_c
        else
            flux_x_in = ρ₀ * cₚ * H_n[3] * ΔT_n[3] * (-u_n[3])  # E neighbor
            flux_x_out = heat_content * (-u_c)
        end
        div_x = (flux_x_out - flux_x_in) / dx

        # y-direction flux divergence (upwind scheme)
        # If v > 0: flow from S (index 2), out to N
        # If v < 0: flow from N (index 1), out to S
        if v_c > zero_vel
            flux_y_in = ρ₀ * cₚ * H_n[2] * ΔT_n[2] * v_n[2]  # S neighbor
            flux_y_out = heat_content * v_c
        else
            flux_y_in = ρ₀ * cₚ * H_n[1] * ΔT_n[1] * (-v_n[1])  # N neighbor
            flux_y_out = heat_content * (-v_c)
        end
        div_y = (flux_y_out - flux_y_in) / dy

        # Update: ∂E/∂t = P - ∇·(flux)
        return max(zero_energy, E_c + timestep * (P_c - div_x - div_y))
    end

    # Copy results back
    state.heat_deficit .= E_next

    return state
end

"""
    update_momentum!(params::KLAM21, state, elevation, roughness_length;
                     timestep, cellsize, regional_wind=nothing)

Update velocity components for one timestep.

Implements Equation 5 (simplified, without diffusion):
    ∂vₕ/∂t = g(ΔT/T₀)∇(h₀ + βH) - c* |vₜ| vₜ / H + regional_wind_term

# Arguments
- `params`: KLAM21 parameters
- `state`: KLAM21State with depth, temperature_disturbance, velocity_x, velocity_y
- `elevation`: Digital elevation model (m)
- `roughness_length`: Roughness length (m), scalar or grid

# Keywords
- `timestep`: Timestep (s or any time unit)
- `cellsize`: Grid cell size (dx, dy) with length units
- `regional_wind`: Regional wind as (u, v) tuple (m/s), or nothing
"""
function update_momentum!(
    params::KLAM21,
    state,
    elevation::AbstractMatrix,
    roughness_length;
    timestep,
    cellsize,
    regional_wind = nothing,
)
    H = state.depth
    ΔT = state.temperature_disturbance
    u = state.velocity_x
    v = state.velocity_y

    g = params.gravitational_acceleration
    T₀ = params.reference_temperature
    β = params.effective_fraction
    l = params.mixing_length
    dx = abs(cellsize[1])
    dy = abs(cellsize[2])

    min_depth = 0.1u"m"
    ν_desired = l^2 / (1.0u"s")
    ν_stable = min(dx, dy)^2 / (4 * timestep)
    ν = min(ν_desired, ν_stable)

    # Regional wind components (scalars for now, grids would need StencilArray)
    u_reg = isnothing(regional_wind) ? nothing : regional_wind[1]
    v_reg = isnothing(regional_wind) ? nothing : regional_wind[2]
    exchange_rate = 0.001u"s^-1"

    # Set up stencil arrays with VonNeumann neighborhood (4 neighbors + center)
    hood = VonNeumann{1}()
    boundary = Remove(zero(eltype(u)))

    u_sa = StencilArray(u, hood; boundary)
    v_sa = StencilArray(v, hood; boundary)
    H_sa = StencilArray(H, hood; boundary = Remove(zero(eltype(H))))
    ΔT_sa = StencilArray(ΔT, hood; boundary = Remove(zero(eltype(ΔT))))
    elev_sa = StencilArray(elevation, hood; boundary = Remove(zero(eltype(elevation))))
    z₀_sa = roughness_length isa AbstractMatrix ?
            StencilArray(roughness_length, hood; boundary = Remove(zero(eltype(roughness_length)))) :
            nothing

    # Double buffer for output
    u_next = similar(u)
    v_next = similar(v)

    # Compute new u velocity using stencil operation
    mapstencil!(u_next, u_sa, v_sa, H_sa, ΔT_sa, elev_sa) do u_hood, v_hood, H_hood, ΔT_hood, elev_hood
        H_c = center(H_hood)
        H_c <= zero(H_c) && return center(u_hood)

        H_eff = β * H_c
        u_c = center(u_hood)
        v_c = center(v_hood)
        ΔT_c = center(ΔT_hood)

        # Get neighbors: VonNeumann order is (N, S, E, W) = indices 1,2,3,4
        u_n = neighbors(u_hood)
        H_n = neighbors(H_hood)
        elev_n = neighbors(elev_hood)

        # Gradient of effective surface (central difference)
        # VonNeumann{1} neighbors: 1=N(0,1), 2=S(0,-1), 3=E(1,0), 4=W(-1,0)
        eff_E = elev_n[3] + β * H_n[3]
        eff_W = elev_n[4] + β * H_n[4]
        grad_x = (eff_E - eff_W) / (2 * dx)

        # Buoyancy forcing
        buoyancy = g * ΔT_c / T₀
        F_gravity_x = -buoyancy * grad_x

        # Friction
        c_star = _friction_coefficient_inner(params.von_karman_constant, H_eff, roughness_length)
        speed = sqrt(u_c^2 + v_c^2)
        H_safe = max(H_eff, min_depth)
        F_friction_x = -c_star * speed * u_c / H_safe

        # Laplacian (diffusion)
        laplacian_u = (u_n[3] - 2*u_c + u_n[4]) / dx^2 + (u_n[1] - 2*u_c + u_n[2]) / dy^2
        diff_u = ν * laplacian_u

        # Regional wind coupling
        F_regional_x = isnothing(u_reg) ? zero(F_gravity_x) : exchange_rate * (u_reg - u_c)

        # New u velocity
        return u_c + timestep * (F_gravity_x + F_friction_x + diff_u + F_regional_x)
    end

    # Compute new v velocity using stencil operation
    mapstencil!(v_next, u_sa, v_sa, H_sa, ΔT_sa, elev_sa) do u_hood, v_hood, H_hood, ΔT_hood, elev_hood
        H_c = center(H_hood)
        H_c <= zero(H_c) && return center(v_hood)

        H_eff = β * H_c
        u_c = center(u_hood)
        v_c = center(v_hood)
        ΔT_c = center(ΔT_hood)

        # Get neighbors
        v_n = neighbors(v_hood)
        H_n = neighbors(H_hood)
        elev_n = neighbors(elev_hood)

        # Gradient of effective surface (central difference)
        eff_N = elev_n[1] + β * H_n[1]
        eff_S = elev_n[2] + β * H_n[2]
        grad_y = (eff_N - eff_S) / (2 * dy)

        # Buoyancy forcing
        buoyancy = g * ΔT_c / T₀
        F_gravity_y = -buoyancy * grad_y

        # Friction
        c_star = _friction_coefficient_inner(params.von_karman_constant, H_eff, roughness_length)
        speed = sqrt(u_c^2 + v_c^2)
        H_safe = max(H_eff, min_depth)
        F_friction_y = -c_star * speed * v_c / H_safe

        # Laplacian (diffusion)
        laplacian_v = (v_n[3] - 2*v_c + v_n[4]) / dx^2 + (v_n[1] - 2*v_c + v_n[2]) / dy^2
        diff_v = ν * laplacian_v

        # Regional wind coupling
        F_regional_y = isnothing(v_reg) ? zero(F_gravity_y) : exchange_rate * (v_reg - v_c)

        # New v velocity
        return v_c + timestep * (F_gravity_y + F_friction_y + diff_v + F_regional_y)
    end

    # Copy results back to state
    state.velocity_x .= u_next
    state.velocity_y .= v_next

    return state
end

# Inner friction calculation without matrix indexing (for use in stencil)
function _friction_coefficient_inner(κ, effective_depth, roughness_length)
    effective_depth <= zero(effective_depth) && return 0.0
    z₀ = roughness_length isa AbstractMatrix ? roughness_length : roughness_length
    height_ratio = 0.25 * effective_depth / z₀
    height_ratio <= 1 && return 1.0
    return (2 * κ / log(ustrip(height_ratio)))^2
end

# =============================================================================
# Main Timestep Function - Unified Interface
# =============================================================================

"""
    cold_air_step!(model::ColdAirModel, state, env::SpatialEnvironment; kwargs...)

Perform one timestep of a cold air drainage model.

Unified interface for all ColdAirModel subtypes. Dispatches to model-specific
implementations while ensuring consistent handling of:
- Drainage conditions checking
- Environmental parameter passing to temperature_disturbance
- Solver control parameters

See method-specific docstrings for details.
"""
function cold_air_step! end

"""
    cold_air_step!(model::KLAM21, state::KLAM21State, env::SpatialEnvironment;
                   timestep, heat_loss_rate,
                   hours_since_sunset=0.0, hours_since_sunrise=12.0,
                   max_velocity=50.0u"m/s", max_depth=500.0u"m")

Perform one timestep of the KLAM_21 model.

Updates heat deficit, diagnoses depth and temperature disturbance,
and updates velocity components. Checks drainage conditions before running.

# Arguments
- `model`: KLAM21 parameters
- `state`: KLAM21State with fields: heat_deficit, depth, temperature_disturbance, velocity_x, velocity_y
- `env`: SpatialEnvironment containing terrain, weather, land_surface

# Keywords
- `timestep`: Timestep (s, min, hr - any time unit)
- `heat_loss_rate`: Local heat loss rate P (W/m²), from energy balance
- `hours_since_sunset`: Hours since sunset (for conditions check)
- `hours_since_sunrise`: Hours since sunrise (for conditions check)
- `max_velocity`: Maximum velocity for numerical stability
- `max_depth`: Maximum depth for numerical stability
"""
function cold_air_step!(
    model::KLAM21,
    state::KLAM21State,
    env::SpatialEnvironment;
    timestep,
    heat_loss_rate::AbstractMatrix,
    hours_since_sunset = 0.0,
    hours_since_sunrise = 12.0,
    max_velocity = 50.0u"m/s",
    max_depth = 500.0u"m",
)
    nx, ny = size(state.heat_deficit)
    elevation = env.terrain.dem
    roughness_length = env.land_surface.roughness_length
    cellsize = env.terrain.cellsize
    wind = regional_wind(env)

    # Get environmental values for conditions check and temperature profile
    # Use scalar or center value for uniform check
    T_ref = env.weather.air_temperature isa AbstractMatrix ?
            env.weather.air_temperature[nx÷2, ny÷2] : env.weather.air_temperature
    ws_ref = env.weather.wind_speed isa AbstractMatrix ?
             env.weather.wind_speed[nx÷2, ny÷2] : env.weather.wind_speed
    ea_ref = vapor_pressure(env.weather.air_temperature isa AbstractMatrix ?
             env.weather.air_temperature[nx÷2, ny÷2] : env.weather.air_temperature,
             env.weather.humidity isa AbstractMatrix ?
             env.weather.humidity[nx÷2, ny÷2] : env.weather.humidity)
    cf_ref = env.weather.cloud_cover isa AbstractMatrix ?
             env.weather.cloud_cover[nx÷2, ny÷2] : env.weather.cloud_cover
    P_ref = env.weather.pressure isa AbstractMatrix ?
            env.weather.pressure[nx÷2, ny÷2] : env.weather.pressure

    # Check drainage conditions
    ε = sky_emissivity(T_ref, ea_ref, cf_ref)
    favorable = is_drainage_favorable(
        model.conditions;
        wind_speed = ws_ref,
        emissivity = ε,
        hours_since_sunset = hours_since_sunset,
        hours_since_sunrise = hours_since_sunrise,
    )

    if !favorable
        # No drainage - just decay existing cold air
        return state
    end

    # 1. Update momentum (velocity)
    update_momentum!(model, state, elevation, roughness_length;
                     timestep, cellsize, regional_wind = wind)

    # 2. Clamp velocities for numerical stability (KA kernel)
    clamp_velocity!(state; max_velocity)

    # 3. Update heat deficit
    update_heat_deficit!(model, state, heat_loss_rate, elevation;
                         timestep, cellsize)

    # 4. Diagnose depth and temperature disturbance from updated heat deficit (KA kernel)
    diagnose_depth!(model, state; max_depth)

    return state
end

# =============================================================================
# Helper Functions
# =============================================================================

"""
    is_drainage_favorable(cond::AlwaysDraining; kwargs...)

Always returns true - drainage is always active.
"""
is_drainage_favorable(::AlwaysDraining; kwargs...) = true

"""
    is_drainage_favorable(cond::KlokOerlemansConditions; wind_speed, emissivity,
                          hours_since_sunset, hours_since_sunrise)

Check if conditions favor drainage using Klok & Oerlemans criteria.
"""
function is_drainage_favorable(
    cond::KlokOerlemansConditions;
    wind_speed,
    emissivity,
    hours_since_sunset,
    hours_since_sunrise,
)
    # Night or within post-dawn window
    is_night_or_dawn = hours_since_sunset >= 0 || hours_since_sunrise <= cond.post_dawn_hours
    is_night_or_dawn || return false

    wind_speed <= cond.wind_speed_threshold || return false
    emissivity <= cond.emissivity_threshold || return false
    return true
end

"""
    is_drainage_favorable(cond::SimpleWindConditions; wind_speed,
                          hours_since_sunset, hours_since_sunrise, kwargs...)

Check if conditions favor drainage using simple wind threshold.
"""
function is_drainage_favorable(
    cond::SimpleWindConditions;
    wind_speed,
    hours_since_sunset,
    hours_since_sunrise,
    kwargs...
)
    is_night_or_dawn = hours_since_sunset >= 0 || hours_since_sunrise <= cond.post_dawn_hours
    is_night_or_dawn || return false

    wind_speed <= cond.wind_speed_threshold || return false
    return true
end

"""
    flow_weight(::LogFlowWeighting, normalized_flow)

Compute log flow weight (Maclean et al. 2019).

Returns -log(F), which is positive for F < 1. This matches the paper's
equation 6 which has a leading negative sign: ΔTK = −Γm × Δzm × log(F).
"""
function flow_weight(::LogFlowWeighting, normalized_flow)
    F = max(normalized_flow, 0.01)  # Avoid log(0)
    return -log(F)
end

"""
    flow_weight(::LinearFlowWeighting, normalized_flow)

Compute linear flow weight.
"""
flow_weight(::LinearFlowWeighting, normalized_flow) = normalized_flow

"""
    flow_weight(::NoFlowWeighting, normalized_flow)

No flow weighting - returns 1.
"""
flow_weight(::NoFlowWeighting, normalized_flow) = 1.0

"""
    sky_emissivity(air_temperature, vapor_pressure, cloud_fraction)

Calculate effective sky emissivity following Klok & Oerlemans (2002).

Empirical formula - coefficients assume specific units.

# Arguments
- `air_temperature`: Air temperature (K)
- `vapor_pressure`: Vapor pressure (kPa)
- `cloud_fraction`: Cloud cover fraction (0-1)

# Returns
Effective sky emissivity (dimensionless, 0-1)
"""
function sky_emissivity(air_temperature, vapor_pressure, cloud_fraction)
    # Empirical formula - strip to expected units
    T = ustrip(u"K", air_temperature)
    ea = ustrip(u"kPa", vapor_pressure)

    # Clear sky emissivity (Klok & Oerlemans 2002)
    ε_clear = 0.23 + 0.433 * (ea / T)^(1/8)

    # Combined emissivity with cloud effect
    n = cloud_fraction
    ε = ε_clear * (1 - n^2) + 0.976 * n^2
    return ε
end

"""
    mixing_ratio(vapor_pressure, atmospheric_pressure)

Compute water vapor mixing ratio (dimensionless, kg water/kg dry air).

    rv = ε × ea / (P - ea)

where ε ≈ 0.622 is the ratio of molecular weights (water/air).

# Arguments
- `vapor_pressure`: Water vapor pressure (any pressure unit)
- `atmospheric_pressure`: Total atmospheric pressure (any pressure unit)

# Returns
Dimensionless mixing ratio (kg/kg)
"""
function mixing_ratio(vapor_pressure, atmospheric_pressure)
    # Convert to consistent units (Pa) for subtraction
    ea = ustrip(u"Pa", vapor_pressure)
    P = ustrip(u"Pa", atmospheric_pressure)
    return 0.622 * ea / (P - ea)  # dimensionless
end

# -----------------------------------------------------------------------------
# Basin Utilities
# -----------------------------------------------------------------------------

"""
    normalize_flow_by_basin(flow_accumulation, basin_labels)

Normalize flow accumulation as proportion of basin maximum.

Returns flow values scaled 0-1 within each basin.
"""
function normalize_flow_by_basin(flow_accumulation::AbstractMatrix, basin_labels::AbstractMatrix)
    nx, ny = size(flow_accumulation)
    normalized = similar(flow_accumulation, Float64)

    # Find max flow per basin
    basin_max = Dict{eltype(basin_labels), eltype(flow_accumulation)}()
    @inbounds for j in 1:ny, i in 1:nx
        label = basin_labels[i,j]
        flow = flow_accumulation[i,j]
        if !haskey(basin_max, label) || flow > basin_max[label]
            basin_max[label] = flow
        end
    end

    # Normalize
    @inbounds for j in 1:ny, i in 1:nx
        label = basin_labels[i,j]
        max_flow = basin_max[label]
        normalized[i,j] = max_flow > 0 ? flow_accumulation[i,j] / max_flow : 0.0
    end

    return normalized
end

"""
    basin_max_elevation(elevation, basin_labels)

Compute maximum elevation for each basin.

Returns a grid with each cell containing its basin's maximum elevation.
"""
function basin_max_elevation(elevation::AbstractMatrix, basin_labels::AbstractMatrix)
    nx, ny = size(elevation)
    result = similar(elevation)

    # Find max elevation per basin
    basin_max = Dict{eltype(basin_labels), eltype(elevation)}()
    @inbounds for j in 1:ny, i in 1:nx
        label = basin_labels[i,j]
        elev = elevation[i,j]
        if !haskey(basin_max, label) || elev > basin_max[label]
            basin_max[label] = elev
        end
    end

    # Fill result
    @inbounds for j in 1:ny, i in 1:nx
        result[i,j] = basin_max[basin_labels[i,j]]
    end

    return result
end

# =============================================================================
# MacleanColdAirDrainage Model (Static/Diagnostic)
# =============================================================================

"""
    MacleanColdAirDrainage <: ColdAirModel

Static cold air drainage model from Maclean et al. (2019).

A diagnostic approach that computes temperature disturbance from:
- Elevation difference below basin maximum (treated as cold air depth)
- Temperature profile (lapse rate based)
- Flow accumulation weighting

Equation: temperature_change = profile(elevation_difference) × flow_weight

The model is composable - each component can be swapped:
- `conditions`: When drainage occurs (DrainageConditions)
- `temperature_profile`: How temperature varies with depth (TemperatureProfile)
- `flow_weighting`: How flow affects temperature (FlowWeighting)

# Default configuration reproduces Maclean et al. (2019) exactly.

# Reference
Maclean, Mosedale & Bennie (2019) "Microclima: An R package for modelling
meso- and microclimate", Methods Ecol Evol, 10, 280-290.
"""
Base.@kwdef struct MacleanColdAirDrainage{C<:DrainageConditions,P<:TemperatureProfile,F<:FlowWeighting} <: ColdAirModel
    conditions::C = KlokOerlemansConditions()
    temperature_profile::P = LinearProfile(HessLapseRate())
    flow_weighting::F = LogFlowWeighting()
end

"""
    MacleanColdAirState

State for MacleanColdAirDrainage model.

# Fields
- `temperature_disturbance`: Temperature difference due to cold air (K), positive = colder than ambient
- `drainage_active`: Whether drainage conditions are met (Bool)
- `favorable_hours`: Consecutive hours of favorable conditions (Int)
"""
struct MacleanColdAirState{T,A,H}
    temperature_disturbance::T  # K, 2D
    drainage_active::A          # Bool, 2D
    favorable_hours::H          # Int, 2D
end

"""
    MacleanColdAirState(nx, ny)

Create zero-initialized MacleanColdAirState.
"""
function MacleanColdAirState(nx::Int, ny::Int)
    MacleanColdAirState(
        fill(0.0u"K", nx, ny),
        fill(false, nx, ny),
        zeros(Int, nx, ny),
    )
end

"""
    cold_air_step!(model::MacleanColdAirDrainage, state::MacleanColdAirState,
                   env::SpatialEnvironment; kwargs...)

Perform one timestep of the Maclean cold air drainage model.

# Arguments
- `model`: MacleanColdAirDrainage parameters
- `state`: MacleanColdAirState to update
- `env`: SpatialEnvironment containing terrain, weather, land_surface

# Keywords
- `hours_since_sunset`: Hours since sunset (negative = before sunset)
- `hours_since_sunrise`: Hours since sunrise (negative = before sunrise)
- `basin_labels`: Basin ID for each cell (from watershed delineation)
- `flow_accumulation`: Raw flow accumulation (cell counts)
"""
function cold_air_step!(
    model::MacleanColdAirDrainage,
    state::MacleanColdAirState,
    env::SpatialEnvironment;
    hours_since_sunset,
    hours_since_sunrise,
    basin_labels::AbstractMatrix,
    flow_accumulation::AbstractMatrix,
)
    elevation = env.terrain.dem
    nx, ny = size(elevation)

    # Precompute basin properties
    basin_max = basin_max_elevation(elevation, basin_labels)
    normalized_flow = normalize_flow_by_basin(flow_accumulation, basin_labels)

    # Get atmospheric pressure for mixing ratio
    P = env.weather.pressure

    # Get consecutive hours requirement
    consec_required = if model.conditions isa KlokOerlemansConditions
        model.conditions.consecutive_hours_required
    else
        1  # No consecutive requirement for simple conditions
    end

    @inbounds for j in 1:ny
        for i in 1:nx
            # Get cell values from environment
            ws = cell_value(env.weather.wind_speed, i, j)
            T = cell_value(env.weather.air_temperature, i, j)
            cf = cell_value(env.weather.cloud_cover, i, j)
            P_cell = cell_value(P, i, j)

            # Compute vapor pressure from humidity
            ea = vapor_pressure(env, i, j)

            # Compute mixing ratio for HessLapseRate
            rv = mixing_ratio(ea, P_cell)

            # Compute emissivity for this cell
            ε = sky_emissivity(T, ea, cf)

            # Check drainage conditions
            favorable = is_drainage_favorable(
                model.conditions;
                wind_speed = ws,
                emissivity = ε,
                hours_since_sunset = hours_since_sunset,
                hours_since_sunrise = hours_since_sunrise,
            )

            if favorable
                state.favorable_hours[i,j] += 1

                if state.favorable_hours[i,j] >= consec_required
                    state.drainage_active[i,j] = true

                    # Elevation difference is treated as "depth" into cold air pool
                    depth = basin_max[i,j] - elevation[i,j]
                    F = normalized_flow[i,j]

                    if depth > 0.0u"m" && F > 0
                        # Compute base temperature disturbance from profile
                        ΔT_base = temperature_disturbance(
                            model, depth;
                            temperature = T, mixing_ratio = rv
                        )
                        # Apply flow weighting
                        fw = flow_weight(model.flow_weighting, F)
                        state.temperature_disturbance[i,j] = ΔT_base * fw
                    else
                        state.temperature_disturbance[i,j] = 0.0u"K"
                    end
                end
            else
                # Reset if conditions not met
                state.favorable_hours[i,j] = 0
                state.drainage_active[i,j] = false
                state.temperature_disturbance[i,j] = 0.0u"K"
            end
        end
    end

    return state
end


# =============================================================================
# Cold Air Coupling to Surface Temperature
# =============================================================================

"""
    ColdAirCouplingMethod

Abstract type for cold air coupling methods.
"""
abstract type ColdAirCouplingMethod end

"""
    KLAM21Coupling <: ColdAirCouplingMethod

Coupling method that uses KLAM_21 temperature profile to modify surface temperature.

The temperature disturbance from the cold air layer is applied to the
near-surface air temperature.
"""
struct KLAM21Coupling <: ColdAirCouplingMethod end

"""
    apply_cold_air_coupling!(state, method::KLAM21Coupling)

Apply KLAM_21 cold air coupling to surface temperature.

Modifies near-surface temperature based on cold_air.temperature_disturbance from state.
Uses the nested state.cold_air (KLAM21State).
"""
function apply_cold_air_coupling!(
    state,  # SpatialMicroState with nested cold_air::KLAM21State
    ::KLAM21Coupling,
)
    nx, ny, _ = size(state)
    cold_air = state.cold_air

    @inbounds for j in 1:ny
        for i in 1:nx
            # temperature_disturbance is positive = colder than ambient
            if cold_air.depth[i,j] > 0.0u"m"
                state.soil_temperature[i,j,1] -= cold_air.temperature_disturbance[i,j]
            end
        end
    end

    return state
end
