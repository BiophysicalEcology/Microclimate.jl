"""
Snow accumulation and melt models for microclimate simulation.

Three implementations with different complexity/data requirement tradeoffs:

1. `DegreeDaySnow` - Simple temperature-index model requiring only air temperature
2. `Snow17` - NWS operational model with seasonal melt factors and rain-on-snow
3. `UtahEnergyBalance` - Full energy balance model requiring radiation, wind, humidity

References:
- DegreeDaySnow: Hock (2003) "Temperature index melt modelling in mountain areas"
- Snow17: Anderson (1973, 2006) "National Weather Service River Forecast System - Snow Accumulation and Ablation Model"
- UtahEnergyBalance: Tarboton & Luce (1996) "Utah Energy Balance Snow Accumulation and Melt Model"
"""

using Unitful: m, mm, kg, K, °C, W, J, s, hr, d, Pa, ustrip, unit, @u_str

# ============================================================================
# Physical constants
# ============================================================================

# Thermodynamic constants
const STEFAN_BOLTZMANN = 5.67e-8u"W/m^2/K^4"
const LATENT_HEAT_FUSION = 3.335e5u"J/kg"
const LATENT_HEAT_SUBLIMATION = 2.838e6u"J/kg"
const SPECIFIC_HEAT_AIR = 1005.0u"J/kg/K"
const SPECIFIC_HEAT_ICE = 2100.0u"J/kg/K"
const GAS_CONSTANT_WATER_VAPOR = 461.0u"J/kg/K"
const VON_KARMAN = 0.4

# Material densities
const DENSITY_WATER = 1000.0u"kg/m^3"
const DENSITY_ICE = 917.0u"kg/m^3"
const DENSITY_AIR = 1.225u"kg/m^3"

# Default snow properties
const FRESH_SNOW_DENSITY_DEFAULT = 100.0u"kg/m^3"
const FRESH_SNOW_ALBEDO_DEFAULT = 0.85

# ============================================================================
# Helper functions
# ============================================================================

"""
    snow_depth(swe, density)

Compute snow depth from snow water equivalent and density.
"""
snow_depth(swe, density) = swe * DENSITY_WATER / density

# ============================================================================
# Environment forcing struct
# ============================================================================

"""
    SnowForcing

Container for atmospheric forcing data required by snow models.

All snow models require at minimum `air_temperature` and `precipitation`.
Energy balance models additionally require radiation, wind, and humidity fields.

# Fields
- `air_temperature`: Air temperature (°C or K)
- `precipitation`: Total precipitation (m or mm water equivalent)
- `shortwave_radiation`: Incoming shortwave radiation (W/m²), default 0
- `longwave_radiation`: Incoming longwave radiation (W/m²), default 0
- `wind_speed`: Wind speed at reference height (m/s), default 0
- `relative_humidity`: Relative humidity (0-1), default 0.5
- `atmospheric_pressure`: Atmospheric pressure (Pa), default 101325
- `ground_heat_flux`: Ground heat flux (W/m²), default 0 (Tarboton & Luce 1996)
- `net_energy_flux`: Net energy flux at snow surface (W/m²), for simplified models
- `day_of_year`: Day of year (1-366), for seasonal melt factor calculation
"""
Base.@kwdef struct SnowForcing{T,P,SW,LW,WS,RH,AP,GHF,NE}
    air_temperature::T
    precipitation::P
    shortwave_radiation::SW = 0.0u"W/m^2"
    longwave_radiation::LW = 0.0u"W/m^2"
    wind_speed::WS = 0.0u"m/s"
    relative_humidity::RH = 0.5
    atmospheric_pressure::AP = 101325.0u"Pa"
    ground_heat_flux::GHF = 0.0u"W/m^2"
    net_energy_flux::NE = 0.0u"W/m^2"
    day_of_year::Int = 1
end

# ============================================================================
# Abstract type and interface
# ============================================================================

abstract type SnowModel end

# ============================================================================
# Snow accumulation methods
# ============================================================================

"""
    SnowAccumulation

Abstract type for precipitation partitioning methods.

Subtypes determine how precipitation is split into rain and snow based on temperature.
"""
abstract type SnowAccumulation end

"""
    NoAccumulation()

All precipitation passes through as rain (no snow accumulation).
"""
struct NoAccumulation <: SnowAccumulation end

"""
    ThresholdAccumulation(; threshold=1.0u"°C")

Hard temperature threshold for rain/snow partitioning.

Below `threshold`: all snow. At or above `threshold`: all rain.
"""
Base.@kwdef struct ThresholdAccumulation{T} <: SnowAccumulation
    threshold::T = 1.0u"°C"
end

"""
    LinearTransitionAccumulation(; snow_threshold=0.0u"°C", rain_threshold=2.0u"°C")

Linear transition between rain and snow over a temperature range.

- At or below `snow_threshold`: all snow
- At or above `rain_threshold`: all rain
- Between: linear interpolation
"""
Base.@kwdef struct LinearTransitionAccumulation{TS,TR} <: SnowAccumulation
    snow_threshold::TS = 0.0u"°C"
    rain_threshold::TR = 2.0u"°C"
end

# Dispatch snow_accumulation on accumulation method
function snow_accumulation(::NoAccumulation, precipitation, air_temperature)
    return (snowfall=zero(precipitation), rainfall=precipitation)
end

function snow_accumulation(model::ThresholdAccumulation, precipitation, air_temperature)
    if air_temperature <= model.threshold
        return (snowfall=precipitation, rainfall=zero(precipitation))
    else
        return (snowfall=zero(precipitation), rainfall=precipitation)
    end
end

function snow_accumulation(model::LinearTransitionAccumulation, precipitation, air_temperature)
    T = air_temperature
    T_snow = model.snow_threshold
    T_rain = model.rain_threshold

    if T <= T_snow
        snow_fraction = 1.0
    elseif T >= T_rain
        snow_fraction = 0.0
    else
        snow_fraction = (T_rain - T) / (T_rain - T_snow)
    end

    snowfall = precipitation * snow_fraction
    rainfall = precipitation * (1.0 - snow_fraction)
    return (snowfall=snowfall, rainfall=rainfall)
end

# ============================================================================
# Snow model types
# ============================================================================

"""
    NoSnow()

Disable snow modeling. This is the default for `MicroProblem` and `SpatialMicroProblem`.
"""
struct NoSnow <: SnowModel end

# NoSnow uses NoAccumulation
snow_accumulation(::NoSnow, precipitation, air_temperature) =
    snow_accumulation(NoAccumulation(), precipitation, air_temperature)

# update_snow_state for NoSnow is defined after SnowState

"""
    snow_accumulation(model, precipitation, air_temperature)

Partition precipitation into rain and snow based on air temperature.

Returns `(snowfall, rainfall)` in the same units as precipitation.
"""
function snow_accumulation end

# Generic method for all SnowModels with an accumulation field
snow_accumulation(model::SnowModel, precipitation, air_temperature) =
    snow_accumulation(model.accumulation, precipitation, air_temperature)

"""
    snow_melt(model, state::SnowState, forcing, timestep)

Compute snowmelt given current snow state and atmospheric forcing.

Returns melt amount (length units, same as SWE).

# Arguments
- `model`: Snow model
- `state`: Current snow state (`SnowState`)
- `forcing`: Atmospheric forcing (`SnowForcing` or compatible object)
- `timestep`: Time step duration
"""
function snow_melt end

# ============================================================================
# Snow state container (defined early for use in function signatures)
# ============================================================================

"""
    SnowState

Container for snow state variables at a single point or grid cell.

# Fields
- `water_equivalent`: Snow water equivalent (m)
- `depth`: Physical snow depth (m)
- `density`: Snow density (kg/m³)
- `temperature`: Snow temperature (°C)
- `albedo`: Snow surface albedo (0-1)
- `age`: Days since last snowfall
- `surface_age`: Dimensionless snow surface age τ for Dickinson albedo (grain growth based)
- `liquid_water`: Liquid water content (m water equivalent)
- `cold_content`: Energy deficit before melt can occur (mm water equivalent)
- `energy_content`: Energy content relative to ice at 0°C (J/m²), used by UEB
- `antecedent_temperature_index`: Weighted average of recent temperatures (°C), used by Snow17
- `lagged_excess_storage`: Storage of lagged water in transit through pack (mm), used by Snow17
"""
Base.@kwdef struct SnowState{SWE,D,RHO,T,A,AGE,SA,LW,CC,EC,ATI,LES}
    water_equivalent::SWE = 0.0u"m"
    depth::D = 0.0u"m"
    density::RHO = FRESH_SNOW_DENSITY_DEFAULT
    temperature::T = 0.0u"°C"
    albedo::A = FRESH_SNOW_ALBEDO_DEFAULT
    age::AGE = 0.0u"d"
    surface_age::SA = 0.0  # dimensionless τ for Dickinson albedo
    liquid_water::LW = 0.0u"m"
    cold_content::CC = 0.0u"mm"
    energy_content::EC = 0.0u"J/m^2"
    antecedent_temperature_index::ATI = 0.0u"°C"
    lagged_excess_storage::LES = 0.0u"mm"
end

# ============================================================================
# Thermal conductivity formulas
# ============================================================================

abstract type ThermalConductivityFormula end

"""
    Djachkova()

Djachkova's exponential formula for snow thermal conductivity.

    k = 0.0442 × exp(5.181 × ρ)

where ρ is density in g/cm³.

# References
- Djachkova, as cited in Anderson (2006) and Kearney (2020)
"""
struct Djachkova <: ThermalConductivityFormula end

"""
    snow_thermal_conductivity(::Djachkova, density)

Compute snow thermal conductivity using Djachkova's formula.
"""
function snow_thermal_conductivity(::Djachkova, density)
    ρ = ustrip(u"g/cm^3", density)
    return 0.0442 * exp(5.181 * ρ) * u"W/m/K"
end

# ============================================================================
# Density evolution formulas
# ============================================================================

"""
    SnowDensityFormula

Abstract type for snow density evolution methods.

# Available Formulas

| Formula | Physics | Performance | Use Case |
|---------|---------|-------------|----------|
| `SimpleMixingDensity` | Mass-weighted mixing only | ~3 ns | Fast screening, when depth not critical |
| `AndersonSnowDensity` | Temperature-dependent fresh snow + mixing | ~5 ns | Better fresh snow density |
| `AndersonDensityEvolution` | Full compaction + metamorphism (Eq 14-17) | ~35 ns | Accurate depth evolution |
| `CompactionSnowDensity` | Exponential relaxation to maximum | ~10 ns | Simple time-dependent compaction |
| `SturmSnowDensity` | Depth + age dependent (empirical) | ~15 ns | Empirical depth-density relationship |

# Formula Selection Guide

**`SimpleMixingDensity`** (default for Snow17, DegreeDaySnow):
- No time-dependent evolution - density only changes with new snowfall
- Fastest option, suitable when depth accuracy is not critical
- Good for: rapid screening, SWE-focused studies

**`AndersonSnowDensity`**:
- Temperature-dependent fresh snow density (50 kg/m³ at -15°C, ~150 kg/m³ at 0°C)
- No time-dependent compaction
- Good for: improved fresh snow representation without computational overhead

**`AndersonDensityEvolution`** (SNOW-17 Eq 14-17):
- Full physics: compaction (overburden), destructive metamorphism, liquid water acceleration
- Density increases over time even without snowfall
- Deep snow compacts faster; wet snow compacts faster; warm snow compacts faster
- Good for: accurate snow depth evolution, insulation studies, long-duration snowpacks

**`CompactionSnowDensity`** (default for UtahEnergyBalance):
- Exponential relaxation toward maximum density
- Simpler than Anderson, but includes time dependence
- Good for: moderate complexity needs

**`SturmSnowDensity`** (default for KearneySnow):
- Empirical formula based on depth and age
- Calibrated to field observations
- Good for: empirically-grounded estimates

# Interface

All formulas implement:
```julia
snow_density(formula, state, snowfall, timestep)
snow_density(formula, state, snowfall, timestep, air_temperature)  # for temperature-dependent
```

# Example

```julia
# Default: simple mixing
model = Snow17()

# Full density evolution
model = Snow17(density_formula=AndersonDensityEvolution())
```
"""
abstract type SnowDensityFormula end

"""
    SturmSnowDensity(; maximum, fresh, depth_coefficient, age_coefficient)

Sturm et al. (2010) density evolution formula.

    ρ = ρ_max - (ρ_max - ρ_fresh) × exp(-k₁h - k₂D)

where h is depth (cm), D is days since snowfall.

# References
- Sturm, M. et al. (2010). Estimating snow water equivalent using snow depth data and climate classes.
"""
Base.@kwdef struct SturmSnowDensity{M,F,DC,AC} <: SnowDensityFormula
    maximum::M = 597.9u"kg/m^3"
    fresh::F = 217.8u"kg/m^3"
    depth_coefficient::DC = 0.001u"cm^-1"
    age_coefficient::AC = 0.0038u"d^-1"
end

"""
    snow_density(formula, state, snowfall, timestep)
    snow_density(formula, state, snowfall, timestep, air_temperature)

Compute/evolve snow density. All formulas use the same interface.

- `state`: SnowState with depth, age, density, water_equivalent
- `snowfall`: new snowfall this timestep
- `timestep`: time step duration
- `air_temperature`: (optional) air temperature for temperature-dependent formulas
"""
function snow_density end

# Default 5-argument interface: delegate to 4-argument (ignore temperature)
function snow_density(formula::SnowDensityFormula, state, snowfall, timestep, air_temperature)
    return snow_density(formula, state, snowfall, timestep)
end

function snow_density(formula::SturmSnowDensity, state, snowfall, timestep)
    h = ustrip(u"cm", state.depth)
    D = ustrip(u"d", state.age)
    k1 = ustrip(u"cm^-1", formula.depth_coefficient)
    k2 = ustrip(u"d^-1", formula.age_coefficient)
    ρ_max = formula.maximum
    ρ_fresh = formula.fresh
    return ρ_max - (ρ_max - ρ_fresh) * exp(-k1 * h - k2 * D)
end

"""
    CompactionSnowDensity(; fresh=100.0u"kg/m^3", maximum=500.0u"kg/m^3", timescale=21.0u"d")

Exponential compaction density evolution formula.

    ρ_new = ρ_max - (ρ_max - ρ) × exp(-Δt / τ)

Fresh snow is mixed by mass-weighted averaging.

# References
- Tarboton, D.G. and C.H. Luce (1996). Utah Energy Balance Snow Accumulation and Melt Model.
"""
Base.@kwdef struct CompactionSnowDensity{F,M,T} <: SnowDensityFormula
    fresh::F = 100.0u"kg/m^3"
    maximum::M = 500.0u"kg/m^3"
    timescale::T = 21.0u"d"
end

function snow_density(formula::CompactionSnowDensity, state, snowfall, timestep)
    ρ = state.density
    swe = state.water_equivalent
    ρ_fresh = formula.fresh
    ρ_max = formula.maximum
    τ = formula.timescale

    if swe <= zero(swe)
        return ρ_fresh
    end

    # Mix with fresh snow
    if snowfall > zero(snowfall)
        swe_old = swe - snowfall
        if swe_old > zero(swe_old)
            ρ = (swe_old * ρ + snowfall * ρ_fresh) / swe
        else
            ρ = ρ_fresh
        end
    end

    # Compaction over time toward maximum density
    return ρ_max - (ρ_max - ρ) * exp(-timestep / τ)
end

"""
    SimpleMixingDensity(; fresh=100.0u"kg/m^3")

Simple mass-weighted density mixing with no compaction.

When snowfall occurs, new density is the mass-weighted average of
existing snow and fresh snow densities.

Used by simple temperature-index models (DegreeDaySnow, Snow17).
"""
Base.@kwdef struct SimpleMixingDensity{F} <: SnowDensityFormula
    fresh::F = FRESH_SNOW_DENSITY_DEFAULT
end

function snow_density(formula::SimpleMixingDensity, state, snowfall, timestep)
    swe = state.water_equivalent
    ρ_fresh = formula.fresh

    if swe <= zero(swe)
        return ρ_fresh
    end

    if snowfall > zero(snowfall)
        swe_old = swe - snowfall
        if swe_old > zero(swe_old)
            return (swe_old * state.density + snowfall * ρ_fresh) / swe
        else
            return ρ_fresh
        end
    end

    return state.density
end

"""
    AndersonSnowDensity(; minimum_density=50.0, coefficient=1.7)

Anderson (1976) temperature-dependent new snow density formula.
Used by SNOW-17 model.

    ρn = 50 kg/m³  when Ta ≤ -15°C
    ρn = 50 + 1.7 × (Ta + 15)^1.5  when Ta > -15°C

Results in ~50 kg/m³ at -15°C, increasing to ~150 kg/m³ at 0°C.

When snowfall occurs, new density is the mass-weighted average of
existing snow and fresh snow densities (same as SimpleMixingDensity).

# References
- Anderson, E.A. (1976). A Point Energy and Mass Balance Model of a Snow Cover.
  NOAA Technical Report NWS 19.
- Anderson, E.A. (2006). Snow Accumulation and Ablation Model – SNOW-17.
  Equations 2a-2b.
"""
Base.@kwdef struct AndersonSnowDensity <: SnowDensityFormula
    minimum_density::Float64 = 50.0   # kg/m³ at Ta ≤ -15°C
    coefficient::Float64 = 1.7        # for (Ta+15)^1.5 term in kg/m³
end

function snow_density(formula::AndersonSnowDensity, state, snowfall, timestep, air_temperature)
    Ta_C = ustrip(u"°C", air_temperature)
    swe = state.water_equivalent

    # Compute fresh snow density (Anderson Eq 2a-2b)
    if Ta_C <= -15.0
        ρ_fresh = formula.minimum_density * u"kg/m^3"
    else
        Ta_shifted = Ta_C + 15.0  # Shift so -15°C → 0
        ρ_fresh = (formula.minimum_density + formula.coefficient * Ta_shifted^1.5) * u"kg/m^3"
    end

    # No existing snow
    if swe <= zero(swe)
        return ρ_fresh
    end

    # Mass-weighted mixing with existing snow
    if snowfall > zero(snowfall)
        swe_old = swe - snowfall
        if swe_old > zero(swe_old)
            return (swe_old * state.density + snowfall * ρ_fresh) / swe
        else
            return ρ_fresh
        end
    end

    return state.density
end

"""
    AndersonDensityEvolution(; kwargs...)

Anderson (1976/2006) time-dependent snow density evolution with compaction
and destructive metamorphism.

Implements Equation 14 from SNOW-17 documentation:
- Compaction due to overburden pressure (B term)
- Destructive metamorphism converting fresh flakes to rounded grains (A term)
- Enhanced settling rate when liquid water is present

The density evolves even without new snowfall, producing realistic
densification of aging snowpacks.

# Parameters
- `c1`: Compaction rate coefficient (0.026 cm⁻¹·hr⁻¹)
- `c2`: Kojima viscosity coefficient (21 cm³·g⁻¹)
- `c3`: Destructive metamorphism rate at 0°C (0.005 hr⁻¹)
- `c4`: Temperature coefficient for metamorphism (0.10 °C⁻¹)
- `c5_wet`: Settling rate multiplier when liquid water present (2.0)
- `cx`: Decay factor for metamorphism above threshold density (23)
- `ρd`: Threshold density for reduced metamorphism (150 kg/m³)
- `ρmax`: Maximum allowed density (600 kg/m³)
- `minimum_density`: Minimum fresh snow density at cold temps (50 kg/m³)
- `density_coefficient`: Fresh snow density temperature coefficient (1.7)

# References
- Anderson, E.A. (1976). A Point Energy and Mass Balance Model of a Snow Cover.
- Anderson, E.A. (2006). Snow Accumulation and Ablation Model – SNOW-17, Eq 14-17.
"""
Base.@kwdef struct AndersonDensityEvolution <: SnowDensityFormula
    # Compaction constants
    c1::Float64 = 0.026      # cm⁻¹·hr⁻¹ - fractional increase due to overburden
    c2::Float64 = 21.0       # cm³·g⁻¹ - Kojima viscosity coefficient

    # Destructive metamorphism constants
    c3::Float64 = 0.005      # hr⁻¹ - settling rate at 0°C for ρ < ρd
    c4::Float64 = 0.10       # °C⁻¹ - temperature coefficient
    c5_wet::Float64 = 2.0    # multiplier when liquid water present
    cx::Float64 = 23.0       # decay factor when ρ > ρd

    # Density thresholds
    ρd::Float64 = 150.0      # kg/m³ - threshold for reduced metamorphism
    ρmax::Float64 = 600.0    # kg/m³ - maximum allowed density

    # Fresh snow density (Anderson Eq 2a-2b)
    minimum_density::Float64 = 50.0   # kg/m³ at Ta ≤ -15°C
    density_coefficient::Float64 = 1.7
end

function snow_density(formula::AndersonDensityEvolution, state, snowfall, timestep, air_temperature)
    Ta_C = ustrip(u"°C", air_temperature)
    dt_hr = ustrip(u"hr", timestep)
    swe = state.water_equivalent
    ρ_old = state.density

    # Fresh snow density (Eq 2a-2b)
    if Ta_C <= -15.0
        ρ_fresh = formula.minimum_density
    else
        Ta_shifted = Ta_C + 15.0
        ρ_fresh = formula.minimum_density + formula.density_coefficient * Ta_shifted^1.5
    end
    ρ_fresh = ρ_fresh * u"kg/m^3"

    # No existing snow - return fresh snow density
    if swe <= zero(swe)
        return ρ_fresh
    end

    # Convert to spec units for Eq 14
    ρx = ustrip(u"g/cm^3", ρ_old)  # density in g/cm³
    Wix = ustrip(u"mm", swe)       # SWE in mm

    # Snow temperature - use state temperature if available, else estimate from air
    Ts = ustrip(u"°C", state.temperature)
    Ts = min(Ts, 0.0)  # Snow can't be above 0°C

    # Check for liquid water presence (c5 factor)
    liquid_water = ustrip(u"mm", state.liquid_water)
    c5 = liquid_water > 0.0 ? formula.c5_wet : 0.0

    # Density threshold in g/cm³
    ρd = formula.ρd / 1000.0  # convert kg/m³ to g/cm³

    # β factor: 0 if ρx ≤ ρd, 1 if ρx > ρd
    β = ρx > ρd ? 1.0 : 0.0

    # Equation 14: Compaction and metamorphism
    # B = c1 × Δt × exp(0.08×Ts - c2×ρx)
    B = formula.c1 * dt_hr * exp(0.08 * Ts - formula.c2 * ρx)

    # A = c3 × c5 × Δt × exp(c4×Ts - cx×β×(ρx - ρd))
    A = formula.c3 * (1.0 + c5) * dt_hr * exp(formula.c4 * Ts - formula.cx * β * (ρx - ρd))

    # ρx2 = ρx1 × (exp(B×0.1×Wix) - 1) / (B×0.1×Wix) × exp(A)
    # Protect against division by zero
    B_term = B * 0.1 * Wix
    if abs(B_term) < 1e-10
        compaction_factor = 1.0
    else
        compaction_factor = (exp(B_term) - 1.0) / B_term
    end

    ρx_new = ρx * compaction_factor * exp(A)

    # Apply maximum density constraint
    ρmax_gcm3 = formula.ρmax / 1000.0
    ρx_new = min(ρx_new, ρmax_gcm3)

    # Convert back to kg/m³
    ρ_evolved = ρx_new * 1000.0 * u"kg/m^3"

    # Mix with fresh snow if there's new snowfall
    if snowfall > zero(snowfall)
        swe_old = swe - snowfall
        if swe_old > zero(swe_old)
            # Mass-weighted average of evolved old snow and fresh snow
            return (swe_old * ρ_evolved + snowfall * ρ_fresh) / swe
        else
            return ρ_fresh
        end
    end

    return ρ_evolved
end

# ============================================================================
# Albedo decay formulas
# ============================================================================

abstract type AlbedoFormula end

"""
    fresh_snow_albedo(formula)

Return fresh snow albedo for the given formula. Used when resetting state.
"""
fresh_snow_albedo(formula::AlbedoFormula) = formula.fresh

"""
    NoAlbedo(; fresh=0.85)

Constant albedo with no decay. Used by simple temperature-index models.
"""
Base.@kwdef struct NoAlbedo{F} <: AlbedoFormula
    fresh::F = FRESH_SNOW_ALBEDO_DEFAULT
end

snow_albedo(formula::NoAlbedo, days_since_snowfall) = formula.fresh

"""
    AndersonAlbedo(; fresh=0.85, minimum=0.45)

Anderson (2006) logarithmic albedo decay formula from Figure A-4.

    α = -9.874 × ln(days) + 78.34 (as percentage)

Clamped to range [minimum, fresh].

# References
- Anderson, E.A. (2006). Snow Accumulation and Ablation Model – SNOW-17, Figure A-4.
"""
Base.@kwdef struct AndersonAlbedo{F,M} <: AlbedoFormula
    fresh::F = FRESH_SNOW_ALBEDO_DEFAULT
    minimum::M = 0.45
end

"""
    snow_albedo(formula::AndersonAlbedo, days_since_snowfall)

Compute snow albedo using Anderson decay formula.
"""
function snow_albedo(formula::AndersonAlbedo, days_since_snowfall)
    days = ustrip(u"d", days_since_snowfall)
    if days < 1.0
        return formula.fresh
    end
    α_percent = -9.874 * log(days) + 78.34
    return clamp(α_percent / 100.0, formula.minimum, formula.fresh)
end

"""
    ExponentialSnowAlbedo(; fresh=0.85, minimum=0.5, decay_rate=0.1u"d^-1")

Exponential albedo decay formula.

    α(t) = α_min + (α_fresh - α_min) × exp(-k × t)

where t is days since snowfall, k is decay rate.

# References
- Tarboton, D.G. and C.H. Luce (1996). Utah Energy Balance Snow Accumulation and Melt Model.
"""
Base.@kwdef struct ExponentialSnowAlbedo{F,M,K} <: AlbedoFormula
    fresh::F = 0.85
    minimum::M = 0.5
    decay_rate::K = 0.1u"d^-1"
end

"""
    snow_albedo(formula::ExponentialSnowAlbedo, days_since_snowfall)

Compute snow albedo using exponential decay formula.
"""
function snow_albedo(formula::ExponentialSnowAlbedo, days_since_snowfall)
    α_fresh = formula.fresh
    α_min = formula.minimum
    k = formula.decay_rate
    return α_min + (α_fresh - α_min) * exp(-k * days_since_snowfall)
end

"""
    DickinsonAlbedo(; αvo=0.85, αiro=0.65, Cv=0.2, Cir=0.5, r3=0.03, τo=1e6, b=2.0,
                     bare_ground_albedo=0.25, shallow_threshold=0.1u"m")

Dickinson et al. (1993) two-band albedo formula with grain growth aging.

Computes albedo as average of visible and near-infrared bands, with aging
driven by grain growth (temperature-dependent), melt/refreeze, and dirt/soot.
Includes illumination angle adjustment and shallow snow interpolation.

# Parameters
- `αvo`: Fresh snow visible band reflectance (0.85)
- `αiro`: Fresh snow near-infrared band reflectance (0.65)
- `Cv`: Visible band sensitivity to aging (0.2)
- `Cir`: Near-infrared band sensitivity to aging (0.5)
- `r3`: Dirt/soot aging factor (0.03, use 0.01 for Antarctica)
- `τo`: Reference timescale for aging (10⁶ s)
- `b`: Illumination angle parameter (2.0)
- `bare_ground_albedo`: Albedo of bare ground for shallow snow interpolation (0.25)
- `shallow_threshold`: Depth below which shallow snow interpolation applies (0.1 m)

# State
Uses `state.surface_age` (dimensionless τ) which must be updated each timestep.

# References
- Dickinson, R.E. et al. (1993). Biosphere-Atmosphere Transfer Scheme (BATS)
  version 1e as coupled to the NCAR Community Climate Model, p21.
- Tarboton, D.G. and C.H. Luce (1996). Utah Energy Balance Snow Accumulation
  and Melt Model (UEB), Equations 17-25.
"""
Base.@kwdef struct DickinsonAlbedo{AVO,AIRO,CV,CIR,R3,TO,B,BG,ST} <: AlbedoFormula
    αvo::AVO = 0.85           # fresh snow visible reflectance
    αiro::AIRO = 0.65         # fresh snow near-infrared reflectance
    Cv::CV = 0.2              # visible band aging sensitivity
    Cir::CIR = 0.5            # near-infrared band aging sensitivity
    r3::R3 = 0.03             # dirt/soot factor (0.01 for Antarctica)
    τo::TO = 1.0e6            # reference timescale (seconds)
    b::B = 2.0                # illumination angle parameter
    bare_ground_albedo::BG = 0.25
    shallow_threshold::ST = 0.1u"m"
end

"""
    update_surface_age(formula::DickinsonAlbedo, τ, snow_temperature, timestep, snowfall)

Update dimensionless snow surface age τ based on grain growth.

Implements Equations 19-22 from Tarboton & Luce (1996):
- r1 = exp[5000 × (1/273.16 - 1/Ts)] (grain growth from vapor diffusion)
- r2 = min(r1^10, 1) (melt/refreeze effect)
- r3 = 0.03 (dirt/soot)
- Δτ = (r1 + r2 + r3) × Δt / τo

Snowfall ≥ 0.01 m resets τ to 0. Smaller snowfall reduces τ by factor (1 - 100×Ps).

Returns new surface age τ.
"""
function update_surface_age(formula::DickinsonAlbedo, τ, snow_temperature, timestep, snowfall)
    # Snowfall reset (spec: 0.01 m restores to new)
    snowfall_m = ustrip(u"m", snowfall)
    if snowfall_m >= 0.01
        return 0.0
    elseif snowfall_m > 0.0
        # Partial snowfall reduces age
        τ = τ * (1.0 - 100.0 * snowfall_m)
    end

    # Snow surface temperature in Kelvin
    Ts_K = ustrip(u"K", uconvert(u"K", snow_temperature))
    Ts_K = min(Ts_K, 273.16)  # Can't exceed freezing

    # r1: grain growth from vapor diffusion (Eq 21)
    r1 = exp(5000.0 * (1.0/273.16 - 1.0/Ts_K))

    # r2: melt/refreeze effect (Eq 22)
    r2 = min(r1^10, 1.0)

    # r3: dirt/soot (parameter)
    r3 = formula.r3

    # Age increment (Eq 20)
    Δt_s = ustrip(u"s", timestep)
    Δτ = (r1 + r2 + r3) * Δt_s / formula.τo

    return τ + Δτ
end

"""
    snow_albedo(formula::DickinsonAlbedo, surface_age; depth=nothing, cos_illumination_angle=nothing)

Compute snow albedo using Dickinson two-band formula.

# Arguments
- `formula`: DickinsonAlbedo parameters
- `surface_age`: Dimensionless snow surface age τ
- `depth`: Snow depth for shallow snow interpolation (optional)
- `cos_illumination_angle`: Cosine of solar zenith angle (optional)

# Returns
Snow albedo (0-1)
"""
function snow_albedo(formula::DickinsonAlbedo, surface_age; depth=nothing, cos_illumination_angle=nothing)
    τ = surface_age

    # Age function (Eq 19)
    Fage = τ / (1.0 + τ)

    # Diffuse reflectances (Eq 17-18)
    αvd = (1.0 - formula.Cv * Fage) * formula.αvo
    αird = (1.0 - formula.Cir * Fage) * formula.αiro

    # Illumination angle adjustment (Eq 23-25)
    if !isnothing(cos_illumination_angle)
        cosψ = cos_illumination_angle
        if cosψ < 0.5
            b = formula.b
            fψ = (1.0/b) * ((b + 1.0)/(1.0 + 2.0*b*cosψ) - 1.0)
        else
            fψ = 0.0
        end
        αv = αvd + 0.4 * fψ * (1.0 - αvd)
        αir = αird + 0.4 * fψ * (1.0 - αird)
    else
        αv = αvd
        αir = αird
    end

    # Average of two bands
    α_snow = (αv + αir) / 2.0

    # Shallow snow interpolation with bare ground
    if !isnothing(depth)
        z = ustrip(u"m", depth)
        h = ustrip(u"m", formula.shallow_threshold)
        if z < h && z > 0.0
            # r = (1 - z/h) × exp(-z/2h)
            r = (1.0 - z/h) * exp(-z/(2.0*h))
            α_snow = r * formula.bare_ground_albedo + (1.0 - r) * α_snow
        elseif z <= 0.0
            α_snow = formula.bare_ground_albedo
        end
    end

    return α_snow
end

# Fresh snow albedo for DickinsonAlbedo (average of visible and infrared bands)
fresh_snow_albedo(formula::DickinsonAlbedo) = (formula.αvo + formula.αiro) / 2.0

# Fallback for DickinsonAlbedo with days_since_snowfall (for compatibility)
# This uses a rough conversion but update_surface_age should be used properly
function snow_albedo(formula::DickinsonAlbedo, days_since_snowfall::Unitful.Time)
    # Rough approximation: convert days to dimensionless age assuming moderate temperature
    # This is for compatibility only; proper usage tracks surface_age in state
    days = ustrip(u"d", days_since_snowfall)
    # Approximate τ from days (assuming r1≈1, r2≈1, r3=0.03 at near-freezing temps)
    τ_approx = days * 86400.0 * 2.03 / formula.τo
    return snow_albedo(formula, τ_approx)
end

# Generic fallback: update_surface_age does nothing for non-Dickinson formulas
update_surface_age(::AlbedoFormula, surface_age, snow_temperature, timestep, snowfall) = surface_age

"""
    compute_albedo(formula, state, snowfall, timestep)

Compute snow albedo, updating surface_age state if needed.
Returns (albedo, new_surface_age).
"""
function compute_albedo(formula::AlbedoFormula, state, snowfall, timestep)
    # Default: use days-based interface, surface_age unchanged
    return snow_albedo(formula, state.age), state.surface_age
end

function compute_albedo(formula::DickinsonAlbedo, state, snowfall, timestep)
    # Update dimensionless surface age based on grain growth
    new_surface_age = update_surface_age(formula, state.surface_age,
                                          state.temperature, timestep, snowfall)
    # Compute albedo with surface_age and depth
    albedo = snow_albedo(formula, new_surface_age; depth=state.depth)
    return albedo, new_surface_age
end

# ============================================================================
# Melt outflow formulas
# ============================================================================

abstract type MeltOutflowFormula end

"""
    MeltOutflowForcing

Container for all inputs needed by melt outflow formulas.

Constructed by the calling code from SnowState, SnowForcing, and pre-computed values.
Each formula extracts what it needs.

# Fields
- `air_temperature`: Air temperature
- `snow_temperature`: Snow surface/pack temperature
- `water_equivalent`: Snow water equivalent (m)
- `density`: Snow density (kg/m³)
- `energy_content`: Energy content relative to 0°C ice (J/m²), for UEB
- `cold_content`: Cold content deficit (mm), for Snow17
- `liquid_water`: Liquid water in pack (m), for Snow17
- `net_energy_flux`: Net energy flux at surface (W/m²), for Kearney
- `rainfall`: Rainfall this timestep (m), after rain/snow partitioning
- `day_of_year`: Day of year (1-366), for seasonal melt factor
- `timestep`: Time step duration
"""
Base.@kwdef struct MeltOutflowForcing{TA,TS,SWE,RHO,EC,CC,LW,NEF,RF,DOY,DT}
    air_temperature::TA
    snow_temperature::TS
    water_equivalent::SWE
    density::RHO
    energy_content::EC
    cold_content::CC
    liquid_water::LW
    net_energy_flux::NEF
    rainfall::RF
    day_of_year::DOY
    timestep::DT
end

"""
    DarcyMeltOutflow(; saturated_hydraulic_conductivity=20.0u"m/hr", capillary_retention=0.05)

Darcy's law melt outflow formula from Tarboton & Luce (1996), Equations 47-48.

Melt outflow rate is proportional to the cube of relative saturation:

    Mr = Ksat × S*³

where S* is the relative saturation in excess of capillary retention:

    S* = (Lf × ρs/ρi - Lc) / ((1-Lf) × (ρw/ρi - ρw/ρs) - Lc)

and Lf = U/(ρw × hf × W) is the liquid fraction of the snowpack.

# Parameters
- `saturated_hydraulic_conductivity`: Ksat (m/hr), calibrated value 20 m/hr (Table 1)
- `capillary_retention`: Lc, liquid retained by capillary forces as fraction of solid (0.05)

# References
- Tarboton, D.G. and C.H. Luce (1996). Utah Energy Balance Snow Accumulation
  and Melt Model (UEB), Equations 47-48.
- Male and Gray (1981), p. 400, Equation 9.45.
"""
Base.@kwdef struct DarcyMeltOutflow{K,L} <: MeltOutflowFormula
    saturated_hydraulic_conductivity::K = 20.0u"m/hr"  # Ksat
    capillary_retention::L = 0.05                       # Lc
end

"""
    melt_outflow(formula::DarcyMeltOutflow, mof::MeltOutflowForcing)

Compute melt outflow using Darcy's law (Equations 47-48).

Uses `mof.energy_content`, `mof.water_equivalent`, `mof.density`, `mof.timestep`.
"""
function melt_outflow(formula::DarcyMeltOutflow, mof::MeltOutflowForcing)
    U = mof.energy_content
    W = mof.water_equivalent
    ρs = mof.density
    timestep = mof.timestep

    zero_melt = zero(W)

    # No melt if energy content is negative (cold snow) or no snow
    if U <= zero(U) || W <= zero(W)
        return zero_melt
    end

    # Energy to melt all snow
    energy_to_melt_all = W * DENSITY_WATER * LATENT_HEAT_FUSION

    # Liquid fraction Lf = U / (ρw × hf × W) (from Eq 4 context)
    Lf = ustrip(U / energy_to_melt_all)
    Lf = clamp(Lf, 0.0, 1.0)

    # Relative saturation S* (Equation 48)
    # S* = (Lf × ρs/ρi - Lc) / ((1-Lf) × (ρw/ρi - ρw/ρs) - Lc)
    ρs_val = ustrip(u"kg/m^3", ρs)
    ρi = ustrip(u"kg/m^3", DENSITY_ICE)
    ρw = ustrip(u"kg/m^3", DENSITY_WATER)
    Lc = formula.capillary_retention

    # Numerator: liquid water volume fraction minus capillary retention
    numerator = Lf * ρs_val / ρi - Lc

    # Denominator: pore volume fraction minus capillary retention
    denominator = (1.0 - Lf) * (ρw / ρi - ρw / ρs_val) - Lc

    if numerator <= 0.0 || denominator <= 0.0
        # No excess liquid above capillary retention
        return zero_melt
    end

    S_star = clamp(numerator / denominator, 0.0, 1.0)

    # Melt outflow rate Mr = Ksat × S*³ (Equation 47)
    Ksat = formula.saturated_hydraulic_conductivity
    Mr = Ksat * S_star^3

    # Melt over timestep, limited by available liquid water
    # Available liquid = Lf × W (liquid fraction of total SWE)
    available_liquid = Lf * W
    return min(Mr * timestep, available_liquid)
end

"""
    EnergyLimitedMeltOutflow()

Simple energy-limited melt outflow (non-UEB).

All energy above zero is immediately converted to melt outflow.
Used by simpler models that don't track liquid water routing.
"""
struct EnergyLimitedMeltOutflow <: MeltOutflowFormula end

"""
    melt_outflow(formula::EnergyLimitedMeltOutflow, mof::MeltOutflowForcing)

Compute melt outflow by direct energy-to-melt conversion.

Uses `mof.energy_content`, `mof.water_equivalent`.
"""
function melt_outflow(::EnergyLimitedMeltOutflow, mof::MeltOutflowForcing)
    U = mof.energy_content
    W = mof.water_equivalent

    if U <= zero(U) || W <= zero(W)
        return zero(W)
    end

    # Direct conversion of positive energy to melt
    melt = U / (LATENT_HEAT_FUSION * DENSITY_WATER)
    return min(melt, W)
end

"""
    Snow17MeltOutflow(; maximum_melt_factor, minimum_melt_factor, melt_base_temperature,
                       rain_threshold, rain_melt_coefficient, wind_function, elevation,
                       daily_ground_melt, latitude)

Snow17 melt outflow formula (Anderson 2006).

Implements Equations 5, 6, 7, and 26 from the SNOW-17 documentation:
- Rain-on-snow melt (Eq 5) when rain intensity > threshold
- Non-rain melt (Eq 6) using degree-day approach with seasonal melt factor
- Seasonal melt factor with Av adjustment (Eq 7)
- Ground melt (Eq 26)

Handles cold content deficit before melt can occur.

# Parameters
- `maximum_melt_factor`: MFMAX - maximum melt factor at summer solstice (mm/K/6hr)
- `minimum_melt_factor`: MFMIN - minimum melt factor at winter solstice (mm/K/6hr)
- `melt_base_temperature`: MBASE - temperature above which melt occurs (0°C)
- `rain_threshold`: Rain intensity threshold for rain-on-snow mode (0.25 mm/hr)
- `rain_melt_coefficient`: Rain-induced melt coefficient (0.0125)
- `wind_function`: UADJ - wind function for turbulent transfer (0.04 mm/mbar/6hr)
- `elevation`: Site elevation for atmospheric pressure calculation (0 m)
- `daily_ground_melt`: Constant ground melt rate (0.0 mm/day)
- `latitude`: Latitude for seasonal melt factor calculation (degrees)

# References
- Anderson, E.A. (2006). Snow Accumulation and Ablation Model – SNOW-17.
"""
Base.@kwdef struct Snow17MeltOutflow{MFX,MFN,MBT,RT,RMC,WF,EL,DGM,LAT} <: MeltOutflowFormula
    maximum_melt_factor::MFX = 1.05u"mm/K/(6hr)"
    minimum_melt_factor::MFN = 0.6u"mm/K/(6hr)"
    melt_base_temperature::MBT = 0.0u"°C"
    rain_threshold::RT = 0.25u"mm/hr"
    rain_melt_coefficient::RMC = 0.0125
    wind_function::WF = 0.04u"mm/mbar/(6hr)"
    elevation::EL = 0.0u"m"
    daily_ground_melt::DGM = 0.0u"mm/d"
    latitude::LAT = 45.0
end

"""
    seasonal_melt_factor(formula::Snow17MeltOutflow, day_of_year)

Compute seasonal melt factor for Snow17 (Anderson Eq 7).
"""
function seasonal_melt_factor(formula::Snow17MeltOutflow, day_of_year::Integer)
    mfmax = formula.maximum_melt_factor
    mfmin = formula.minimum_melt_factor
    lat = formula.latitude

    # Seasonal variation using sine curve
    if lat >= 0  # Northern hemisphere
        phase = 2π * (day_of_year - 81) / 365
    else  # Southern hemisphere
        phase = 2π * (day_of_year - 81 + 182) / 365
    end

    # Compute Av factor for high latitudes (Anderson Eq 7)
    Av = compute_av_factor(lat, day_of_year)

    # Mf = (mfmax + mfmin)/2 + Av × (mfmax - mfmin)/2 × sin(phase)
    mf = (mfmax + mfmin) / 2 + Av * (mfmax - mfmin) / 2 * sin(phase)
    return mf
end

"""
    melt_outflow(formula::Snow17MeltOutflow, mof::MeltOutflowForcing)

Compute melt outflow using Snow17 approach (Anderson 2006, Eq 5, 6, 26).

Uses `mof.air_temperature`, `mof.rainfall`, `mof.day_of_year`, `mof.cold_content`,
`mof.water_equivalent`, `mof.liquid_water`, `mof.timestep`.
"""
function melt_outflow(formula::Snow17MeltOutflow, mof::MeltOutflowForcing)
    swe = mof.water_equivalent

    if swe <= zero(swe)
        return zero(swe)
    end

    # Compute seasonal melt factor from day of year
    melt_factor = seasonal_melt_factor(formula, mof.day_of_year)

    T_melt = formula.melt_base_temperature
    timestep = mof.timestep
    rainfall = mof.rainfall
    cold_content = mof.cold_content
    liquid_water = mof.liquid_water

    # Scale melt factor to timestep (base is 6-hourly)
    timestep_hr = ustrip(u"hr", timestep)
    mf_scaled = melt_factor * (timestep_hr / 6.0)

    # Get temperature values
    Ta_C = ustrip(u"°C", mof.air_temperature)
    T_melt_C = ustrip(u"°C", T_melt)

    # Rain intensity in mm/hr
    rain_intensity = ustrip(u"mm/hr", rainfall / timestep)
    rain_threshold = ustrip(u"mm/hr", formula.rain_threshold)

    # Determine if rain-on-snow conditions (Eq 5) or non-rain melt (Eq 6)
    if rain_intensity > rain_threshold
        # Rain-on-snow melt (Anderson 2006, Equation 5)
        σ = 6.12e-10
        Ta_K = Ta_C + 273.0
        longwave_melt = σ * timestep_hr * (Ta_K^4 - 273.4^4)

        Tr_C = max(Ta_C, 0.0)
        rainfall_mm = ustrip(u"mm", rainfall)
        rain_heat_melt = formula.rain_melt_coefficient * rainfall_mm * Tr_C

        UADJ = ustrip(u"mm/mbar/(6hr)", formula.wind_function)
        esat = snow17_saturation_vapor_pressure(Ta_C)
        Pa = snow17_atmospheric_pressure(ustrip(u"m", formula.elevation))

        vapor_term = 0.9 * esat - 6.11
        sensible_term = 0.00057 * Pa * Ta_C

        turbulent_melt = 8.5 * UADJ * (timestep_hr / 6.0) * (vapor_term + sensible_term)

        melt = (longwave_melt + rain_heat_melt + turbulent_melt) * u"mm"
    else
        # Non-rain melt (Anderson 2006, Equation 6)
        if Ta_C > T_melt_C
            ΔT = Ta_C - T_melt_C
            potential_melt = ustrip(u"mm", mf_scaled) * ΔT

            if rainfall > zero(rainfall)
                rainfall_mm = ustrip(u"mm", rainfall)
                Tr_C = max(Ta_C, 0.0)
                potential_melt += formula.rain_melt_coefficient * rainfall_mm * Tr_C
            end

            cc_mm = ustrip(u"mm", cold_content)
            if cc_mm > 0.0
                if potential_melt >= cc_mm
                    potential_melt -= cc_mm
                else
                    potential_melt = 0.0
                end
            end

            melt = potential_melt * u"mm"
        else
            melt = zero(swe)
        end
    end

    # Add ground melt (Equation 26)
    ground_melt = formula.daily_ground_melt * ustrip(u"d", timestep)
    melt += ground_melt

    # Can't melt more than available
    max_melt = swe + liquid_water
    melt = min(melt, max_melt)

    return melt
end

"""
    EnergyFluxMeltOutflow(; rain_melt_factor=0.0125, temperature_threshold=-0.5,
                           minimum_snow_depth=0.025u"m")

Energy flux driven melt outflow formula matching Kearney/NicheMapR approach.

Melt is computed from:
1. Net energy flux when snow temperature ≥ threshold
2. Rain-on-snow enhancement (Anderson 2006)

# Parameters
- `rain_melt_factor`: Melt enhancement per unit rainfall per degree (0.0125)
- `temperature_threshold`: Snow temperature above which energy melt occurs (-0.5°C)
- `minimum_snow_depth`: Below this depth, melt remaining snow completely (0.025 m)

# References
- Kearney, M. (2020). NicheMapR microclimate model.
- Anderson, E.A. (2006). Snow Accumulation and Ablation Model – SNOW-17.
"""
Base.@kwdef struct EnergyFluxMeltOutflow{RMF,TT,MSD} <: MeltOutflowFormula
    rain_melt_factor::RMF = 0.0125
    temperature_threshold::TT = -0.5  # °C
    minimum_snow_depth::MSD = 0.025u"m"
end

"""
    melt_outflow(formula::EnergyFluxMeltOutflow, mof::MeltOutflowForcing)

Compute melt outflow using energy flux approach (Kearney/NicheMapR).

Uses `mof.net_energy_flux`, `mof.snow_temperature`, `mof.air_temperature`,
`mof.rainfall`, `mof.water_equivalent`, `mof.density`, `mof.timestep`.
"""
function melt_outflow(formula::EnergyFluxMeltOutflow, mof::MeltOutflowForcing)
    swe = mof.water_equivalent
    snow_temperature = mof.snow_temperature
    air_temperature = mof.air_temperature
    rainfall = mof.rainfall
    density = mof.density
    timestep = mof.timestep
    net_energy_flux = mof.net_energy_flux

    if swe <= zero(swe)
        return zero(swe)
    end

    melt = zero(swe)

    # Get temperature in °C for comparisons
    T_snow_C = ustrip(u"°C", snow_temperature)

    # Energy-driven melt
    if net_energy_flux > 0.0u"W/m^2" && T_snow_C >= formula.temperature_threshold
        # Energy available for melting
        energy_available = net_energy_flux * timestep

        # Mass that can be melted
        melt_mass = energy_available / LATENT_HEAT_FUSION
        melt = melt_mass / DENSITY_WATER  # convert to m water equivalent
    end

    # Rain-on-snow melt (Anderson 2006)
    T_air_C = ustrip(u"°C", air_temperature)
    if rainfall > zero(rainfall) && T_air_C > 0
        rain_melt = rainfall * T_air_C * formula.rain_melt_factor
        melt += rain_melt
    end

    # Can't melt more than available
    melt = min(melt, swe)

    # Check minimum depth threshold
    remaining = swe - melt
    min_swe = formula.minimum_snow_depth * density / DENSITY_WATER
    if remaining < min_swe
        # Melt remaining snow completely
        melt = swe
    end

    return melt
end

# ============================================================================
# Temperature evolution formulas
# ============================================================================

abstract type SnowTemperatureFormula end

"""
    snow_temperature(formula, state, forcing, melt, timestep)

Compute new snow temperature. All formulas use the same interface.

- `state`: Current SnowState
- `forcing`: SnowForcing with air_temperature, net_energy_flux, etc.
- `melt`: Melt amount this timestep
- `timestep`: Time step duration
"""
function snow_temperature end

"""
    NoTemperatureEvolution()

No temperature tracking - returns current temperature unchanged.
Used by simple degree-day models.
"""
struct NoTemperatureEvolution <: SnowTemperatureFormula end

function snow_temperature(::NoTemperatureEvolution, state, forcing, melt, timestep)
    return state.temperature
end

"""
    RelaxationTemperature(; relaxation_rate=0.1, use_energy_flux=false)

Relaxation toward target temperature.

- If warming (melt > 0, or net_energy_flux > 0 if use_energy_flux=true): target is 0°C
- Otherwise: target is min(0°C, air_temperature)

# Parameters
- `relaxation_rate`: Rate of relaxation toward target (0-1)
- `use_energy_flux`: If true, use net_energy_flux > 0 as warming condition (Kearney style).
                     If false, use melt > 0 as warming condition (Utah style).

# References
- Tarboton, D.G. and C.H. Luce (1996). Utah Energy Balance Snow Accumulation and Melt Model.
- Kearney, M.R. (2020). How will snow alter exposure of organisms to cold stress.
"""
Base.@kwdef struct RelaxationTemperature{R,E} <: SnowTemperatureFormula
    relaxation_rate::R = 0.1
    use_energy_flux::E = false
end

function snow_temperature(formula::RelaxationTemperature, state, forcing, melt, timestep)
    if state.water_equivalent <= zero(state.water_equivalent)
        return state.temperature
    end

    T_snow_C = ustrip(u"°C", state.temperature)
    T_air_C = ustrip(u"°C", forcing.air_temperature)
    r = formula.relaxation_rate

    # Determine if warming (relaxing toward 0°C)
    if formula.use_energy_flux
        warming = forcing.net_energy_flux > 0.0u"W/m^2"
    else
        warming = melt > zero(melt)
    end

    if warming
        T_target = 0.0
    else
        T_target = min(0.0, T_air_C)
    end

    return (T_snow_C + r * (T_target - T_snow_C)) * u"°C"
end

"""
    ColdContentTemperature()

Compute temperature from cold content (Snow17 approach).

    T = -cold_content × 160 / snow_water_equivalent

# References
- Anderson, E.A. (2006). Snow Accumulation and Ablation Model – SNOW-17.
"""
struct ColdContentTemperature <: SnowTemperatureFormula end

function snow_temperature(::ColdContentTemperature, state, forcing, melt, timestep)
    if state.water_equivalent <= zero(state.water_equivalent)
        return 0.0u"°C"
    end

    if state.cold_content > 0.0u"mm"
        swe_mm = ustrip(u"mm", state.water_equivalent)
        cc_mm = ustrip(u"mm", state.cold_content)
        T_est = -cc_mm * 160.0 / max(swe_mm, 1.0)
        return max(-40.0, T_est) * u"°C"
    else
        return 0.0u"°C"
    end
end

"""
    EnergyContentTemperature(; soil_depth=0.4u"m", soil_density=1500.0u"kg/m^3", soil_heat_capacity=840.0u"J/kg/K")

Derive temperature from energy content (Utah Energy Balance approach).

Temperature is diagnosed from energy content U and water equivalent W:
- U < 0: T = U / (ρw W Cs + ρg De Cg)  (all solid)
- 0 ≤ U ≤ ρw W hf: T = 0°C  (solid-liquid mixture)
- U > ρw W hf: T = (U - ρw W hf) / (ρg De Cg + ρw W Cw)  (all liquid, W=0)

# Parameters
- `soil_depth`: Depth of soil interacting thermally with snowpack (De)
- `soil_density`: Soil density (ρg)
- `soil_heat_capacity`: Soil specific heat (Cg)

# References
- Tarboton, D.G. and C.H. Luce (1996). Utah Energy Balance Snow Accumulation and Melt Model.
  Equations 3-5.
"""
Base.@kwdef struct EnergyContentTemperature{SD,SRHO,SHC} <: SnowTemperatureFormula
    soil_depth::SD = 0.4u"m"
    soil_density::SRHO = 1500.0u"kg/m^3"
    soil_heat_capacity::SHC = 840.0u"J/kg/K"
end

function snow_temperature(formula::EnergyContentTemperature, state, forcing, melt, timestep)
    U = state.energy_content
    W = state.water_equivalent

    if W <= zero(W)
        return 0.0u"°C"
    end

    # Heat capacities
    De = formula.soil_depth
    ρg = formula.soil_density
    Cg = formula.soil_heat_capacity
    Cs = SPECIFIC_HEAT_ICE
    Cw = 4180.0u"J/kg/K"  # specific heat of water

    # Heat capacity of snow (per unit area)
    C_snow = DENSITY_WATER * W * Cs

    # Heat capacity of soil layer (per unit area)
    C_soil = ρg * De * Cg

    # Energy to melt all snow
    E_melt = DENSITY_WATER * W * LATENT_HEAT_FUSION

    if U < zero(U)
        # All solid phase (equation 3)
        # T = U / (ρw W Cs + ρg De Cg) gives temperature relative to 0°C
        T = U / (C_snow + C_soil)
        return max(-40.0u"°C", uconvert(u"°C", T))
    elseif U <= E_melt
        # Solid-liquid mixture at 0°C (equation 4)
        return 0.0u"°C"
    else
        # All liquid - practically W=0 here (equation 5)
        return 0.0u"°C"
    end
end

# ============================================================================
# Thermal properties
# ============================================================================

"""
    snow_specific_heat(density, temperature)

Compute snow specific heat capacity, including latent heat near freezing.

# Arguments
- `density`: Snow density (kg/m³)
- `temperature`: Snow temperature (°C or K)

# Returns
Specific heat capacity (J/kg/K)
"""
function snow_specific_heat(density, temperature)
    # Weighted average of ice and air specific heats
    ice_fraction = density / DENSITY_ICE
    c_snow = SPECIFIC_HEAT_ICE * ice_fraction + SPECIFIC_HEAT_AIR * (1 - ice_fraction)

    # Add latent heat of fusion near freezing point
    # This represents the phase change energy sink/source
    T_celsius = ustrip(u"°C", temperature)
    if -0.5 < T_celsius < 0.5
        c_snow += LATENT_HEAT_FUSION / u"K"
    end

    return c_snow
end

# ============================================================================
# Melt models
# ============================================================================

abstract type MeltModel end

"""
    DegreeDayMelt(; melt_factor=3.0u"mm/K/d", melt_threshold=0.0u"°C")

Simple degree-day (temperature-index) melt model.

Melt is proportional to positive degree-days above a threshold temperature.
Suitable for data-sparse regions or rapid screening studies.

# Parameters
- `melt_factor`: Melt rate per degree-day (typical range 1.6-6.0 mm/K/d)
- `melt_threshold`: Temperature above which melt occurs

# References
- Hock, R. (2003). Temperature index melt modelling in mountain areas.
  Journal of Hydrology, 282(1-4), 104-115.
"""
Base.@kwdef struct DegreeDayMelt{MF,TT} <: MeltModel
    melt_factor::MF = 3.0u"mm/K/d"  # K is used because melt is proportional to temperature *difference*
    melt_threshold::TT = 0.0u"°C"
end

function snow_melt(model::DegreeDayMelt, state::SnowState, forcing, timestep)
    (; air_temperature) = forcing
    swe = state.water_equivalent

    # Hock (2003) Equation 2: M = 0 if Td ≤ T0
    if air_temperature <= model.melt_threshold || swe <= zero(swe)
        return zero(swe)
    end

    # Hock (2003) Equation 2: M = fm × (Td - T0)
    # Scaled by timestep to generalize beyond daily intervals
    degree_days = (air_temperature - model.melt_threshold) * timestep
    melt = model.melt_factor * degree_days

    # Can't melt more than available
    return min(melt, swe)
end

# ============================================================================
# GenericSnowModel - Composable snow model
# ============================================================================

"""
    GenericSnowModel(; melt, accumulation, albedo_formula, density_formula, temperature_formula)

Composable snow model that combines independent melt, accumulation, and evolution formulas.

This allows mixing and matching components from different sources. For example,
use `DegreeDayMelt` (Hock 2003) for melt with `SturmSnowDensity` (Sturm 2010)
for density evolution.

# Parameters
- `melt`: Melt model (e.g., `DegreeDayMelt()`)
- `accumulation`: Snow accumulation method (default: `ThresholdAccumulation()`)
- `albedo_formula`: Albedo decay formula (default: `NoAlbedo()`)
- `density_formula`: Density evolution formula (default: `SimpleMixingDensity()`)
- `temperature_formula`: Temperature evolution formula (default: `NoTemperatureEvolution()`)

# Example
```julia
model = GenericSnowModel(
    melt = DegreeDayMelt(melt_factor=4.0u"mm/K/d"),
    accumulation = LinearTransitionAccumulation(),
    density_formula = SturmSnowDensity(),
)
```
"""
Base.@kwdef struct GenericSnowModel{M<:MeltModel,A<:SnowAccumulation,AF<:AlbedoFormula,DF<:SnowDensityFormula,TF<:SnowTemperatureFormula} <: SnowModel
    melt::M = DegreeDayMelt()
    accumulation::A = ThresholdAccumulation()
    albedo_formula::AF = NoAlbedo()
    density_formula::DF = SimpleMixingDensity()
    temperature_formula::TF = NoTemperatureEvolution()
end

# Delegate melt calculation to the melt model
function snow_melt(model::GenericSnowModel, state::SnowState, forcing, timestep)
    return snow_melt(model.melt, state, forcing, timestep)
end

# ============================================================================
# DegreeDaySnow - Convenience alias for GenericSnowModel with DegreeDayMelt
# ============================================================================

"""
    DegreeDaySnow(; melt_factor=3.0u"mm/K/d", melt_threshold=0.0u"°C", ...)

Simple degree-day snow model. Convenience constructor for `GenericSnowModel`
with `DegreeDayMelt`.

# Parameters
- `melt_factor`: Melt rate per degree-day (typical range 1.6-6.0 mm/K/d)
- `melt_threshold`: Temperature above which melt occurs
- `accumulation`: Snow accumulation method (default: `ThresholdAccumulation()`)
- `albedo_formula`: Albedo decay formula (default: `NoAlbedo()`)
- `density_formula`: Density evolution formula (default: `SimpleMixingDensity()`)
- `temperature_formula`: Temperature evolution formula (default: `NoTemperatureEvolution()`)

# References
- Hock, R. (2003). Temperature index melt modelling in mountain areas.
  Journal of Hydrology, 282(1-4), 104-115.
"""
function DegreeDaySnow(;
    melt_factor = 3.0u"mm/K/d",
    melt_threshold = 0.0u"°C",
    accumulation = ThresholdAccumulation(),
    albedo_formula = NoAlbedo(),
    density_formula = SimpleMixingDensity(),
    temperature_formula = NoTemperatureEvolution(),
)
    melt = DegreeDayMelt(; melt_factor, melt_threshold)
    return GenericSnowModel(; melt, accumulation, albedo_formula, density_formula, temperature_formula)
end

# ============================================================================
# Snow17 - NWS Temperature Index Model
# ============================================================================

"""
    Snow17(; maximum_melt_factor=1.05u"mm/K/(6hr)", minimum_melt_factor=0.6u"mm/K/(6hr)", ...)

NWS SNOW-17 temperature index snow accumulation and ablation model.

Uses temperature as primary driver with seasonal variation in melt factors,
separate rain-on-snow and non-rain melt equations, and cold content tracking.
Requires only temperature and precipitation inputs.

# Implementation Status

**Fully implemented:**
- Form of precipitation / rain-snow partitioning (Eq 1)
- New snow density, temperature-dependent with `AndersonSnowDensity` (Eq 2a-2b)
- Cold content from snowfall (Eq 4)
- Rain-on-snow melt with longwave, turbulent transfer, rain heat (Eq 5)
- Non-rain melt with seasonal melt factor (Eq 6-7)
- Alaska/high-latitude Av factor adjustment for ≥54°N/S (Eq 7)
- Antecedent Temperature Index (ATI) for cold content evolution (Eq 8)
- Cold content change from temperature gradient (Eq 9)
- Liquid water holding capacity based on ice portion (Eq 10)
- Water transmission lag and attenuation through snowpack (Eq 21-24)
- Ground melt (Eq 26)
- Snowfall correction factor (SCF) for gage undercatch (Eq 1)

**Available via formula selection:**
- Detailed density evolution with compaction/metamorphism (Eq 14-17):
  use `density_formula=AndersonDensityEvolution()`

**Not implemented:**
- Areal extent of snow cover / depletion curves (Eq 29-30)
- User-specified seasonal melt factor variation curves

# Parameters
- `maximum_melt_factor`: Maximum melt factor at summer solstice (MFMAX, mm/K/6hr)
- `minimum_melt_factor`: Minimum melt factor at winter solstice (MFMIN, mm/K/6hr)
- `negative_melt_factor`: Maximum negative melt factor for cold content (NMF, mm/K/6hr)
- `temperature_index_parameter`: Antecedent temperature index parameter (TIPM, 0-1)
- `liquid_water_holding_capacity`: Fraction of ice held as liquid (PLWHC, 0-0.4)
- `latitude`: Latitude for seasonal melt factor and Av factor calculation (degrees)
- `snowfall_correction_factor`: Gage undercatch multiplier (SCF, typically 1.0-1.5)
- `melt_outflow_formula`: Melt calculation formula (default: `Snow17MeltOutflow()`)
- `density_formula`: Density evolution formula (default: `SimpleMixingDensity()`)
  - `AndersonSnowDensity()`: Temperature-dependent fresh snow density (Eq 2a-2b)
  - `AndersonDensityEvolution()`: Full compaction + metamorphism (Eq 14-17)
- `accumulation`: Snow accumulation method (default: `LinearTransitionAccumulation`)

# References
- Anderson, E.A. (1973). National Weather Service River Forecast System -
  Snow Accumulation and Ablation Model. NOAA Technical Memorandum NWS HYDRO-17.
- Anderson, E.A. (2006). Snow Accumulation and Ablation Model – SNOW-17.
  NWS Technical Documentation.
"""
Base.@kwdef struct Snow17{MFX,MFN,NMF,TIPM,LWC,LAT,SCF,MF<:MeltOutflowFormula,AF<:AlbedoFormula,DF<:SnowDensityFormula,TF<:SnowTemperatureFormula,A<:SnowAccumulation} <: SnowModel
    maximum_melt_factor::MFX = 1.05u"mm/K/(6hr)"  # MFMAX
    minimum_melt_factor::MFN = 0.6u"mm/K/(6hr)"   # MFMIN
    negative_melt_factor::NMF = 0.15u"mm/K/(6hr)" # NMF - maximum negative melt factor
    temperature_index_parameter::TIPM = 0.1       # TIPM - ATI weighting (0-1)
    liquid_water_holding_capacity::LWC = 0.04     # PLWHC - fraction of ice portion
    latitude::LAT = 45.0                          # degrees
    snowfall_correction_factor::SCF = 1.0         # SCF - gage undercatch multiplier (Eq 1)
    melt_outflow_formula::MF = Snow17MeltOutflow()
    albedo_formula::AF = NoAlbedo()
    density_formula::DF = SimpleMixingDensity()
    temperature_formula::TF = ColdContentTemperature()
    accumulation::A = LinearTransitionAccumulation(snow_threshold=0.0u"°C", rain_threshold=2.0u"°C")
end

"""
    compute_av_factor(latitude, day_of_year)

Compute the Av factor for high-latitude seasonal melt factor adjustment.
Anderson (2006) Equation 7.

For latitudes < 54°N/S, returns 1.0 (full seasonal variation).
For latitudes ≥ 54°N/S (Alaska/high-latitude):
- Av = 0.0 from Sep 24 to Mar 18 (no seasonal variation in winter)
- Av = 1.0 from Apr 27 to Aug 15 (full seasonal variation in summer)
- Linear interpolation in transition periods
"""
function compute_av_factor(latitude, day_of_year)
    if abs(latitude) < 54.0
        return 1.0
    end

    # Day-of-year thresholds (approximate, non-leap year)
    # Spring transition: Mar 18 (day 77) to Apr 27 (day 117) - 40 days
    # Summer full variation: Apr 27 (day 117) to Aug 15 (day 227)
    # Fall transition: Aug 15 (day 227) to Sep 24 (day 267) - 40 days
    # Winter no variation: Sep 24 (day 267) to Mar 18 (day 77 next year)

    doy = day_of_year

    # Adjust for southern hemisphere (6 month offset)
    if latitude < 0
        doy = mod(doy + 182, 365)
        if doy == 0
            doy = 365
        end
    end

    if doy >= 117 && doy <= 227
        # Summer: full variation
        return 1.0
    elseif doy >= 267 || doy <= 77
        # Winter: no variation (Av = 0 means melt factor stays at average)
        return 0.0
    elseif doy > 77 && doy < 117
        # Spring transition: linear interpolation from 0 to 1
        return (doy - 77) / 40.0
    else  # doy > 227 && doy < 267
        # Fall transition: linear interpolation from 1 to 0
        return 1.0 - (doy - 227) / 40.0
    end
end

"""
    seasonal_melt_factor(model::Snow17, day_of_year)

Compute seasonally-varying melt factor based on day of year.
Includes Av factor adjustment for high latitudes (≥54°N/S).
Anderson (2006) Equation 7.
"""
function seasonal_melt_factor(model::Snow17, day_of_year::Integer)
    mfmax = model.maximum_melt_factor
    mfmin = model.minimum_melt_factor
    lat = model.latitude

    # Seasonal variation using sine curve
    # Maximum at summer solstice (day ~172 in NH), minimum at winter solstice
    if lat >= 0  # Northern hemisphere
        phase = 2π * (day_of_year - 81) / 365  # 81 = spring equinox
    else  # Southern hemisphere
        phase = 2π * (day_of_year - 81 + 182) / 365
    end

    # Compute Av factor for high latitudes (Anderson Eq 7)
    Av = compute_av_factor(lat, day_of_year)

    # Mf = (mfmax + mfmin)/2 + Av × (mfmax - mfmin)/2 × sin(phase)
    mf = (mfmax + mfmin) / 2 + Av * (mfmax - mfmin) / 2 * sin(phase)
    return mf
end

"""
    snow17_saturation_vapor_pressure(Ta_C)

Compute saturation vapor pressure at temperature Ta (°C) in millibars.
Anderson (2006) Equation 5 formula.
"""
function snow17_saturation_vapor_pressure(Ta_C)
    return 2.7489e8 * exp(-4278.63 / (Ta_C + 242.792))
end

"""
    snow17_atmospheric_pressure(elevation_m)

Compute atmospheric pressure at given elevation using standard atmosphere.
Anderson (2006) formula below Equation 5. Returns pressure in millibars.
"""
function snow17_atmospheric_pressure(elevation_m)
    He = elevation_m
    return 33.86 * (29.9 - 0.335 * He / 100 + 0.00022 * (He / 100)^2.4)
end

"""
    compute_water_lag(excess_mm, ice_mm, timestep_hr)

Compute lag time for excess water through the snowpack.
Anderson (2006) Equation 21.

Returns lag time in hours (maximum 5.33 hours).
"""
function compute_water_lag(excess_mm, ice_mm, timestep_hr)
    if excess_mm <= 0 || ice_mm <= 0
        return 0.0
    end
    # L = 5.33 × [1 - exp((-0.03 × (Δtp/6) × Wi) / E)]
    return 5.33 * (1.0 - exp((-0.03 * (timestep_hr / 6.0) * ice_mm) / excess_mm))
end

"""
    compute_withdrawal_rate(lagged_water_inches, ice_inches)

Compute one-hour withdrawal rate for water draining from storage.
Anderson (2006) Equation 22.

Parameters are in inches per the original formulation.
Returns withdrawal rate as a fraction per hour.
"""
function compute_withdrawal_rate(lagged_water_inches, ice_inches)
    if ice_inches <= 0
        return 1.0  # Drain immediately if no snow
    end
    if lagged_water_inches <= 0
        return 0.0  # Nothing to withdraw
    end
    # R1 = 1.0 / (1.0 + 5.0 × exp((-500 × Els) / Wis))^1.3
    ratio = (500.0 * lagged_water_inches) / ice_inches
    return 1.0 / (1.0 + 5.0 * exp(-ratio))^1.3
end

"""
    route_lagged_water(excess_water_mm, ice_mm, storage_mm, timestep_hr)

Route excess water through lagged storage with attenuation.
Anderson (2006) Equations 21-24.

Returns (outflow_mm, new_storage_mm).
"""
function route_lagged_water(excess_water_mm, ice_mm, storage_mm, timestep_hr)
    # Handle edge cases
    if ice_mm <= 0
        # No snow - all water flows through immediately
        return excess_water_mm + storage_mm, 0.0
    end

    # Compute lag time for new excess water
    if excess_water_mm > 0
        lag_hours = compute_water_lag(excess_water_mm, ice_mm, timestep_hr)
    else
        lag_hours = 0.0
    end

    # Break excess water into increments based on lag
    # For simplicity, we'll use a single lag value and add water gradually
    # The model operates hourly internally for attenuation

    # Convert to inches for withdrawal rate calculation
    mm_to_inches = 1.0 / 25.4
    ice_inches = ice_mm * mm_to_inches

    # Sub-step hourly for attenuation (Eq 22-24)
    current_storage = storage_mm
    total_outflow = 0.0

    # Number of hours to process (at least 1)
    num_hours = max(1, Int(floor(timestep_hr)))

    # Distribute excess water input over the timestep, accounting for lag
    # Simplified: add excess water at the start of the period
    # In full implementation, would spread based on lag time
    if lag_hours < 1.0 && excess_water_mm > 0
        # Short lag - add most water to storage immediately
        current_storage += excess_water_mm
    elseif excess_water_mm > 0
        # Longer lag - add water gradually over the lag period
        water_per_hour = excess_water_mm / min(lag_hours, Float64(num_hours))
        current_storage += water_per_hour  # First hour's input
    end

    for hour in 1:num_hours
        # Add lagged excess water input for subsequent hours
        if hour > 1 && excess_water_mm > 0 && lag_hours >= 1.0
            water_per_hour = excess_water_mm / min(lag_hours, Float64(num_hours))
            if hour <= ceil(lag_hours)
                current_storage += water_per_hour
            end
        end

        # Compute withdrawal rate (Eq 22)
        lagged_water_inches = current_storage * mm_to_inches
        R1 = compute_withdrawal_rate(lagged_water_inches, ice_inches)

        # Hourly outflow (Eq 23): O_mr1 = S × R1
        hourly_outflow = current_storage * R1

        # Update storage (Eq 24): S2 = S1 - O_mr1
        current_storage = max(0.0, current_storage - hourly_outflow)
        total_outflow += hourly_outflow
    end

    # Handle sub-hourly remainder
    remainder_hours = timestep_hr - num_hours
    if remainder_hours > 0 && current_storage > 0
        lagged_water_inches = current_storage * mm_to_inches
        R1 = compute_withdrawal_rate(lagged_water_inches, ice_inches)
        partial_outflow = current_storage * R1 * remainder_hours
        current_storage = max(0.0, current_storage - partial_outflow)
        total_outflow += partial_outflow
    end

    return total_outflow, current_storage
end

function snow_melt(model::Snow17, state::SnowState, forcing, timestep)
    (; air_temperature, precipitation, day_of_year) = forcing

    # Compute rainfall (precipitation minus snowfall)
    rainfall = precipitation - snow_accumulation(model, precipitation, air_temperature).snowfall

    # Create forcing struct for melt outflow
    mof = MeltOutflowForcing(
        air_temperature=air_temperature,
        snow_temperature=state.temperature,
        water_equivalent=state.water_equivalent,
        density=state.density,
        energy_content=state.energy_content,
        cold_content=state.cold_content,
        liquid_water=state.liquid_water,
        net_energy_flux=0.0u"W/m^2",  # Snow17 doesn't use energy flux
        rainfall=rainfall,
        day_of_year=day_of_year,
        timestep=timestep,
    )
    return melt_outflow(model.melt_outflow_formula, mof)
end

# ============================================================================
# UtahEnergyBalance - Full energy balance model
# ============================================================================

"""
    UtahEnergyBalance(; albedo_formula=DickinsonAlbedo(), density_formula=CompactionSnowDensity(), ...)

Utah Energy Balance (UEB) snow model.

Full energy balance approach using physically-based calculations of radiative,
sensible, latent, and advective heat exchanges. Requires radiation, wind, and
humidity inputs.

# Parameters
- `albedo_formula`: Albedo formula (default: DickinsonAlbedo - two-band with grain growth)
- `density_formula`: Density evolution formula (default: CompactionSnowDensity)
- `melt_outflow_formula`: Melt outflow formula (default: DarcyMeltOutflow - Eq 47-48)
- `snow_roughness_height`: Aerodynamic roughness length zo (m), calibrated value 0.005 m
- `snow_emissivity`: Longwave emissivity of snow surface εs (0.99)
- `accumulation`: Snow accumulation method (default: LinearTransitionAccumulation)

# References
- Tarboton, D.G. and C.H. Luce (1996). Utah Energy Balance Snow Accumulation
  and Melt Model (UEB). Utah Water Research Laboratory and USDA Forest Service
  Intermountain Research Station.
- Dickinson, R.E. et al. (1993). Biosphere-Atmosphere Transfer Scheme (BATS).
- You, J. (2004). Snow Hydrology: The Parameterization of Subgrid Processes
  within a Physically Based Snow Energy and Mass Balance Model. PhD Dissertation.
"""
Base.@kwdef struct UtahEnergyBalance{AF<:AlbedoFormula,DF<:SnowDensityFormula,TF<:SnowTemperatureFormula,MF<:MeltOutflowFormula,SRH,SE,A<:SnowAccumulation} <: SnowModel
    albedo_formula::AF = DickinsonAlbedo()
    density_formula::DF = CompactionSnowDensity()
    temperature_formula::TF = EnergyContentTemperature(soil_density=1700.0u"kg/m^3", soil_heat_capacity=2090.0u"J/kg/K")
    melt_outflow_formula::MF = DarcyMeltOutflow()
    snow_roughness_height::SRH = 0.005u"m"
    snow_emissivity::SE = 0.99
    accumulation::A = LinearTransitionAccumulation(snow_threshold=-1.0u"°C", rain_threshold=3.0u"°C")
end

"""
    snow_energy_balance(model::UtahEnergyBalance, state::SnowState, forcing, timestep)

Compute snow energy balance and resulting melt.

# Arguments
- `model`: UtahEnergyBalance model with parameters
- `state`: Current snow state (SnowState)
- `forcing`: Atmospheric forcing (SnowForcing)
- `timestep`: Time step duration

# Returns
Named tuple with energy fluxes and melt rate.
"""
function snow_energy_balance(model::UtahEnergyBalance, state::SnowState, forcing, timestep)
    (; shortwave_radiation, longwave_radiation, air_temperature,
       wind_speed, relative_humidity, atmospheric_pressure) = forcing
    (; water_equivalent, temperature, albedo, density) = state

    zero_flux = 0.0u"W/m^2"
    zero_melt = zero(water_equivalent)

    if water_equivalent <= zero_melt
        return (
            net_radiation=zero_flux,
            sensible_heat=zero_flux,
            latent_heat=zero_flux,
            ground_heat=zero_flux,
            total_flux=zero_flux,
            melt=zero_melt,
            sublimation=zero_melt,
        )
    end

    # Convert temperatures to K (Unitful handles °C offset correctly)
    T_snow_K = uconvert(u"K", temperature)
    T_air_K = uconvert(u"K", air_temperature)

    ε = model.snow_emissivity
    z0 = model.snow_roughness_height

    # Net shortwave radiation
    Q_sw = (1.0 - albedo) * shortwave_radiation

    # Net longwave radiation
    Q_lw_in = longwave_radiation
    Q_lw_out = ε * STEFAN_BOLTZMANN * T_snow_K^4
    Q_lw = Q_lw_in - Q_lw_out

    # Net radiation
    Q_net = Q_sw + Q_lw

    # Turbulent fluxes using bulk aerodynamic approach
    z_ref = 2.0u"m"

    # Aerodynamic resistance (neutral stability approximation)
    u = max(wind_speed, 0.1u"m/s")
    r_a = (log(z_ref / z0))^2 / (VON_KARMAN^2 * u)

    # Sensible heat flux (positive = warming snow)
    H = DENSITY_AIR * SPECIFIC_HEAT_AIR * (T_air_K - T_snow_K) / r_a

    # Latent heat flux (sublimation/deposition) - simplified
    T_snow_C = ustrip(u"K", T_snow_K) - 273.15
    T_air_C = ustrip(u"K", T_air_K) - 273.15
    e_sat_snow = 611.2u"Pa" * exp(17.67 * T_snow_C / (T_snow_K - 29.65u"K"))
    e_sat_air = 611.2u"Pa" * exp(17.67 * T_air_C / (T_air_K - 29.65u"K"))
    e_air = relative_humidity * e_sat_air

    # Water vapor density gradient
    ρ_v_snow = e_sat_snow / (GAS_CONSTANT_WATER_VAPOR * T_snow_K)
    ρ_v_air = e_air / (GAS_CONSTANT_WATER_VAPOR * T_air_K)

    E = (ρ_v_air - ρ_v_snow) / r_a  # kg/m²/s (positive = deposition)
    L_E = E * LATENT_HEAT_SUBLIMATION

    # Ground heat flux (input, default 0 when not known per Tarboton & Luce 1996)
    Q_ground = forcing.ground_heat_flux

    # Total energy flux (equation 1 from Tarboton & Luce 1996)
    # dU/dt = Qsn + Qli + Qp + Qg - Qle + Qh + Qe - Qm
    # Here Q_net includes Qsn + Qli - Qle, H is Qh, L_E is Qe
    Q_total = Q_net + H + L_E + Q_ground

    # Sublimation mass flux (equation 2: dW/dt includes -E)
    # E is in kg/m²/s, convert to m water equivalent over timestep
    # Negative E means sublimation (mass loss), positive means deposition
    sublimation = -E * timestep / DENSITY_WATER
    sublimation = max(zero_melt, min(sublimation, water_equivalent))  # clamp to available

    # Compute energy content after this timestep's fluxes
    energy_content = state.energy_content
    energy_change = Q_total * timestep
    new_energy = energy_content + energy_change
    swe_after_sublimation = water_equivalent - sublimation

    # Melt outflow using composable formula (default: Darcy's law, Eq 47-48)
    mof = MeltOutflowForcing(
        air_temperature=air_temperature,
        snow_temperature=temperature,
        water_equivalent=swe_after_sublimation,
        density=density,
        energy_content=new_energy,
        cold_content=state.cold_content,
        liquid_water=state.liquid_water,
        net_energy_flux=Q_total,
        rainfall=zero(water_equivalent),  # UEB doesn't use rainfall in Darcy outflow
        day_of_year=forcing.day_of_year,
        timestep=timestep,
    )
    melt = melt_outflow(model.melt_outflow_formula, mof)

    return (;
        net_radiation=Q_net,
        sensible_heat=H,
        latent_heat=L_E,
        ground_heat=Q_ground,
        total_flux=Q_total,
        melt,
        sublimation,
    )
end

# ============================================================================
# KearneySnow - NicheMapR implementation
# ============================================================================

"""
    KearneySnow(; rain_melt_factor=0.0125, ...)

Snow model matching the NicheMapR microclimate model exactly.

This is a Julia port of the snow algorithms from NicheMapR (Kearney 2020),
implementing the same equations and default parameters for compatibility
and validation against the R/Fortran implementation.

Key algorithms:
- Sturm et al. (2010) density evolution with depth and age
- Aggarwal (2009) thermal conductivity from density
- Anderson (2006) albedo decay
- Heat-budget driven melt with rain-on-snow enhancement

# Parameters
- `rain_melt_factor`: Rain-on-snow melt coefficient (dimensionless)
  (NicheMapR default: rainmelt=0.0125)
- `maximum`: Maximum snow density
  (NicheMapR densfun[1]=0.5979 Mg/m³)
- `fresh`: Fresh snow density
  (NicheMapR densfun[2]=0.2178 Mg/m³)
- `density_depth_coefficient`: k1 in Sturm formula (1/cm)
  (NicheMapR densfun[3]=0.001)
- `density_age_coefficient`: k2 in Sturm formula (1/day)
  (NicheMapR densfun[4]=0.0038)
- `minimum_snow_depth`: Minimum snow depth for numerical stability
  (NicheMapR minsnow=2.5 cm)
- `undercatch`: Undercatch multiplier for converting rainfall to snow
- `snowmelt_fraction`: Proportion of calculated snowmelt that doesn't refreeze
- `accumulation`: Snow accumulation method (default: ThresholdAccumulation at 1.5°C)

# References
- Kearney, M.R. (2020). How will snow alter exposure of organisms to cold stress
  under climate warming? Global Ecology and Biogeography, 29, 1246-1256.
- Sturm, M. et al. (2010). Estimating snow water equivalent using snow depth data
  and climate classes. Journal of Hydrometeorology, 11, 1380-1394.
- Aggarwal, R. (2009). Thermal conductivity of snow. Defence Science Journal, 59, 126-130.
- Anderson, E.A. (2006). Snow Accumulation and Ablation Model – SNOW-17.
"""
Base.@kwdef struct KearneySnow{MF<:MeltOutflowFormula,UC,SMF,DF<:SnowDensityFormula,AF<:AlbedoFormula,TF<:SnowTemperatureFormula,A<:SnowAccumulation} <: SnowModel
    melt_outflow_formula::MF = EnergyFluxMeltOutflow()
    undercatch::UC = 1.0
    snowmelt_fraction::SMF = 1.0
    density_formula::DF = SturmSnowDensity()
    albedo_formula::AF = AndersonAlbedo()
    temperature_formula::TF = RelaxationTemperature(use_energy_flux=true)
    accumulation::A = ThresholdAccumulation(threshold=1.5u"°C")
end

"""
    snow_melt(model::KearneySnow, ...)

Compute snowmelt using Kearney/NicheMapR approach.

Delegates to `melt_outflow` with the model's `melt_outflow_formula`.
"""
function snow_melt(model::KearneySnow, state::SnowState, forcing, timestep)
    (; air_temperature, precipitation, net_energy_flux, day_of_year) = forcing

    # Compute rainfall (liquid precipitation after partitioning)
    rainfall = precipitation - snow_accumulation(model, precipitation, air_temperature).snowfall

    # Create forcing struct for melt outflow
    mof = MeltOutflowForcing(
        air_temperature=air_temperature,
        snow_temperature=state.temperature,
        water_equivalent=state.water_equivalent,
        density=state.density,
        energy_content=state.energy_content,
        cold_content=state.cold_content,
        liquid_water=state.liquid_water,
        net_energy_flux=net_energy_flux,
        rainfall=rainfall,
        day_of_year=day_of_year,
        timestep=timestep,
    )
    return melt_outflow(model.melt_outflow_formula, mof)
end

# ============================================================================
# Refreezing
# ============================================================================

"""
    snow_refreeze(model, state)

Refreeze liquid water in snowpack when cold content is positive.

Returns updated state with liquid water converted back to ice.
Used by Snow17 which tracks liquid water storage and cold content.
"""
function snow_refreeze(::Snow17, state)
    if state.cold_content > 0.0u"mm" && state.liquid_water > 0.0u"m"
        refreeze = min(state.liquid_water, state.cold_content)
        @set! state.liquid_water -= refreeze
        @set! state.water_equivalent += refreeze
        @set! state.cold_content -= refreeze
    end
    return state
end

# ============================================================================
# State update functions
# ============================================================================

"""
    update_snow_state(model::SnowModel, state::SnowState, forcing, timestep)

Update snow state for one timestep.

Handles accumulation, melt, density evolution, and albedo decay.
Returns `(new_state, water_output)` where water_output is melt + rain for routing to surface water.

# Arguments
- `model`: Snow model (`NoSnow`, `GenericSnowModel`, `Snow17`, `UtahEnergyBalance`, `KearneySnow`)
- `state`: Current snow state (`SnowState`)
- `forcing`: Atmospheric forcing data. Can be a `SnowForcing` struct or any object with
  the required fields. All models require `air_temperature` and `precipitation`.
  Energy balance models additionally require radiation, wind, and humidity fields.
- `timestep`: Time step duration

# Returns
- `(new_state, water_output)`: Updated snow state and water available for runoff (melt + rain)

See also: [`SnowForcing`](@ref), [`SnowState`](@ref)
"""
function update_snow_state end

# NoSnow: returns precipitation unchanged
function update_snow_state(::NoSnow, state::SnowState, forcing, timestep)
    return state, forcing.precipitation
end

function update_snow_state(model::GenericSnowModel, state::SnowState, forcing, timestep)
    (; precipitation, air_temperature) = forcing

    # Partition precipitation
    (; snowfall, rainfall) = snow_accumulation(model, precipitation, air_temperature)

    # Add snowfall to pack
    if snowfall > zero(snowfall)
        @set! state.water_equivalent += snowfall
        @set! state.age = 0.0u"d"
    end

    # Update albedo, density and depth
    @set! state.albedo = snow_albedo(model.albedo_formula, state.age)
    @set! state.density = snow_density(model.density_formula, state, snowfall, timestep, air_temperature)
    if state.water_equivalent > zero(state.water_equivalent)
        @set! state.depth = snow_depth(state.water_equivalent, state.density)
    end

    # Compute melt
    melt = snow_melt(model, state, forcing, timestep)
    outflow = melt

    # Remove melt from pack
    @set! state.water_equivalent = max(zero(state.water_equivalent), state.water_equivalent - melt)
    if state.water_equivalent > zero(state.water_equivalent)
        @set! state.depth = snow_depth(state.water_equivalent, state.density)
    else
        @set! state.depth = 0.0u"m"
        @set! state.density = model.density_formula.fresh
        @set! state.albedo = model.albedo_formula.fresh
        @set! state.age = 0.0u"d"
    end

    @set! state.age += timestep
    @set! state.temperature = snow_temperature(model.temperature_formula, state, forcing, melt, timestep)

    return state, outflow + rainfall
end

function update_snow_state(model::UtahEnergyBalance, state::SnowState, forcing, timestep)
    (; precipitation, air_temperature, shortwave_radiation, longwave_radiation,
       wind_speed, relative_humidity, atmospheric_pressure) = forcing

    # Partition precipitation
    (; snowfall, rainfall) = snow_accumulation(model, precipitation, air_temperature)

    # Add snowfall and its energy content (Qp in equation 1)
    if snowfall > zero(snowfall)
        @set! state.water_equivalent += snowfall
        @set! state.age = 0.0u"d"

        # Energy content of snowfall (cold snow adds negative energy)
        T_snow_C = min(0.0, ustrip(u"°C", air_temperature))
        snowfall_energy = snowfall * DENSITY_WATER * SPECIFIC_HEAT_ICE * T_snow_C * u"K"
        @set! state.energy_content += snowfall_energy
    end

    # Update density and depth
    @set! state.density = snow_density(model.density_formula, state, snowfall, timestep, air_temperature)
    if state.water_equivalent > zero(state.water_equivalent)
        @set! state.depth = snow_depth(state.water_equivalent, state.density)
    end

    # Update albedo (handles both Dickinson grain-growth and simple decay formulas)
    albedo, new_surface_age = compute_albedo(model.albedo_formula, state, snowfall, timestep)
    @set! state.albedo = albedo
    @set! state.surface_age = new_surface_age

    # Energy balance (equation 1: dU/dt) and mass balance (equation 2: dW/dt)
    eb = snow_energy_balance(model, state, forcing, timestep)

    # Update energy content (equation 1)
    @set! state.energy_content += eb.total_flux * timestep

    # Mass balance: subtract sublimation and melt (equation 2)
    # Sublimation removes mass without phase change energy cost to snowpack
    # (latent heat is already in the energy balance as Qe)
    sublimation = eb.sublimation
    melt = eb.melt
    outflow = melt  # sublimation goes to atmosphere, not outflow

    # Remove melt energy from energy content
    melt_energy = melt * DENSITY_WATER * LATENT_HEAT_FUSION
    @set! state.energy_content -= melt_energy

    # Update water equivalent (equation 2: dW/dt = Pr + Ps - Mr - E)
    @set! state.water_equivalent = max(zero(state.water_equivalent),
                                        state.water_equivalent - melt - sublimation)

    # Update depth after melt
    if state.water_equivalent > zero(state.water_equivalent)
        @set! state.depth = snow_depth(state.water_equivalent, state.density)
    else
        # No snow - reset state
        @set! state.depth = 0.0u"m"
        @set! state.density = model.density_formula.fresh
        @set! state.albedo = fresh_snow_albedo(model.albedo_formula)
        @set! state.age = 0.0u"d"
        @set! state.surface_age = 0.0
        @set! state.energy_content = 0.0u"J/m^2"
    end

    @set! state.age += timestep

    # Temperature derived from energy content (equations 3-5)
    @set! state.temperature = snow_temperature(model.temperature_formula, state, forcing, melt, timestep)

    return state, outflow + rainfall
end

function update_snow_state(model::Snow17, state::SnowState, forcing, timestep)
    (; precipitation, air_temperature, day_of_year) = forcing

    # Partition precipitation
    (; snowfall, rainfall) = snow_accumulation(model, precipitation, air_temperature)

    # Apply snowfall correction factor for gage undercatch (Anderson Eq 1)
    snowfall = snowfall * model.snowfall_correction_factor

    T_air_C = ustrip(u"°C", air_temperature)
    timestep_hr = ustrip(u"hr", timestep)
    snowfall_mm = ustrip(u"mm", snowfall)

    # Update Antecedent Temperature Index (Anderson Equation 8)
    # ATI₂ = ATI₁ + TIPM^Δt × (Ta - ATI₁)
    # ATI is capped at 0°C and reset to new snow temperature on significant snowfall
    ATI_C = ustrip(u"°C", state.antecedent_temperature_index)
    TIPM = model.temperature_index_parameter

    # Significant snowfall threshold: 1.5 mm per 6 hours, scaled to timestep
    significant_snowfall = 1.5 * (timestep_hr / 6.0)

    if snowfall_mm > significant_snowfall
        # Reset ATI to temperature of new snow
        ATI_C = min(0.0, T_air_C)
    else
        # Exponential weighting: TIPM^(Δt/6hr) gives proper time scaling
        tipm_scaled = TIPM^(timestep_hr / 6.0)
        ATI_C = ATI_C + tipm_scaled * (T_air_C - ATI_C)
        # ATI cannot exceed 0°C
        ATI_C = min(0.0, ATI_C)
    end
    @set! state.antecedent_temperature_index = ATI_C * u"°C"

    # Add snowfall to pack
    if snowfall > zero(snowfall)
        @set! state.water_equivalent += snowfall
        @set! state.age = 0.0u"d"

        # Increase cold content from cold snowfall (Anderson Equation 4)
        # ΔDp = -(Tn × Pn) / (Lf/ci) where Lf/ci = 160 K
        T_snow_C = min(0.0, T_air_C)
        if T_snow_C < 0.0
            cold_content_increase = -T_snow_C * snowfall_mm / 160.0
            @set! state.cold_content += cold_content_increase * u"mm"
        end
    end

    # Update albedo, density and depth
    @set! state.albedo = snow_albedo(model.albedo_formula, state.age)
    @set! state.density = snow_density(model.density_formula, state, snowfall, timestep, air_temperature)
    if state.water_equivalent > zero(state.water_equivalent)
        @set! state.depth = snow_depth(state.water_equivalent, state.density)
    end

    # Update cold content based on temperature gradient (Anderson Equation 9)
    # ΔDt = NMf × (ATI - Tsur)
    # where NMf = NMF × (Δtp/6) × (Mf/MFMAX)
    # and Tsur = min(Ta, 0°C) is snow surface temperature
    T_surface = min(0.0, T_air_C)
    mf = seasonal_melt_factor(model, day_of_year)

    # Negative melt factor varies seasonally like the melt factor (Eq 9)
    # NMf = NMF × (Δtp/6) × (Mf/MFMAX)
    nmf = ustrip(u"mm/K/(6hr)", model.negative_melt_factor)
    mfmax = ustrip(u"mm/K/(6hr)", model.maximum_melt_factor)
    mf_ratio = ustrip(u"mm/K/(6hr)", mf) / mfmax
    nmf_scaled = nmf * (timestep_hr / 6.0) * mf_ratio

    # Change in cold content from temperature gradient using ATI
    # ΔDt = NMf × (ATI - Tsur)
    if state.water_equivalent > zero(state.water_equivalent)
        delta_cold = nmf_scaled * (ATI_C - T_surface)
        new_cold_content = ustrip(u"mm", state.cold_content) + delta_cold
        @set! state.cold_content = max(0.0, new_cold_content) * u"mm"
    end

    # Compute melt
    melt = snow_melt(model, state, forcing, timestep)

    # Handle liquid water storage (Snow17 tracks this)
    # Liquid water capacity is based on ICE portion only (Anderson Equation 10)
    # Wqx = PLWHC × Wi
    # Note: state.water_equivalent here represents ice (Wi), liquid_water is separate
    ice_portion = state.water_equivalent  # Wi
    liquid_water_capacity = model.liquid_water_holding_capacity * ice_portion

    # Get current lagged storage
    current_storage_mm = ustrip(u"mm", state.lagged_excess_storage)
    ice_mm = ustrip(u"mm", state.water_equivalent)

    if melt > zero(melt)
        # Add melt to liquid water storage
        @set! state.liquid_water += melt

        # Compute excess water above liquid water capacity
        if state.liquid_water > liquid_water_capacity
            excess_water = state.liquid_water - liquid_water_capacity
            @set! state.liquid_water = liquid_water_capacity
            excess_water_mm = ustrip(u"mm", excess_water)
        else
            excess_water_mm = 0.0
        end

        # Route excess water through lag/attenuation (Anderson Eq 21-24)
        outflow_mm, new_storage_mm = route_lagged_water(
            excess_water_mm, ice_mm, current_storage_mm, timestep_hr
        )
        outflow = outflow_mm * u"mm"
        @set! state.lagged_excess_storage = new_storage_mm * u"mm"

        # Remove melt from ice portion
        @set! state.water_equivalent = max(zero(state.water_equivalent), state.water_equivalent - melt)
    else
        # No new melt, but still drain from lagged storage
        outflow_mm, new_storage_mm = route_lagged_water(
            0.0, ice_mm, current_storage_mm, timestep_hr
        )
        outflow = outflow_mm * u"mm"
        @set! state.lagged_excess_storage = new_storage_mm * u"mm"

        state = snow_refreeze(model, state)
    end

    # Update albedo, density and depth
    if state.water_equivalent > zero(state.water_equivalent)
        @set! state.depth = snow_depth(state.water_equivalent, state.density)
    else
        # No snow - reset all state
        @set! state.depth = 0.0u"m"
        @set! state.density = model.density_formula.fresh
        @set! state.albedo = model.albedo_formula.fresh
        @set! state.liquid_water = 0.0u"m"
        @set! state.cold_content = 0.0u"mm"
        @set! state.antecedent_temperature_index = 0.0u"°C"
        @set! state.lagged_excess_storage = 0.0u"mm"
    end

    @set! state.age += timestep
    @set! state.temperature = snow_temperature(model.temperature_formula, state, forcing, melt, timestep)

    return state, outflow + rainfall
end

function update_snow_state(model::KearneySnow, state::SnowState, forcing, timestep)
    (; precipitation, air_temperature, net_energy_flux) = forcing

    # Partition precipitation
    (; snowfall, rainfall) = snow_accumulation(model, precipitation, air_temperature)

    # Add snowfall
    if snowfall > zero(snowfall)
        @set! state.water_equivalent += snowfall
        @set! state.age = 0.0u"d"
    end

    # Update albedo, density and depth
    @set! state.albedo = snow_albedo(model.albedo_formula, state.age)
    if state.water_equivalent > zero(state.water_equivalent)
        @set! state.density = snow_density(model.density_formula, state, snowfall, timestep, air_temperature)
        @set! state.depth = snow_depth(state.water_equivalent, state.density)
    end

    # Compute melt
    melt = snow_melt(model, state, forcing, timestep)

    if melt > zero(melt)
        @set! state.water_equivalent = max(zero(state.water_equivalent), state.water_equivalent - melt)
    end

    outflow = melt

    # Update albedo, density and depth
    if state.water_equivalent > zero(state.water_equivalent)
        @set! state.depth = snow_depth(state.water_equivalent, state.density)
    else
        @set! state.depth = 0.0u"m"
        @set! state.density = model.density_formula.fresh
        @set! state.albedo = model.albedo_formula.fresh
        @set! state.age = 0.0u"d"
    end

    @set! state.age += timestep
    @set! state.temperature = snow_temperature(model.temperature_formula, state, forcing, melt, timestep)

    return state, outflow + rainfall
end
