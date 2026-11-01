"""
    AbstractPedotransferModel

Supertype for soil texture to Campbell hydraulic parameter regressions:
[`CosbyUnivariate`](@ref), [`CosbyMultivariate`](@ref) or [`CampbellParticleSize`](@ref).
Variants implement [`pedotransfer`](@ref).
"""
abstract type AbstractPedotransferModel end

"""
    CosbyUnivariate()

Cosby et al. (1984) Table 5: `campbell_b_parameter` from clay, air-entry potential
and saturated conductivity from sand. Ignores bulk density.

# References
- Cosby, B. J., Hornberger, G. M., Clapp, R. B. & Ginn, T. R. (1984). A statistical exploration
  of the relationships of soil moisture characteristics to the physical properties of soils.
  *Water Resources Research*, 20(6), 682-690.
"""
struct CosbyUnivariate <: AbstractPedotransferModel end

"""
    CosbyMultivariate()

Cosby et al. (1984) Table 4 multivariate regressions on clay, silt and sand. Ignores bulk density.
"""
struct CosbyMultivariate <: AbstractPedotransferModel end

"""
    CampbellParticleSize()

Campbell (1985, ch. 5) from the geometric mean and standard deviation of the
particle-size distribution, with air-entry potential and saturated conductivity
scaled from their values at a bulk density of 1.3 Mg/m³.

# References
- Campbell, G. S. (1985). *Soil Physics with BASIC*. Elsevier.
"""
struct CampbellParticleSize <: AbstractPedotransferModel end

"""
    pedotransfer(model::AbstractPedotransferModel, clay_percent, silt_percent, sand_percent, bulk_density)

Campbell hydraulic parameters of a soil from its texture (percent clay <2 µm,
silt 2-50 µm and sand >50 µm) and `bulk_density`, as
`(; campbell_b_parameter, air_entry_water_potential, saturated_hydraulic_conductivity)`,
named as in [`CampbellHydraulicProfile`](@ref). `air_entry_water_potential` is a positive magnitude.
"""
function pedotransfer end

const _COSBY_CONDUCTIVITY = 0.0007196666u"kg*s/m^3"  # Cosby's inches/hour, to kg s m⁻³

_cosby_parameters(b_exponent, air_entry_exponent, conductivity_exponent) = (;
    campbell_b_parameter = log10(10^b_exponent / 10.2),
    air_entry_water_potential = 10^air_entry_exponent / 10.2 * u"J/kg",
    saturated_hydraulic_conductivity = 10^conductivity_exponent * _COSBY_CONDUCTIVITY,
)

pedotransfer(::CosbyUnivariate, clay_percent, silt_percent, sand_percent, bulk_density) =
    _cosby_parameters(0.1590 * clay_percent + 2.910, -0.0131 * sand_percent + 1.880, 0.0153 * sand_percent - 0.884)
pedotransfer(::CosbyMultivariate, clay_percent, silt_percent, sand_percent, bulk_density) =
    _cosby_parameters(0.1570 * clay_percent - 0.0030 * sand_percent + 3.10,
        -0.0095 * sand_percent + 0.0063 * silt_percent + 1.54,
        0.0126 * sand_percent - 0.0064 * clay_percent - 0.60)

function pedotransfer(::CampbellParticleSize, clay_percent, silt_percent, sand_percent, bulk_density)
    # Representative diameters (mm) of the clay, silt and sand classes
    fractions = (clay_percent, silt_percent, sand_percent) ./ 100
    log_diameters = log.((0.001, 0.026, 1.05))
    log_mean = sum(fractions .* log_diameters)
    log_std = sqrt(sum(fractions .* log_diameters .^ 2) - log_mean^2)
    air_entry_at_reference = 0.5 * exp(log_mean)^(-1 / 2)  # J/kg at 1.3 Mg/m³
    campbell_b_parameter = 2 * air_entry_at_reference + 0.2 * exp(log_std)
    density_ratio = ustrip(u"Mg/m^3", bulk_density) / 1.3
    return (;
        campbell_b_parameter,
        air_entry_water_potential = air_entry_at_reference * density_ratio^(0.67 * campbell_b_parameter) * u"J/kg",
        saturated_hydraulic_conductivity = 0.004u"kg*s/m^3" * density_ratio^(-1.3 * campbell_b_parameter) *
            exp(-6.9 * fractions[1] - 3.7 * fractions[2]),
    )
end

"""
    CampbellHydraulicProfile(model::AbstractPedotransferModel; clay_percent, silt_percent, sand_percent,
        bulk_density, root_density)

Per-depth Campbell hydraulic parameters from per-depth texture and bulk density through [`pedotransfer`](@ref).
"""
function CampbellHydraulicProfile(model::AbstractPedotransferModel;
    clay_percent, silt_percent, sand_percent, bulk_density, root_density,
)
    parameters = pedotransfer.(Ref(model), clay_percent, silt_percent, sand_percent, bulk_density)
    return CampbellHydraulicProfile(;
        air_entry_water_potential = [p.air_entry_water_potential for p in parameters],
        saturated_hydraulic_conductivity = [p.saturated_hydraulic_conductivity for p in parameters],
        campbell_b_parameter = [p.campbell_b_parameter for p in parameters],
        root_density,
    )
end
