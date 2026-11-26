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

# Cosby regresses b itself, and the logarithms of air-entry potential (cm) and conductivity
_cosby_parameters(b_parameter, air_entry_exponent, conductivity_exponent) = (;
    campbell_b_parameter = b_parameter,
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
    # Representative diameters (mm) of the clay, silt and sand classes. Sand is
    # NicheMapR's 1.05, not the class mean of 1.025
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

"""
    pedotransfer_soil_profile(model::AbstractPedotransferModel, depths; texture_depths,
        clay_percent, silt_percent, sand_percent, bulk_density, root_density, soil_properties...)

A [`SoilProfile`](@ref) on `depths` from texture and bulk density sampled at `texture_depths`, with
Campbell hydraulic parameters from [`pedotransfer`](@ref). The samples are interpolated linearly
between `texture_depths` and held constant beyond them. `root_density` defaults to
`example_campbell_hydraulic_profile`'s, interpolated the same way. The other keywords go to
`example_soil_profile`.
"""
function pedotransfer_soil_profile(model::AbstractPedotransferModel, depths;
    texture_depths, clay_percent, silt_percent, sand_percent, bulk_density,
    root_density = onto_depths(DEFAULT_DEPTHS, example_campbell_hydraulic_profile().root_density, depths),
    soil_properties...,
)
    from_samples(values) = onto_depths(texture_depths, values, depths)
    bulk_density = from_samples(bulk_density)
    hydraulics = CampbellHydraulicProfile(model; clay_percent = from_samples(clay_percent),
        silt_percent = from_samples(silt_percent), sand_percent = from_samples(sand_percent), bulk_density, root_density)
    return example_soil_profile(depths; bulk_density, hydraulics, soil_properties...)
end

# Values known at `known_depths`, linear between them and constant beyond them, at `depths`
function onto_depths(known_depths, values, depths)
    unit = oneunit(first(values))
    profile = extrapolate(interpolate((ustrip.(u"m", known_depths),), values ./ unit, Gridded(Interpolations.Linear())), Flat())
    return profile.(ustrip.(u"m", depths)) .* unit
end

"""
    international_to_usda_texture(; clay_percent, silt_percent, sand_percent)

Texture in the International Society of Soil Science (Atterberg) particle-size classes, silt 2-20 µm
and sand 20-2000 µm, as in the Soil and Landscape Grid of Australia, converted to the United States Department of Agriculture (USDA)
classes the pedotransfer models take, silt 2-50 µm and sand 50-2000 µm. The cumulative particle-size
distribution is taken as linear in log diameter across the sand class (Nemes et al. 1999), so
log(50/20) / log(2000/20) of the sand moves to silt.

# References
- Nemes, A., Wösten, J. H. M., Lilly, A. & Oude Voshaar, J. H. (1999). Evaluation of different
  procedures to interpolate particle-size distributions to achieve compatibility within soil
  databases. *Geoderma*, 90(3-4), 187-202.
"""
function international_to_usda_texture(; clay_percent, silt_percent, sand_percent)
    fine_sand = sand_percent .* (log(50 / 20) / log(2000 / 20))
    return (; clay_percent, silt_percent = silt_percent .+ fine_sand, sand_percent = sand_percent .- fine_sand)
end
