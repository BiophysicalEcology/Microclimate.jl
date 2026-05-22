"""
    KearneyShortwave()

NicheMapR's (Kearney et al.) shortwave budget at the surface:
- Ångström scaling: global = sunshine_fraction(sunshine_fraction_model, cloud) × (diffuse_clear_sky + direct_clear_sky)
- Diffuse/direct split via a Maxwell (1987) clearness index passed to the diffuse-fraction model.

References
- Maxwell, E. L. (1987). A Quasi-Physical Model for Converting Hourly Global
  Horizontal to Direct Normal Insolation. SERI/TR-215-3087. Golden, CO.
"""
struct KearneyShortwave <: AbstractShortwaveScheme end

function shortwave_radiation!(::KearneyShortwave, output, cloud::AbstractArray, diffuse_clear_sky, direct_clear_sky, zenith::AbstractArray, doy;
    diffuse_fraction_model::AbstractDiffuseFractionModel=ErbsDiffuseFraction(),
    sunshine_fraction_model::AbstractSunshineFractionModel=Angstrom(),
)
    (; global_horizontal) = output.solar_radiation
    global_radiation = global_horizontal
    diffuse_fraction = output.diffuse_fraction
    solar_constant = 1367.0u"W/m^2"
    ϵ = 1e-9u"W/m^2"
    @inbounds for i in eachindex(global_radiation)
        # 1) Extraterrestrial horizontal irradiance for this hour's day
        d = doy[i]
        θ1 = 360.0 * (d - 1) / 365.0
        θ2 = 2.0 * θ1
        ec = 1.00011 + 0.034221*cosd(θ1) + 0.00128*sind(θ1) +
                       0.000719*cosd(θ2) + 0.000077*sind(θ2)
        eth = solar_constant * ec
        # 2) scale by sunshine fraction (cloud-cover adjustment)
        sf = sunshine_fraction(sunshine_fraction_model, cloud[i])
        gcs = diffuse_clear_sky[i] + direct_clear_sky[i]
        gr = max(sf * gcs, 0.0u"W/m^2")
        global_radiation[i] = gr
        # 3) Split global into diffuse/direct via clearness index
        ci = clamp(gr / max(eth, ϵ), 0.0, 1.2)
        diffuse_fraction[i] = clamp(calc_diffuse_fraction(diffuse_fraction_model, ci), 0.0, 1.0)
    end
    return (; global_radiation, diffuse_fraction)
end
