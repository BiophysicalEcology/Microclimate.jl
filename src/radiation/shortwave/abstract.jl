abstract type AbstractShortwaveScheme end

"""
    shortwave_radiation!(scheme, output, cloud, diffuse_clear_sky, direct_clear_sky, zenith, doy; sunshine_fraction_model, diffuse_fraction_model)

Full shortwave budget: starting from clear-sky diffuse and direct components,
apply cloud-cover scaling and split the resulting global radiation into a
diffuse fraction. Writes into `output.solar_radiation.global_horizontal` and
`output.diffuse_fraction`; returns `(; global_radiation, diffuse_fraction)`.
"""
function shortwave_radiation! end
