abstract type AbstractCloudAdjustModel end

"""
    sunshine_fraction(model, cloud_cover)

Fraction of clear-sky radiation that reaches the ground given cloud cover
(0–1). Combined with the model's a/b coefficients to give the Ångström–type
scaling.
"""
function sunshine_fraction end
