"""
    ErbsDiffuseFraction()

Estimate the diffuse fraction of global solar radiation from the clearness index
using the piecewise model of Erbs et al. (1982).

# References
Erbs, D. G., Klein, S. A., & Duffie, J. A. (1982). Estimation of the diffuse
radiation fraction for hourly, daily and monthly-average global radiation.
*Solar Energy*, 28(4), 293–302.
"""
struct ErbsDiffuseFraction <: AbstractDiffuseFractionModel end

function calc_diffuse_fraction(::ErbsDiffuseFraction, clearness_index)
    if clearness_index <= 0.22
        return 1 - 0.09 * clearness_index
    elseif clearness_index <= 0.80
        return 0.9511 - 0.1604 * clearness_index + 4.388 * clearness_index^2 -
               16.638 * clearness_index^3 + 12.336 * clearness_index^4
    else
        return 0.165
    end
end
