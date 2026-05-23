"""
    Angstrom(; a=0.36, b=0.64, gamma=1.0)

Ångström–Prescott cloud scaling: clear-sky radiation is multiplied by
`a + b * (1 - cloud_cover)^gamma`. The defaults match Linacre's "Climate Data
and Resources" (1992, eq. 5.33).
"""
@kwdef struct Angstrom{A,B,G} <: AbstractSunshineFractionModel
    a::A = 0.36
    b::B = 0.64
    gamma::G = 1.0
end

sunshine_fraction(m::Angstrom, cloud_cover) = m.a + m.b * (1.0 - cloud_cover)^m.gamma
