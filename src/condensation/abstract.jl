abstract type AbstractCondensationModel end

"""
    condensation_energy_flux(model::AbstractCondensationModel; kwargs...)

Signed surface energy flux (W/m²) available for dew/frost formation this
hour -- Monteith convention, positive = evaporative loss, negative =
condensation. Ground variants (on `MicroModel.condensation_model`):
[`NoCondensation`](@ref) (off), [`BulkTransferCondensation`](@ref),
[`GarrattSegalCondensation`](@ref) (default). Leaf variant (on
`MultilayerCanopy.condensation_model`): [`MonteithLeafCondensation`](@ref).
Each expects different keyword arguments -- ground and leaf variants are
not interchangeable.
"""
function condensation_energy_flux end
