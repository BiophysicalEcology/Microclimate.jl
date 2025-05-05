Base.@kwdef struct MicroParams
    soilprops::Matrix{Union{Unitful.AbstractQuantity, Float64}}
    dep::Vector{<:Unitful.AbstractQuantity}
    refhyt::Quantity
    ruf::Quantity
    d0::Quantity
    zh::Quantity
    slope::Quantity
    shade::Float64
    viewf::Float64
    elev::Quantity
    refl::Float64
    sle::Float64
    slep::Float64
    pctwet::Float64
    nodes::Vector{Float64}
    tdeep::Quantity
    θ_soil::Vector{Float64}
    runmoist::Bool
    runsnow::Bool
end

Base.@kwdef struct MicroForcing
    SOLRt::Interpolations.AbstractInterpolation
    ZENRt::Interpolations.AbstractInterpolation
    ZSLt::Interpolations.AbstractInterpolation
    TAIRt::Interpolations.AbstractInterpolation
    VELt::Interpolations.AbstractInterpolation
    RHt::Interpolations.AbstractInterpolation
    CLDt::Interpolations.AbstractInterpolation
end

Base.@kwdef struct MicroInput
    params::MicroParams
    forcing::MicroForcing
end