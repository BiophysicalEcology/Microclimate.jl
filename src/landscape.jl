Base.@kwdef struct SoilLayers
    depp::Vector{typeof(0.0u"m")}
    wc::Vector{typeof(1.0u"J/K/m^2")}
    c::Vector{typeof(1.0u"W/K/m^2")}
end

Base.@kwdef struct MoistLayers
    P::Vector{typeof(1.0u"J/kg")} 
    Z::Vector{typeof(1.0u"m")}
    V::Vector{typeof(1.0u"kg/m^2")}
    W::Vector{typeof(1.0u"m^3/m^3")}
    WN::Vector{typeof(1.0u"m^3/m^3")}
    K::Vector{typeof(1.0u"kg*s/m^3")}
    H::Vector{Float64}
    T::Vector{typeof(1.0u"K")}
    rh_soil::Vector{Float64}
    ψ_soil::Vector{typeof(1.0u"J/kg")}
    ψ_root::Vector{typeof(1.0u"J/kg")}
    PR::Vector{typeof(1.0u"J/kg")}
    PP::Vector{typeof(1.0u"J/kg")}
    B1::Vector{Float64}
    N::Vector{Float64}
    N1::Vector{Float64}
    WS::Vector{Float64}
    RR::Vector{typeof(1.0u"m^4/kg/s")}
    BZ::Vector{typeof(1.0u"m")}
    JV::Vector{typeof(1.0u"kg/m^2/s")}
    DJ::Vector{typeof(1.0u"kg*s/m^4")}
    CP::Vector{typeof(1.0u"kg*s/m^4")}
    A::Vector{typeof(1.0u"kg*s/m^4")}
    B::Vector{typeof(1.0u"kg*s/m^4")}
    C::Vector{typeof(1.0u"kg*s/m^4")}
    C2::Vector{Float64}
    F::Vector{typeof(1.0u"kg/m^2/s")}
    F2::Vector{typeof(1.0u"J/kg")}
    DP::Vector{typeof(1.0u"J/kg")}
end


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

Base.@kwdef struct MicroInputs
    params::MicroParams
    forcing::MicroForcing
    soillayers::SoilLayers
end

#Base.@kwdef struct MicroInputs{MP::MicroParams,MF<:MicroForcing,SL<:SoilLayers}
#    params::MP
#    forcing::MF
#    soillayers::SL
#end

function init_soillayers(N)
    depp = fill(0.0u"cm", N + 1)
    wc = fill(1.0u"J/K/m^2", N)
    c = fill(1.0u"W/K/m^2", N)
    return SoilLayers(depp, wc, c)
end
function init_moistlayers(M)
    MoistLayers(
        P = zeros(M+1) .* u"J/kg",
        Z = zeros(M+1) .* u"m",
        V = zeros(M+1) .* u"kg/m^2",
        W = zeros(M+1) .* u"m^3/m^3",
        WN = zeros(M+1) .* u"m^3/m^3",
        K = zeros(M+1) .* u"kg*s/m^3",
        H = zeros(M+1),
        T = zeros(M+1) .* u"K",
        rh_soil = zeros(M),
        ψ_soil = zeros(M) .* u"J/kg",
        ψ_root = zeros(M) .* u"J/kg",
        PR = zeros(M+1) .* u"J/kg",
        PP = zeros(M+1) .* u"J/kg",
        B1 = zeros(M+1),
        N = zeros(M+1),
        N1 = zeros(M+1),
        WS = zeros(M+1),
        RR = zeros(M+1) .* u"m^4/kg/s",
        BZ = zeros(M+1) .* u"m",
        JV = zeros(M+1) .* u"kg/m^2/s",
        DJ = zeros(M+1) .* u"kg*s/m^4",
        CP = zeros(M+1) .* u"kg*s/m^4",
        A = zeros(M+1) .* u"kg*s/m^4",
        B = zeros(M+1) .* u"kg*s/m^4",
        C = zeros(M+1) .* u"kg*s/m^4",
        C2 = zeros(M+1),
        F = zeros(M+1) .* u"kg/m^2/s",
        F2 = zeros(M+1) .* u"J/kg",
        DP = zeros(M+1) .* u"J/kg"
    )
end