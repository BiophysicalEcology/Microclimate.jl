@compound H2O
@compound O2
@compound CO2
@compound N2
"""
    get_pressure(h::Quantity;
                 h_ref::Quantity = 0u"m",
                 P_ref::Quantity = 101325u"Pa",
                 L_ref::Quantity = -0.0065u"K/m",
                 T_ref::Quantity = 288u"K",
                 g_0::Quantity = 9.80665u"m/s^2",
                 M::Quantity = 0.0289644u"kg/mol") -> Quantity

Computes atmospheric pressure at a given altitude using the barometric formula,
assuming a constant temperature lapse rate (standard tropospheric approximation).

# Arguments
- `h`: Elevation at which to compute pressure (with length units, e.g. `u"m"`).
- `h_ref`: Reference altitude (default: `0u"m"`).
- `P_ref`: Pressure at the reference altitude (default: `101325u"Pa"`, standard sea-level pressure).
- `L_ref`: Temperature lapse rate (default: `-0.0065u"K/m"`).
- `T_ref`: Temperature at the reference altitude (default: `288u"K"`).
- `g_0`: Standard gravitational acceleration (default: `9.80665u"m/s^2"`).
- `M`: Molar mass of dry air (default: `0.0289644u"kg/mol"`).

# Returns
- `P_a`: Atmospheric pressure at altitude `h` (with pressure units, e.g. `u"Pa"`).

# Notes
- Uses the simplified barometric formula assuming a linear lapse rate and ideal gas behavior.
- The universal gas constant `R` is used from the `Unitful` package.

# Example
```julia
using Unitful

P = get_pressure(1500u"m")
"""
function get_pressure(h::Quantity;
    P_ref::Quantity = 101325u"Pa",
    L_ref::Quantity = -0.0065u"K/m",
    T_ref::Quantity = 288.0u"K",
    g_0::Quantity = 9.80665u"m/s^2",
    M::Quantity = 0.0289644u"kg/mol")
    R=Unitful.R
P_a = P_ref * (1 + (L_ref / T_ref) * h) ^ ((-g_0 * M) / (R * L_ref))#5.2553026003237262u"kg*m^2*J^-1*s^-2"#

return P_a
end

"""
    vapour_pressure(T)

Calculates saturation vapour pressure (Pa) for a given air temperature.

# Arguments
- `T`: air temperature in K.
"""
function vapour_pressure(T)
    T = Unitful.ustrip(T) + 0.01 # triple point of water is 273.16
    logP_vap = T
    if T <= 273.16
        logP_vap = -9.09718 * (273.16 / T - 1) - 3.56654 * log10(273.16 / T) + 0.876793 * (1 - T / 273.16) + log10(6.1071)
    else
        logP_vap = -7.90298 * (373.16 / T - 1) + 5.02808 * log10(373.16 / T) - 1.3816E-07 * (10^(11.344 * (1 - T / 373.16)) - 1) + 8.1328E-03 * (10^(-3.49149 * (373.16 / T - 1)) - 1) + log10(1013.246)
    end
    (10^logP_vap) * 100u"Pa"
end

"""
    wet_air(T_drybulb, T_wetbulb, rh, T_dew, P_atmos, fO2, fCO2, fN2)
    wet_air(T_drybulb; kw...)

Calculates several properties of humid air as output variables below. The program
is based on equations from List, R. J. 1971. Smithsonian Meteorological Tables. Smithsonian
Institution Press. Washington, DC. wet_air must be used in conjunction with function vapour_pressure.

Input variables are shown below. The user must supply known values for T_drybulb and P (P at one standard
atmosphere is 101 325 pascals). Values for the remaining variables are determined by whether the user has
either (1) psychrometric data (T_wetbulb or rh), or (2) hygrometric data (T_dew)

# TODO fix this desctiption
(1) Psychrometric data:
If T_wetbulb is known but not rh, then set rh=-1 and dp=999
If rh is known but not T_wetbulb then set T_wetbulb=0 and dp=999

(2) Hygrometric data:
If T_dew is known then set T_wetublb = 0 and rh = 0.

# Arguments

- `T_drybulb`: Dry bulb temperature (K)
- `T_wetbulb`: Wet bulb temperature (K)
- `rh`: Relative humidity (%)
- `T_dew`: Dew point temperature (K)
- `P`: Barometric pressure (Pa)
- `fO2`; fractional O2 concentration in atmosphere, -
- `fCO2`; fractional CO2 concentration in atmosphere, -
- `fN2`; fractional N2 concentration in atmosphere, -
# - `P_vap`: Vapour pressure (Pa)
# - `P_vap_sat`: Saturation vapour pressure (Pa)
# - `ρ_vap`: Vapour density (kg m-3)
# - `r_w Mixing`: ratio (kg kg-1)
# - `T_vir`: Virtual temperature (K)
# - `T_vinc`: Virtual temperature increment (K)
# - `ρ_air`: Density of the air (kg m-3)
# - `c_p`: Specific heat of air at constant pressure (J kg-1 K-1)
# - `ψ`: Water potential (Pa)
# - `rh`: Relative humidity (%)

"""
function wet_air(T_drybulb; 
    T_wetbulb=T_drybulb, 
    rh=0, 
    T_dew=nothing, 
    P_atmos=101325u"Pa",
    fO2=0.2095,
    fCO2=0.0004,
    fN2=0.79
)
    return wet_air(T_drybulb, T_wetbulb, rh, T_dew, P_atmos, fO2, fCO2, fN2)
end
function wet_air(T_drybulb, T_wetbulb, rh, T_dew, P_atmos, fO2, fCO2, fN2)
    c_p_H2O_vap = 1864.40u"J/K/kg"
    c_p_dry_air = 1004.84u"J/K/kg" # should be 1006?
    f_w = 1.0053 # (-) correction factor for the departure of the mixture of air and water vapour from ideal gas laws
    M_w = (1molH₂O |> u"kg")/1u"mol" # molar mass of water
    M_a = (fO2*molO₂ + fCO2*molCO₂ + fN2*molN₂)/1u"mol" # molar mass of air
    P_vap_sat = vapour_pressure(T_drybulb)
    if T_dew !== nothing
        P_vap = vapour_pressure(T_dew)
        rh = (P_vap / P_vap_sat) * 100
    else
        if rh !== nothing
            P_vap = P_vap_sat * rh / 100
        else
            δ_bulb = T_drybulb - T_wetbulb
            δ_P_vap = (0.000660 * (1 + 0.00115 * (Unitful.ustrip(T_wetbulb)-273.15)) * Unitful.ustrip(P) * Unitful.ustrip(δ_bulb))u"Pa"
            P_vap = vapour_pressure(T_wetbulb) - δ_P_vap
            rh = (P_vap / P_vap_sat) * 100
        end
    end
    r_w = ((M_w / M_a) * f_w * P_vap) / (P_atmos - f_w * P_vap)
    ρ_vap = P_vap * M_w / (0.998 * Unitful.R * T_drybulb) # 0.998 a correction factor?
    ρ_vap = Unitful.uconvert(u"kg/m^3",ρ_vap) # simplify units
    T_vir = T_drybulb * ((1.0 + r_w / (M_w / M_a)) / (1 + r_w))
    T_vinc = T_vir - T_drybulb
    ρ_air = (M_a / Unitful.R) * P_atmos / (0.999 * T_vir) # 0.999 a correction factor?
    ρ_air = Unitful.uconvert(u"kg/m^3",ρ_air) # simplify units
    c_p = (c_p_dry_air + (r_w * c_p_H2O_vap)) / (1 + r_w)
    ψ = if min(rh) <= 0
        -999u"Pa"
    else
        (4.615e+5 * Unitful.ustrip(T_drybulb) * log(rh / 100))u"Pa"
    end

    return (;P_vap, P_vap_sat, ρ_vap, r_w, T_vinc, ρ_air, c_p, ψ, rh)
end

"""
    dry_air(T_drybulb; kw...)
    dry_air(T_drybulb, P_atmos, elev, fO2, fCO2, fN2)

"""
dry_air(T_drybulb; P_atmos=nothing, elev=0m, fO2 = 0.2095, fCO2 = 0.0004, fN2 = 0.79) = dry_air(T_drybulb, P_atmos, elev, fO2, fCO2, fN2)
function dry_air(T_drybulb, P_atmos, elev, fO2, fCO2, fN2)
    σ = Unitful.uconvert(u"W/m^2/K^4",Unitful.σ) # Stefan-Boltzmann constant, W/m^2/K^4, extract σ when calling Unitful when units issue is fixed in Unitful
    M_a = ((fO2*molO₂ + fCO2*molCO₂ + fN2*molN₂) |> u"kg")/1u"mol" # molar mass of air
    if isnothing(P_atmos)
        P_atmos = get_pressure(elev)#P_std * ((1 - (0.0065 * elev / 288m))^(1 / 0.190284))
    end
    ρ_air = (M_a / Unitful.R) * P_atmos / (T_drybulb)
    ρ_air = Unitful.uconvert(u"kg/m^3",ρ_air) # simplify units
    vis_not = 1.8325e-5u"kg/m/s"
    T_not = 296.16u"K"
    c = 120u"K"
    μ = (vis_not * (T_not + c) / (T_drybulb + c)) * (T_drybulb / T_not)^1.5 # kg / m.s
    ν = μ / ρ_air # m2 / s or J.s/kg
    dif_vpr = 2.26e-5u"m^2/s" * ((T_drybulb / 273.15u"K")^1.81) * (1.e5u"Pa" / P_atmos) # m2 / s
    k_fluid = (0.02425 + (7.038e-5 * (Unitful.ustrip(T_drybulb) - 273.15)))u"W/m/K"
    L_v = (2.5012e6 - 2.3787e3 * (Unitful.ustrip(T_drybulb) - 273.15))u"J/kg"

    tcoeff = 1 / T_drybulb
    ggroup = 0.0980616u"m/s^2" * tcoeff / (ν^2) # 1 / m3.K
    bbemit = σ * ((T_drybulb)^4) # W/m2
    emtmax = 2.897e-3u"K*m" / (T_drybulb) # m

    return (;P_atmos, ρ_air, μ, ν, dif_vpr, k_fluid, L_v, tcoeff, ggroup, bbemit, emtmax)
end

function get_λ_evap(T)
    Tw = Unitful.ustrip(u"°C"(T))
    if Tw > 0
        return (2500.8 - 2.36 * Tw + 0.0016 * Tw^2 - 0.00006 * Tw^3)*1000u"J/kg"
    else
        return (834.1 - 0.29 * Tw - 0.004 * Tw^2)*1000u"J/kg"
    end
end

function waterprop(T::Quantity)
    # Ensure temperature is in °C
    T = ustrip(u"°C", T)  # Convert to Float64 in °C

    # Specific heat capacity (J/kg·K)
    c_p = (4220.02 - 4.5531 * T + 0.182958 * T^2 -
         0.00310614 * T^3 + 1.89399e-5 * T^4) * u"J/kg/K"

    # Density (kg/m^3)
    if T < 30
        ρ = 1000.0 * u"kg/m^3"
    elseif T <= 60
        ρ = (1017.0 - 0.6 * T) * u"kg/m^3"
    else
        # Clamp to 60°C
        ρ = (1017.0 - 0.6 * 60.0) * u"kg/m^3"
    end

    # Thermal conductivity (W/m·K)
    K = (0.551666 + 0.00282144 * T - 2.02383e-5 * T^2) * u"W/m/K"

    # Dynamic viscosity (kg/m·s)
    μ = (0.0017515 - 4.31502e-5 * T + 3.71431e-7 * T^2) * u"kg/m/s"

    return (
        c_p_H2O = c_p,
        ρ_H2O = ρ,
        k_H2O = K,
        μ_H2O = μ,
    )
end

function phase_transition(
    T::Vector,       # current temps at nodes
    T_past::Vector,  # temps at previous step
    ∑phase::Vector,  # accumulated latent heat
    θ::Vector,       # soil moisture by layer
    dep::Vector      # soil depth boundaries (cm)
)
    HTOFN = 333500.0u"J/kg" # latent heat of fusion of waterper unit mass
    c_p = 4186.0u"J/kg/K" # specific heat of water
    nodes = length(dep)
    layermass = zeros(Float64, nodes)u"kg"
    qphase = zeros(Float64, nodes)u"J"
    meanT = similar(T)
    meanTpast = similar(T)
    for j in 1:nodes
        if θ[j] > 0.0
            if j < nodes
                meanT[j] = (T[j] + T[j+1]) / 2.0
                meanTpast[j] = (T_past[j] + T_past[j+1]) / 2.0
            else
                meanT[j] = T[j]
                meanTpast[j] = T_past[j]
            end

            if meanTpast[j] > 273.15u"K" && meanT[j] <= 273.15u"K"
                if j < nodes
                    layermass[j] = u"m"(dep[j+1] - dep[j]) * 1000.0u"kg/m" * θ[j]
                else
                    layermass[j] = u"m"(dep[j] + 100.0u"cm" - dep[j]) * 1000.0u"kg/m" * θ[j]
                end

                qphase[j] = (meanTpast[j] - meanT[j]) * layermass[j] * c_p
                ∑phase[j] += qphase[j]

                if ∑phase[j] > HTOFN * layermass[j]
                    # Fully frozen
                    T[j] = 273.14u"K"
                    if j < nodes
                        T[j+1] = 273.14u"K"
                    end
                    ∑phase[j] = 0.0u"J"
                else
                    # In the process of freezing
                    T[j] = 273.16u"K"
                    if j < nodes
                        T[j+1] = 273.16u"K"
                    end
                end
            end
        end
    end
    return (;
        ∑phase,
        qphase,
        T
    )
end