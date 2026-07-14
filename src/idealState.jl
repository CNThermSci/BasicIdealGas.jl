# idealState.jl - Ideal Gas States

# Structure (type) definition
# ---------------------------

struct IdealState{ℙ <: FLOAT}
    𝐺::IdealGas{ℙ}
    𝑃::ℙ
    𝑇::ℙ
    # Internal, validating constructors
    function IdealState(
            G::IdealGas{ℙ},
            P::ℙ,
            T::ℙ,
        ) where {ℙ <: FLOAT}
        return new{ℙ}(G, P, T)
    end
end

# Base.getproperty
# ----------------

import Base: getproperty, propertynames

function Base.getproperty(st::IdealState, sy::Symbol)
    if sy in (:𝐺, :𝑃, :𝑇)
        return getfield(st, sy)
    elseif sy == :gas
        return getfield(st, :𝐺)
    elseif sy == :M
        return getfield(st, :𝐺).hmod.𝑀 * u"kg/kmol"
    elseif sy in (:RMO, :RU, :Ru)
        return getfield(st, :𝐺).hmod.𝑅 * u"kJ/kmol/K"
    elseif sy == :P
        return getfield(st, :𝑃) * u"kPa"
    elseif sy == :T
        return getfield(st, :𝑇) * u"K"
    end
    GAS, P, T = map(sy -> getfield(st, sy), (:𝐺, :𝑃, :𝑇))
    if sy == :R
        return R(GAS, :MA) * u"kJ/kg/K"
    elseif sy in (:γ, :ga)
        return ga(GAS, T)
    elseif sy == :v
        return _v(GAS, P, T, :MA) * u"m^3/kg"
    elseif sy == :vMO
        return _v(GAS, P, T, :MO) * u"m^3/kmol"
    elseif sy == :ρ
        return _ρ(GAS, P, T, :MA) * u"kg/m^3"
    elseif sy == :ρMO
        return _ρ(GAS, P, T, :MO) * u"kmol/m^3"
    elseif sy == :cp
        return cp(GAS, T, :MA) * u"kJ/kg/K"
    elseif sy == :cpMO
        return cp(GAS, T, :MO) * u"kJ/kmol/K"
    elseif sy == :cv
        return cv(GAS, T, :MA) * u"kJ/kg/K"
    elseif sy == :cvMO
        return cv(GAS, T, :MO) * u"kJ/kmol/K"
    elseif sy == :u
        return u(GAS, T, :MA) * u"kJ/kg"
    elseif sy == :uMO
        return u(GAS, T, :MO) * u"kJ/kmol"
    elseif sy == :h
        return h(GAS, T, :MA) * u"kJ/kg"
    elseif sy == :hMO
        return h(GAS, T, :MO) * u"kJ/kmol"
    elseif sy == :s0
        return s0(GAS, T, :MA) * u"kJ/kg/K"
    elseif sy == :s0MO
        return s0(GAS, T, :MO) * u"kJ/kmol/K"
    end
end

Base.propertynames(::IdealState) = (
    :gas, :M, :R, :Ru, :RU, :P, :T, :γ, :ga,
    :v, :vMO, :ρ, :ρMO, :cp, :cpMO, :cv, :cvMO,
    :u, :uMO, :h, :hMO, :s0, :s0MO,
)

# Export
# ------

export IdealState

# User-facing functions
# ---------------------

function Base.show(io::IO, st::IdealState{ℙ}) where {ℙ <: FLOAT}
    return print(io, "$(st.gas) at ($(st.P), $(st.T))")
end
