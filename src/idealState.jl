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

function Base.getproperty(st::IdealState, s::Symbol)
    if s in (:𝐺, :𝑃, :𝑇)
        return getfield(st, s)
    elseif s == :gas
        return getfield(st, :𝐺)
    elseif s == :P
        return getfield(st, :𝑃) * u"kPa"
    elseif s == :T
        return getfield(st, :𝑇) * u"K"
    end
    GAS, P, T = map(fi -> getfield(st, fi), (:𝐺, :𝑃, :𝑇))
    if s == :v
        return _v(GAS, P, T, :MA) * u"m^3/kg"
    elseif s == :vMO
        return _v(GAS, P, T, :MO) * u"m^3/kmol"
    elseif s == :ρ
        return _ρ(GAS, P, T, :MA) * u"kg/m^3"
    elseif s == :ρMO
        return _ρ(GAS, P, T, :MO) * u"kmol/m^3"
    elseif s == :cp
        return cp(GAS, T, :MA) * u"kJ/kg/K"
    elseif s == :cpMO
        return cp(GAS, T, :MO) * u"kJ/kmol/K"
    elseif s == :cv
        return cv(GAS, T, :MA) * u"kJ/kg/K"
    elseif s == :cvMO
        return cv(GAS, T, :MO) * u"kJ/kmol/K"
    elseif s == :u
        return u(GAS, T, :MA) * u"kJ/kg"
    elseif s == :uMO
        return u(GAS, T, :MO) * u"kJ/kmol"
    elseif s == :h
        return h(GAS, T, :MA) * u"kJ/kg"
    elseif s == :hMO
        return h(GAS, T, :MO) * u"kJ/kmol"
    end
end

Base.propertynames(::IdealState) = (
    :gas, :P, :T,
    :v, :vMO, :ρ, :ρMO, :cp, :cpMO, :cv, :cvMO,
    :u, :uMO, :h, :hMO, :s0, 
)

# Export
# ------

export IdealState

# User-facing functions
# ---------------------

function Base.show(io::IO, st::IdealState{ℙ}) where {ℙ <: FLOAT}
    return print(io, "$(st.gas) at ($(st.P), $(st.T))")
end
