# idealState.jl - Ideal Gas States

# Structure (type) definition
# ---------------------------

struct IdealState{ℙ <: FLOAT}
    𝐺::IdealGas{ℙ}
    𝑃::ℙ                # kPa
    𝑇::ℙ                # K
    # Internal, validating constructors
    function IdealState(
            G::IdealGas{ℙ},
            P::ℙ,
            T::ℙ,
        ) where {ℙ <: FLOAT}
        @assert(G.hmod.Tmin <= T <= G.hmod.Tmax, "T = $(T) out of range for $(G.hmod)")
        return new{ℙ}(G, P, T)
    end
end

# External constructors
# ---------------------

# Set type conversion / 1 indirection
IdealState{ℙ}(
    G::IdealGas,
    P::Real,
    T::Real,
) where {ℙ} = IdealState(ℙ.((G, P, T))...)

# Ideal gas model type conversion / 2 indirections
IdealState(
    G::IdealGas{ℙ},
    P::Real,
    T::Real,
) where {ℙ} = IdealState{ℙ}(G, P, T)

# Set type with unit conversion and stripping / 2 indirections
function IdealState{ℙ}(
        G::IdealGas,
        P::Quantity{<:Real, dimension(u"kPa")},
        T::Quantity{<:Real, dimension(u"K")},
    ) where {ℙ <: FLOAT}
    return IdealState{ℙ}(
        G,
        P isa Quantity ? uconvert(u"kPa", P).val : P,
        T isa Quantity ? uconvert(u"K", T).val : T,
    )
end

# Heat model type with unit conversion and stripping / 3 indirections
IdealState(
    G::IdealGas{ℙ},
    P::Quantity{<:Real, dimension(u"kPa")},
    T::Quantity{<:Real, dimension(u"K")},
) where {ℙ <: FLOAT} = IdealState{ℙ}(G, P, T)

# Conversions
# -----------

import Base: convert

convert(::Type{IdealState{ℙ}}, ξ::IdealState{ℙ}) where {ℙ <: FLOAT} = ξ

function convert(::Type{IdealState{ℙ}}, ξ::IdealState{ℚ}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return IdealState{ℙ}(ξ.𝐺, ξ.𝑃, ξ.𝑇)
end

import Base: Float16, Float32, Float64

Float16(ξ::IdealState) = convert(IdealState{Float16}, ξ)
Float32(ξ::IdealState) = convert(IdealState{Float32}, ξ)
Float64(ξ::IdealState) = convert(IdealState{Float64}, ξ)

# Promotions
# ----------

import Base: promote_rule

function promote_rule(::Type{IdealState{ℙ}}, ::Type{IdealState{ℚ}}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return IdealState{promote_type(ℙ, ℚ)}
end

# Export
# ------

export IdealState

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
    elseif sy == :Pr
        return Pr(GAS, T)
    elseif sy == :vr
        return vr(GAS, T)
    end
end

Base.propertynames(::IdealState) = (
    :gas, :M, :R, :Ru, :RU, :P, :T, :γ, :ga,
    :v, :vMO, :ρ, :ρMO, :cp, :cpMO, :cv, :cvMO,
    :u, :uMO, :h, :hMO, :s0, :s0MO, :Pr, :vr,
)

# User-facing functions
# ---------------------

function Base.show(io::IO, st::IdealState{ℙ}) where {ℙ <: FLOAT}
    return print(io, "$(st.gas) @($(st.P), $(st.T))")
end
