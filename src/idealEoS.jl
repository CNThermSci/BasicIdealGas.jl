# idealEoS.jl - Ideal Gas Equation of State

# Structure (type) definition
# ---------------------------

struct IdealGas{ℙ <: FLOAT}
    form::String            # formula
    name::String            # name
    hmod::SpecificHeat{ℙ}   # heat model
    Pref::ℙ                 # reference pressure
    function IdealGas(
            FORM::AbstractString,
            NAME::AbstractString,
            HMOD::SpecificHeat{ℙ},
            PREF::ℙ = one(ℙ)
        ) where {ℙ <: FLOAT}
        @assert("Error: Pref <= 0", PREF > 0)
        @assert("Error: Empty formula", length(FORM) > 0)
        @assert("Error: Empty name", length(NAME) > 0)
        new{ℙ}(String(FORM), String(NAME), HMOD, PREF)
    end
end

# External constructors
# ---------------------

# Set type conversion / 1 indirection
IdealGas{ℙ}(
    FORM::AbstractString,
    NAME::AbstractString,
    HMOD::SpecificHeat,
    PREF::Real = one(ℙ)
) = IdealGas(FORM, NAME, ℙ(HMOD), ℙ(PREF))

# Conversions
# -----------

import Base: convert

convert(::Type{IdealGas{ℙ}}, ξ::IdealGas{ℙ}) where {ℙ <: FLOAT} = ξ

function convert(::Type{IdealGas{ℙ}}, ξ::IdealGas{ℚ}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return IdealGas{ℙ}(
        ξ.FN, ξ.MM, ξ.Tmin, ξ.Tref, ξ.Tmax, ξ.uref, ξ.sref, ξ.RU, :MO
    )
end

import Base: Float16, Float32, Float64

Float16(ξ::IdealGas) = convert(IdealGas{Float16}, ξ)
Float32(ξ::IdealGas) = convert(IdealGas{Float32}, ξ)
Float64(ξ::IdealGas) = convert(IdealGas{Float64}, ξ)

# Promotions
# ----------

import Base: promote_rule

function promote_rule(::Type{IdealGas{ℙ}}, ::Type{IdealGas{ℚ}}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return IdealGas{promote_type(ℙ, ℚ)}
end

# Export
# ------

export IdealGas

function Base.show(io::IO, G::IdealGas)
    print(io, "$(G.form) gas with $(G.hmod)")
end

for FUNC in (:cp, :cv, :u, :h, :s0)
    @eval begin
        $FUNC(G::IdealGas, T::Real, B::Symbol) = $FUNC(G.hmod, T, B)
    end
end

for FUNC in (:gamma, )
    @eval begin
        $FUNC(G::IdealGas, T::Real) = $FUNC(G.hmod, T)
    end
end

for FUNC in (:R, )
    @eval begin
        $FUNC(G::IdealGas, B::Symbol) = $FUNC(G.hmod, B)
    end
end

# Internal, fast, positional, EoS functions

_P(G::IdealGas, T::Real, v::Real, B::Symbol) = R(G, B) * Float64(T) / Float64(v)

_T(G::IdealGas, P::Real, v::Real, B::Symbol) = Float64(P) * Float64(v) / R(G, B)

_v(G::IdealGas, P::Real, T::Real, B::Symbol) = R(G, B) * Float64(T) / Float64(P)

_r(G::IdealGas, P::Real, T::Real, B::Symbol) = inv(_v(G, P, T, B))

# Keyworded, user-facing counterparts

P(G::IdealGas; T::Real, v::Real, B::Symbol = :MA) = _P(G, T, v, B)

T(G::IdealGas; P::Real, v::Real, B::Symbol = :MA) = _T(G, P, v, B)

v(G::IdealGas; P::Real, T::Real, B::Symbol = :MA) = _v(G, P, T, B)

r(G::IdealGas; P::Real, T::Real, B::Symbol = :MA) = _r(G, P, T, B)

export P, T, v, r

# Internal, fast, positional, entropy function

function _s(G::IdealGas, P::Real, T::Real, B::Symbol)
    return s0(G, T, B) - R(G, B) * log(Float64(P) / G.Pref)
end

# Keyworded, user-facing entropy

function s(G::IdealGas;
           P::Union{Missing,Real} = missing,
           T::Union{Missing,Real} = missing,
           v::Union{Missing,Real} = missing,
           B::Symbol = :MA)
    @assert(count(x -> isa(x, Real), (P, T, v)) == 2,
        "exactly two P-T-v state functions must be specified!")
    return if ismissing(P)
        _s(G, _P(G, T, v, B), T, B)
    elseif ismissing(T)
        _s(G, P, _T(G, P, v, B), B)
    else
        _s(G, P, T, B)
    end
end

export s
