# idealEoS.jl - Ideal Gas Equation of State

# Structure (type) definition
# ---------------------------

struct IdealGas{ℙ <: FLOAT}
    form::String            # formula
    name::String            # name
    hmod::SpecificHeat{ℙ}   # heat model
    Pref::ℙ                 # reference pressure, kPa
    function IdealGas(
            FORM::AbstractString,
            NAME::AbstractString,
            HMOD::SpecificHeat{ℙ},
            PREF::ℙ = one(ℙ)
        ) where {ℙ <: FLOAT}
        @assert(length(FORM) > 0, "Error: Empty formula")
        @assert(length(NAME) > 0, "Error: Empty name")
        @assert(PREF > 0, "Error: Pref <= 0")
        return new{ℙ}(String(FORM), String(NAME), HMOD, PREF)
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
) where {ℙ} = IdealGas(FORM, NAME, ℙ(HMOD), ℙ(PREF))

# Heat model type conversion / 2 indirections
IdealGas(
    FORM::AbstractString,
    NAME::AbstractString,
    HMOD::SpecificHeat{ℙ},
    PREF::Real = one(ℙ)
) where {ℙ} = IdealGas{ℙ}(FORM, NAME, HMOD, PREF)

# Set type with unit conversion and stripping / 2 indirections
function IdealGas{ℙ}(
        FORM::AbstractString,
        NAME::AbstractString,
        HMOD::SpecificHeat,
        PREF::Quantity{<:Real, dimension(u"kPa")} = one(ℙ) * u"kPa"
    ) where {ℙ <: FLOAT}
    return IdealGas{ℙ}(
        FORM, NAME, HMOD,
        PREF isa Quantity ? uconvert(u"kPa", Pref).val : Pref,
    )
end

# Heat model type with unit conversion and stripping / 3 indirections
IdealGas(
    FORM::AbstractString,
    NAME::AbstractString,
    HMOD::SpecificHeat{ℙ},
    PREF::Quantity{<:Real, dimension(u"kPa")}
) where {ℙ} = IdealGas{ℙ}(FORM, NAME, HMOD, PREF)

# Conversions
# -----------

import Base: convert

convert(::Type{IdealGas{ℙ}}, ξ::IdealGas{ℙ}) where {ℙ <: FLOAT} = ξ

function convert(::Type{IdealGas{ℙ}}, ξ::IdealGas{ℚ}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return IdealGas{ℙ}(ξ.form, ξ.name, ξ.hmod, ξ.Pref)
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

# User-facing functions
# ---------------------

function Base.show(io::IO, G::IdealGas{ℙ}) where {ℙ <: FLOAT}
    return print(io, "$(G.form) gas, $(G.hmod)")
end

for FUNC in (:cp, :cv, :u, :h, :s0)
    @eval begin
        ($FUNC(G::IdealGas{ℙ}, T::Real, B::Symbol)::ℙ) where {ℙ} = $FUNC(G.hmod, T, B)
    end
end

for FUNC in (:gamma,)
    @eval begin
        ($FUNC(G::IdealGas{ℙ}, T::Real)::ℙ) where {ℙ} = $FUNC(G.hmod, T)
    end
end

for FUNC in (:R,)
    @eval begin
        ($FUNC(G::IdealGas{ℙ}, B::Symbol)::ℙ) where {ℙ} = $FUNC(G.hmod, B)
    end
end

# Internal, fast, positional, EoS functions
(_P(G::IdealGas{ℙ}, T::Real, v::Real, B::Symbol)::ℙ) where {ℙ} = R(G, B) * ℙ(T / v)

(_T(G::IdealGas{ℙ}, P::Real, v::Real, B::Symbol)::ℙ) where {ℙ} = ℙ(P * v) / R(G, B)

(_v(G::IdealGas{ℙ}, P::Real, T::Real, B::Symbol)::ℙ) where {ℙ} = R(G, B) * ℙ(T / P)

(_r(G::IdealGas{ℙ}, P::Real, T::Real, B::Symbol)::ℙ) where {ℙ} = inv(_v(G, P, T, B))

# Keyworded, user-facing counterparts
(P(G::IdealGas{ℙ}; T::Real, v::Real, B::Symbol = :MA)::ℙ) where {ℙ} = _P(G, T, v, B)

(T(G::IdealGas{ℙ}; P::Real, v::Real, B::Symbol = :MA)::ℙ) where {ℙ} = _T(G, P, v, B)

(v(G::IdealGas{ℙ}; P::Real, T::Real, B::Symbol = :MA)::ℙ) where {ℙ} = _v(G, P, T, B)

(r(G::IdealGas{ℙ}; P::Real, T::Real, B::Symbol = :MA)::ℙ) where {ℙ} = _r(G, P, T, B)

export P, T, v, r

# Internal, fast, positional, entropy function
function _s(G::IdealGas{ℙ}, P::Real, T::Real, B::Symbol)::ℙ where {ℙ}
    return s0(G, T, B) - R(G, B) * log(ℙ(P) / G.Pref)
end

# Keyworded, user-facing entropy

function s(
        G::IdealGas{ℙ};
        P::Union{Missing, Real} = missing,
        T::Union{Missing, Real} = missing,
        v::Union{Missing, Real} = missing,
        B::Symbol = :MA
    ) where {ℙ}
    @assert(
        count(x -> isa(x, Real), (P, T, v)) == 2,
        "exactly two P-T-v state functions must be specified!"
    )
    return if ismissing(P)
        _s(G, _P(G, T, v, B), T, B)
    elseif ismissing(T)
        _s(G, P, _T(G, P, v, B), B)
    else
        _s(G, P, T, B)
    end
end

export s
