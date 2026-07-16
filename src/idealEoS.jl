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
    return print(io, "$(G.form)$(pDeco(ℙ)) gas")
end

for FUNC in (:R,)
    @eval begin
        $FUNC(G::IdealGas, B::Symbol = :MA) = $FUNC(G.hmod, B)
    end
end

for FUNC in (:cp┆R, :cv┆R, :ga, :∫cp┆R, :∫cv┆R, :u┆R, :h┆R, :∫cp┆RT, :s0┆R, :Pr, :vr)
    @eval begin
        $FUNC(G::IdealGas, T::Real) = $FUNC(G.hmod, T)
    end
end

for FUNC in (:cp, :cv, :u, :h, :s0)
    @eval begin
        $FUNC(G::IdealGas, T::Real, B::Symbol = :MA) = $FUNC(G.hmod, T, B)
    end
end

# Internal, fast, positional, EoS functions
_P(G::IdealGas{ℙ}, T::Real, v::Real, B::Symbol = :MA) where {ℙ} = R(G, B) * ℙ(T / v)
_T(G::IdealGas{ℙ}, P::Real, v::Real, B::Symbol = :MA) where {ℙ} = ℙ(P * v) / R(G, B)
_v(G::IdealGas{ℙ}, P::Real, T::Real, B::Symbol = :MA) where {ℙ} = R(G, B) * ℙ(T / P)
_ρ(G::IdealGas{ℙ}, P::Real, T::Real, B::Symbol = :MA) where {ℙ} = inv(_v(G, P, T, B))

# NamedTuple arg, return type inferrable, user-facing EoS functions
function P(
        G::IdealGas,
        ζ::@NamedTuple{T::ℚ, v::ℝ, B::Symbol} where {ℚ <: Real, ℝ <: Real}
    )
    return _P(G, ζ.T, ζ.v, B)
end

function T(
        G::IdealGas,
        ζ::@NamedTuple{P::ℚ, v::ℝ, B::Symbol} where {ℚ <: Real, ℝ <: Real}
    )
    return _T(G, ζ.P, ζ.v, B)
end

function v(
        G::IdealGas,
        ζ::@NamedTuple{P::ℚ, T::ℝ, B::Symbol} where {ℚ <: Real, ℝ <: Real}
    )
    return _v(G, ζ.P, ζ.T, B)
end

function ρ( # "ρ" can be typed by \rho<tab>
        G::IdealGas,
        ζ::@NamedTuple{P::ℚ, T::ℝ, B::Symbol} where {ℚ <: Real, ℝ <: Real}
    )
    return _ρ(G, ζ.P, ζ.T, B)
end

# Keyworded user-facing EoS functions
P(G::IdealGas; T::Real, v::Real, B::Symbol = :MA) = _P(G, T, v, B)
T(G::IdealGas; P::Real, v::Real, B::Symbol = :MA) = _T(G, P, v, B)
v(G::IdealGas; P::Real, T::Real, B::Symbol = :MA) = _v(G, P, T, B)
ρ(G::IdealGas; P::Real, T::Real, B::Symbol = :MA) = _ρ(G, P, T, B)

# Internal, fast, positional, entropy function
function _s(G::IdealGas{ℙ}, P::Real, T::Real, B::Symbol = :MA)::ℙ where {ℙ}
    return s0(G, T, B) - R(G, B) * log(ℙ(P) / G.Pref)
end

# NamedTuple arg, return type inferrable, user-facing entropy functions
function s(
        G::IdealGas,
        ζ::@NamedTuple{P::ℚ, T::ℝ, B::Symbol} where {ℚ <: Real, ℝ <: Real}
    )
    return _s(G, ζ.P, ζ.T, B)
end

function s(
        G::IdealGas,
        ζ::@NamedTuple{P::ℚ, v::ℝ, B::Symbol} where {ℚ <: Real, ℝ <: Real}
    )
    return _s(G, ζ.P, _T(G, ζ.P, ζ.v, B), B)
end

function s(
        G::IdealGas,
        ζ::@NamedTuple{T::ℚ, v::ℝ, B::Symbol} where {ℚ <: Real, ℝ <: Real}
    )
    return _s(G, _P(G, ζ.T, ζ.v, B), ζ.T, B)
end

# Keyworded, user-facing entropy functions
function s(
        G::IdealGas;
        P::Union{Missing, Real} = missing,
        T::Union{Missing, Real} = missing,
        v::Union{Missing, Real} = missing,
        B::Symbol = :MA
    )
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

# User-facing exports
export P, T, v, ρ, s

# Base.getproperty
# ----------------

import Base: getproperty, propertynames

function Base.getproperty(ξ::IdealGas, sy::Symbol)
    # Raw fields
    if sy in fieldnames(IdealGas) return getfield(ξ, sy) end
    # Short-circuit SpecificHeat model accessors
    if sy in propertynames(getfield(ξ, :hmod))
        return getproperty(getfield(ξ, :hmod), sy)
    end
end

Base.propertynames(ξ::IdealGas) = (
    :form, :name, :hmod, :Pref,
    propertynames(getfield(ξ, :hmod))...
)
