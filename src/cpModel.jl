# cpModel.jl - General Specific Heat Model

# Type aliasing
# -------------

# IEEE-754 normalized floating point types of half, single, and double precision
FLOAT = Base.IEEEFloat

# Precision Composition Simplification
⊚(p::Type{ℙ}, f::Function) where ℙ <: FLOAT = f(1) isa ℙ ? f : p ∘ f

# Chained Precision Composition Simplification
⊚(p::Type{ℙ}, c::ComposedFunction{Type{ℚ}}) where {ℙ <: FLOAT, ℚ <: FLOAT} = ⊚(p, c.inner)

# Structure (type) definition
# ---------------------------

struct SpecificHeat{ℙ <: FLOAT}
    ID::Symbol
    FN::Function
    M::ℙ        # kg/kmol
    Tmin::ℙ     # K
    Tref::ℙ     # K
    Tmax::ℙ     # K
    uref::ℙ     # kJ/kmol
    sref::ℙ     # kJ/kmol⋅K
    RU::ℙ       # kJ/kmol⋅K
    # Internal constructors
    # Validating
    SpecificHeat(
        ID::Symbol, FN::Function, M::ℙ,
        Tmin::ℙ, Tref::ℙ, Tmax::ℙ,
        uref::ℙ, sref::ℙ, RU::ℙ,
        B::Symbol
    ) where {ℙ <: FLOAT} = begin
        @assert(ID != Symbol(""), "Error: Empty ID")
        @assert(M > zero(ℙ), "Error: M <= 0")
        @assert(zero(ℙ) <= Tmin <= Tref < Tmax, "Error: Temperature values")
        @assert(RU > zero(ℙ), "Error: RU <= 0")
        @assert(B in (:MA, :MO), "Error: B should be either :MA or :MO")
        return B == :MA ? (
            new{ℙ}(ID, ℙ ⊚ T -> FN(T) * M, M, Tmin, Tref, Tmax, uref * M, sref * M, RU)
        ) : (
            new{ℙ}(ID, ℙ ⊚ FN, M, Tmin, Tref, Tmax, uref, sref, RU)
        )
    end
end

# External constructors
# ---------------------

# Set type conversion / 1 indirection
function SpecificHeat{ℙ}(
        ID::Symbol, FN::Function, M::Real,
        Tmin::Real, Tref::Real, Tmax::Real,
        uref::Real, sref::Real, RU::Real,
        B::Symbol
    ) where {ℙ <: FLOAT}
    return SpecificHeat(ID, ℙ ⊚ FN, ℙ.((M, Tmin, Tref, Tmax, uref, sref, RU))..., B)
end

# Promotion type conversion / 2 indirections
function SpecificHeat(
        ID::Symbol, FN::Function, M::Real,
        Tmin::Real, Tref::Real, Tmax::Real,
        uref::Real, sref::Real, RU::Real,
        B::Symbol
    )
    ℙ = promote_type(typeof.((M, Tmin, Tref, Tmax, uref, sref, RU))...)
    ℙ = ℙ <: FLOAT ? ℙ : Float64
    return SpecificHeat{ℙ}(ID, FN, M, Tmin, Tref, Tmax, uref, sref, RU, B)
end

# Set type with unit conversion and stripping / 2 indirections
function SpecificHeat{ℙ}(
        ID::Symbol,
        FN::Function,
        M::Union{Real, Quantity{<:Real, dimension(u"kg/kmol")}},
        Tmin::Union{Real, Quantity{<:Real, dimension(u"K")}},
        Tref::Union{Real, Quantity{<:Real, dimension(u"K")}},
        Tmax::Union{Real, Quantity{<:Real, dimension(u"K")}},
        uref::Union{
            Quantity{<:Real, dimension(u"kJ/kmol")},
            Quantity{<:Real, dimension(u"kJ/kg")},
        },
        sref::Union{
            Quantity{<:Real, dimension(u"kJ/kmol/K")},
            Quantity{<:Real, dimension(u"kJ/kg/K")},
        },
        RU::Union{Real, Quantity{<:Real, dimension(u"kJ/kmol/K")}}
    ) where {ℙ <: FLOAT}
    _M = M isa Quantity ? uconvert(u"kg/kmol", M).val : M
    _uMO = uref isa Quantity{<:Real, dimension(u"kJ/kmol")} ? (
        uconvert(u"kJ/kmol", uref).val) : (
        uconvert(u"kJ/kg", uref).val * _M)
    _sMO = sref isa Quantity{<:Real, dimension(u"kJ/kmol/K")} ? (
        uconvert(u"kJ/kmol/K", sref).val) : (
        uconvert(u"kJ/kg/K", sref).val * _M)
    return SpecificHeat{ℙ}(
        ID, FN, _M,
        Tmin isa Quantity ? uconvert(u"K", Tmin).val : Tmin,
        Tref isa Quantity ? uconvert(u"K", Tref).val : Tref,
        Tmax isa Quantity ? uconvert(u"K", Tmax).val : Tmax,
        _uMO, _sMO,
        RU isa Quantity ? uconvert(u"kJ/kmol/K", RU).val : RU, :MO,
    )
end

# Promotion type with unit conversion and stripping / 3 indirections
function SpecificHeat(
        ID::Symbol,
        FN::Function,
        M::Union{𝕄, Quantity{<:𝕄, dimension(u"kg/kmol")}},
        Tmin::Union{𝕀, Quantity{<:𝕀, dimension(u"K")}},
        Tref::Union{𝔸, Quantity{<:𝔸, dimension(u"K")}},
        Tmax::Union{𝔼, Quantity{<:𝔼, dimension(u"K")}},
        uref::Union{
            Quantity{<:𝕌, dimension(u"kJ/kmol")},
            Quantity{<:𝕌, dimension(u"kJ/kg")},
        },
        sref::Union{
            Quantity{<:𝕊, dimension(u"kJ/kmol/K")},
            Quantity{<:𝕊, dimension(u"kJ/kg/K")},
        },
        RU::Union{ℝ, Quantity{<:ℝ, dimension(u"kJ/kmol/K")}}
    ) where {𝕄 <: Real, 𝕀 <: Real, 𝔸 <: Real, 𝔼 <: Real, 𝕌 <: Real, 𝕊 <: Real, ℝ <: Real}
    ℙ = promote_type(𝕄, 𝕀, 𝔸, 𝔼, 𝕌, 𝕊, ℝ)
    ℙ = ℙ <: FLOAT ? ℙ : Float64
    return SpecificHeat{ℙ}(ID, FN, M, Tmin, Tref, Tmax, uref, sref, RU)
end

# Conversions
# -----------

import Base: convert

convert(::Type{SpecificHeat{ℙ}}, ξ::SpecificHeat{ℙ}) where {ℙ <: FLOAT} = ξ

function convert(::Type{SpecificHeat{ℙ}}, ξ::SpecificHeat{ℚ}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return SpecificHeat{ℙ}(
        ξ.ID, ξ.FN, ξ.M, ξ.Tmin, ξ.Tref, ξ.Tmax, ξ.uref, ξ.sref, ξ.RU, :MO
    )
end

import Base: Float16, Float32, Float64

Float16(ξ::SpecificHeat) = convert(SpecificHeat{Float16}, ξ)
Float32(ξ::SpecificHeat) = convert(SpecificHeat{Float32}, ξ)
Float64(ξ::SpecificHeat) = convert(SpecificHeat{Float64}, ξ)

# Promotions
# ----------

import Base: promote_rule

function promote_rule(::Type{SpecificHeat{ℙ}}, ::Type{SpecificHeat{ℚ}}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return SpecificHeat{promote_type(ℙ, ℚ)}
end

# Export
# ------

export SpecificHeat

# Show
# ----

# Auxiliary methods
function subscript(x::Int)
    asSub(c::Char) = Char(Int(c) - Int('0') + Int('₀'))
    map(asSub, "$(x)")
end

pDeco(::Type{Float16}) = subscript(16)
pDeco(::Type{Float32}) = subscript(32)
pDeco(::Type{Float64}) = subscript(64)

function Base.show(io::IO, S::SpecificHeat{ℙ}) where {ℙ <: FLOAT}
    return print(io, "$(S.ID) cp$(pDeco(ℙ))(T) model, $(S.Tmin) <= T <= $(S.Tmax)")
end

# User-facing functions
# ---------------------

import Base: cp

function cp(C::SpecificHeat{ℙ}, T::Real, B::Symbol)::ℙ where {ℙ <: FLOAT}
    @assert B in (:MA, :MO)
    @assert C.Tmin <= T <= C.Tmax
    divisor = B == :MA ? C.M : one(ℙ)
    return C.FN(T) / divisor
end

function R(C::SpecificHeat, B::Symbol)
    @assert B in (:MA, :MO)
    divisor = B == :MA ? C.M : 1.0
    return C.RU / divisor
end

function cv(C::SpecificHeat, T::Real, B::Symbol)
    return cp(C, T, B) - R(C, B)
end

gamma(C::SpecificHeat, T::Real) = cp(C, T, :MO) / cv(C, T, :MO)

function u(C::SpecificHeat, T::Real, B::Symbol)
    T = Float64(T)
    @assert B in (:MA, :MO)
    divisor = B == :MA ? C.M : 1.0
    IE = quadgk(T -> cv(C, T, :MO), C.Tref, T, rtol = 1.0e-5)
    return (IE[1] + C.uref) / divisor
end

h(C::SpecificHeat, T::Real, B::Symbol) = u(C, T, B) + R(C, B) * Float64(T)

function s0(C::SpecificHeat, T::Real, B::Symbol)
    T = Float64(T)
    @assert B in (:MA, :MO)
    divisor = B == :MA ? C.M : 1.0
    IE = quadgk(T -> cp(C, T, :MO) / T, C.Tref, T, rtol = 1.0e-5)
    return (IE[1] + C.sref) / divisor
end

export cp, cv, R, gamma, u, h, s0
