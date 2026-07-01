# cpModel.jl - General Specific Heat Model

# Type aliasing
# -------------

# IEEE-754 normalized floating point types of half, single, and double precision
FLOAT = Base.IEEEFloat


const Ru = 8.31447  # Universal R, kJ/kmol·K

# Structure (type) definition
# ---------------------------

struct SpecificHeat{ℙ <: FLOAT}
    id::Symbol
    cp_f::Function
    M::ℙ        # kg/kmol
    Tmin::ℙ     # K
    Tmax::ℙ     # K
    Tref::ℙ     # K
    uref::ℙ     # kJ/kmol
    sref::ℙ     # kJ/kmol⋅K
    RU::ℙ       # kJ/kmol⋅K
    # Internal constructors
    # Validating
    SpecificHeat(
        ID::Symbol, CP_F::Function, M::ℙ,
        Tmin::ℙ, Tmax::ℙ, Tref::ℙ,
        uref::ℙ, sref::ℙ, B::Symbol,
        RU::ℙ
    ) where {ℙ <: FLOAT} = begin
        @assert(ID != Symbol(""), "Error: Empty ID")
        @assert(RU > zero(ℙ), "Error: RU <= 0")
        @assert(M > zero(ℙ), "Error: M <= 0")
        @assert(zero(ℙ) <= Tmin <= Tref < Tmax, "Error: Temperature values")
        @assert(B in (:MA, :MO), "Error: B should be either :MA or :MO")
        mult = B == :MA ? M : one(ℙ)
        new{ℙ}(ID, T -> ℙ(CP_F(T) * mult), M, Tmin, Tmax, Tref, uref * mult, sref * mult, RU)
    end
end

# External constructors
# ---------------------

# Set type conversion / 1 indirection
function SpecificHeat{ℙ}(
        ID::Symbol, CP_F::Function, M::Real,
        Tmin::Real, Tmax::Real, Tref::Real,
        uref::Real, sref::Real, B::Symbol,
        RU::Real = ℙ(Ru),
    ) where {ℙ <: FLOAT}
    return SpecificHeat(ID, CP_F, P.((M, Tmin, Tmax, Tref, uref, sref))..., B, ℙ(RU))
end

# Promotion type conversion / 2 indirections
function SpecificHeat(
        ID::Symbol, CP_F::Function, M::Real,
        Tmin::Real, Tmax::Real, Tref::Real,
        uref::Real, sref::Real, B::Symbol,
        RU::Real = Ru
    )
    ℙ = promote_type(typeof.((M, Tmin, Tmax, Tref, uref, sref))...) # RU left out of promotion
    ℙ = ℙ <: FLOAT ? ℙ : Float64
    return SpecificHeat{ℙ}(ID, CP_F, M, Tmin, Tmax, Tref, uref, sref, B, ℙ(RU))
end

# Set type with unit conversion and stripping / 2 indirections
function SpecificHeat{ℙ}(
        ID::Symbol,
        CP_F::Function,
        M::Union{Real, Quantity{<:Real, dimension(u"kg/kmol")}},
        Tmin::Union{Real, Quantity{<:Real, dimension(u"K")}},
        Tmax::Union{Real, Quantity{<:Real, dimension(u"K")}},
        Tref::Union{Real, Quantity{<:Real, dimension(u"K")}},
        uref::Union{
            Quantity{<:Real, dimension(u"kJ/kmol")},
            Quantity{<:Real, dimension(u"kJ/kg")},
        },
        sref::Union{
            Quantity{<:Real, dimension(u"kJ/kmol/K")},
            Quantity{<:Real, dimension(u"kJ/kg/K")},
        },
        RU::Union{Real, Quantity{<:Real, dimension(u"kJ/kmol/K")}} = ℙ(Ru) * u"kJ/kmol/K",
    ) where {ℙ <: FLOAT}
    _M = M isa Quantity ? uconvert(u"kg/kmol", M).val : M
    _uMO = uref isa Quantity{<:Real, dimension(u"kJ/kmol")} ? (
        uconvert(u"kJ/kmol", uref).val) : (
        uconvert(u"kJ/kg", uref).val * _M)
    _sMO = sref isa Quantity{<:Real, dimension(u"kJ/kmol/K")} ? (
        uconvert(u"kJ/kmol/K", sref).val) : (
        uconvert(u"kJ/kg/K", sref).val * _M)
    return SpecificHeat{ℙ}(
        ID, CP_F, _M,
        Tmin isa Quantity ? uconvert(u"K", Tmin).val : Tmin,
        Tmax isa Quantity ? uconvert(u"K", Tmax).val : Tmax,
        Tref isa Quantity ? uconvert(u"K", Tref).val : Tref,
        _uMO, _sMO, :MO,
        RU isa Quantity ? uconvert(u"kJ/kmol/K", RU).val : RU,
    )
end

# Promotion type with unit conversion and stripping / 3 indirections
function SpecificHeat(
        ID::Symbol,
        CP_F::Function,
        M::Union{𝕄, Quantity{<:𝕄, dimension(u"kg/kmol")}},
        Tmin::Union{𝕀, Quantity{<:𝕀, dimension(u"K")}},
        Tmax::Union{𝔸, Quantity{<:𝔸, dimension(u"K")}},
        Tref::Union{𝔼, Quantity{<:𝔼, dimension(u"K")}},
        uref::Union{
            Quantity{<:𝕌, dimension(u"kJ/kmol")},
            Quantity{<:𝕌, dimension(u"kJ/kg")},
        },
        sref::Union{
            Quantity{<:𝕊, dimension(u"kJ/kmol/K")},
            Quantity{<:𝕊, dimension(u"kJ/kg/K")},
        },
        RU::Union{ℝ, Quantity{<:ℝ, dimension(u"kJ/kmol/K")}} = ℙ(Ru) * u"kJ/kmol/K",
    ) where {𝕄 <: Real, 𝕀 <: Real, 𝔸 <: Real, 𝔼 <: Real, 𝕌 <: Real, 𝕊 <: Real, ℝ <: Real}
    ℙ = promote_type(𝕄, 𝕀, 𝔸, 𝔼, 𝕌, 𝕊, ℝ)
    ℙ = ℙ <: FLOAT ? ℙ : Float64
    return SpecificHeat{ℙ}(ID, CP_F, M, Tmin, Tmax, Tref, uref, sref, RU)
end

# Export
# ------

export SpecificHeat

# Show
# ----

function Base.show(io::IO, S::SpecificHeat)
    return print(io, "$(S.id) cp(T) model, $(S.Tmin) <= T <= $(S.Tmax)")
end

# User-facing functions
# ---------------------

import Base: cp

function cp(C::SpecificHeat, T::Real, B::Symbol)
    T = Float64(T)
    @assert B in (:MA, :MO)
    @assert C.Tmin <= T <= C.Tmax
    divisor = B == :MA ? C.M : 1.0
    return C.cp_f(T) / divisor
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
