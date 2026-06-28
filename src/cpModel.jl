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
    M::ℙ
    Tmin::ℙ
    Tmax::ℙ
    Tref::ℙ
    uref::ℙ
    sref::ℙ
    RU::ℙ
    SpecificHeat(
        ID::Symbol, CP_F::Function, M::ℙ,
        Tmin::ℙ, Tmax::ℙ, Tref::ℙ,
        uref::ℙ, sref::ℙ, B::Symbol,
        RU::ℙ
    ) where {ℙ <: FLOAT} = begin
        @assert ID != Symbol("")
        @assert RU > zero(ℙ)
        @assert M > zero(ℙ)
        @assert zero(ℙ) <= Tmin <= Tref < Tmax
        @assert B in (:MA, :MO)
        mult = B == :MA ? M : one(ℙ)
        new{ℙ}(ID, T -> ℙ(CP_F(T) * mult), M, Tmin, Tmax, Tref, uref * mult, sref * mult, RU)
    end
end

export SpecificHeat

function Base.show(io::IO, S::SpecificHeat)
    return print(io, "$(S.id) cp(T) model, $(S.Tmin) <= T <= $(S.Tmax)")
end

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
