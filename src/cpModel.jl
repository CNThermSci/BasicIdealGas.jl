# cpModel.jl - General Specific Heat Model

const Ru = 8.31447  # Universal R, kJ/kmol·K

struct SpecificHeat
    id::Symbol
    cp_f::Function
    M::Float64
    Tmin::Float64
    Tmax::Float64
    Tref::Float64
    uref::Float64
    sref::Float64
    RU::Float64
    SpecificHeat(
        ID::Symbol, CP_F::Function, M::Real,
        Tmin::Real, Tmax::Real, Tref::Real,
        uref::Real, sref::Real, B::Symbol,
        RU::Real = Ru
    ) = begin
        @assert ID != Symbol("")
        @assert RU > 0
        @assert M > 0
        @assert 0 <= Tmin <= Tref < Tmax
        @assert B in (:MA, :MO)
        mult = B == :MA ? M : 1.0
        new(
            ID, T -> Float64(CP_F(T) * mult), Float64(M),
            Float64(Tmin), Float64(Tmax), Float64(Tref),
            Float64(uref * mult), Float64(sref * mult), RU
        )
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
