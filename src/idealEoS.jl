# idealEoS.jl - Ideal Gas Equation of State

# Structure (type) definition
# ---------------------------

struct IdealGas{ℙ <: FLOAT}
    form::String
    name::String
    cpMod::SpecificHeat
    Pref::Float64
    function IdealGas(FORM::String, NAME::String,
                      CPMD::SpecificHeat, PREF::Real = 1.0)
        Pref = Float64(PREF)
        @assert Pref > 0.0
        @assert FORM != ""
        @assert NAME != ""
        new(FORM, NAME, CPMD, Pref)
    end
end

export IdealGas

function Base.show(io::IO, G::IdealGas)
    print(io, "$(G.form) gas with $(G.cpMod)")
end

for FUNC in (:cp, :cv, :u, :h, :s0)
    @eval begin
        $FUNC(G::IdealGas, T::Real, B::Symbol) = $FUNC(G.cpMod, T, B)
    end
end

for FUNC in (:gamma, )
    @eval begin
        $FUNC(G::IdealGas, T::Real) = $FUNC(G.cpMod, T)
    end
end

for FUNC in (:R, )
    @eval begin
        $FUNC(G::IdealGas, B::Symbol) = $FUNC(G.cpMod, B)
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
