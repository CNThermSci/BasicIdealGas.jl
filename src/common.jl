# IEEE-754 normalized floating point types of half, single, and double precision
# ------------------------------------------------------------------------------

FLOAT = Base.IEEEFloat

# Abstract types
# --------------

abstract type ThermModel{ℙ <: FLOAT} end

# Type aliasing
# -------------

# Thermodynamic state function Quantity type alias - dimension set (for arguments)
const PRES = Quantity{ℙ, dimension(u"kPa")} where {ℙ <: Real}
const TEMP = Quantity{ℙ, dimension(u"K")} where {ℙ <: Real}
const MOLW = Quantity{ℙ, dimension(u"kg/kmol")} where {ℙ <: Real}
const VOLU = Union{
    Quantity{ℙ, dimension(u"m^3/kg")},
    Quantity{ℙ, dimension(u"m^3/kmol")},
} where {ℙ <: Real}
const ENER = Union{
    Quantity{ℙ, dimension(u"kJ/kg")},
    Quantity{ℙ, dimension(u"kJ/kmol")},
} where {ℙ <: Real}
const ENTR = Union{
    Quantity{ℙ, dimension(u"kJ/kg/K")},
    Quantity{ℙ, dimension(u"kJ/kmol/K")},
} where {ℙ <: Real}
const DENS = Union{
    Quantity{ℙ, dimension(u"kg/m^3")},
    Quantity{ℙ, dimension(u"kmol/m^3")},
} where {ℙ <: Real}

# Termodynamic base Unions
const MASS = Union{
    Quantity{ℙ, dimension(u"m^3/kg")},
    Quantity{ℙ, dimension(u"kJ/kg")},
    Quantity{ℙ, dimension(u"kJ/kg/K")},
    Quantity{ℙ, dimension(u"kg/m^3")},
} where {ℙ <: Real}
const MOLR = Union{
    Quantity{ℙ, dimension(u"m^3/kmol")},
    Quantity{ℙ, dimension(u"kJ/kmol")},
    Quantity{ℙ, dimension(u"kJ/kmol/K")},
    Quantity{ℙ, dimension(u"kmol/m^3")},
} where {ℙ <: Real}

# Unit-dressing helper functions
𝑃(P::PRES) = uconvert(u"kPa", P)
𝑃(P::Real) = P * u"kPa"
𝑃(ξ::ThermModel{ℙ}, 𝜋::Union{Real, PRES}) where {ℙ} = ℙ(𝑃(𝜋))

𝑇(T::TEMP) = uconvert(u"K", T)
𝑇(T::Real) = T * u"K"
𝑇(ξ::ThermModel{ℙ}, θ::Union{Real, TEMP}) where {ℙ} = ℙ(𝑇(θ))

𝑀(M::MOLW) = uconvert(u"kg/kmol", M)
𝑀(T::Real) = M * u"kg/kmol"
𝑀(ξ::ThermModel{ℙ}, μ::Union{Real, TEMP}) where {ℙ} = ℙ(𝑀(μ))

𝑣(v::VOLU) = v isa MASS ? uconvert(u"m^3/kg", v) : uconvert(u"m^3/kmol", v)
𝑣(v::Tuple{Real, Symbol}) = v[1] * (v[2] == :MA ? u"m^3/kg" : u"m^3/kmol")
𝑣(ξ::ThermModel{ℙ}, υ::Union{Real, VOLU}) where {ℙ} = ℙ(𝑣(υ))

𝑒(e::ENER) = e isa MASS ? uconvert(u"kJ/kg", e) : uconvert(u"kJ/kmol", e)
𝑒(e::Tuple{Real, Symbol}) = e[1] * (e[2] == :MA ? u"kJ/kg" : u"kJ/kmol")
𝑒(ξ::ThermModel{ℙ}, ϵ::Union{Real, ENER}) where {ℙ} = ℙ(𝑒(ϵ))

𝑠(s::ENTR) = s isa MASS ? uconvert(u"kJ/kg/K", s) : uconvert(u"kJ/kmol/K", s)
𝑠(s::Tuple{Real, Symbol}) = s[1] * (s[2] == :MA ? u"kJ/kg/K" : u"kJ/kmol/K")
𝑠(ξ::ThermModel{ℙ}, ς::Union{Real, ENTR}) where {ℙ} = ℙ(𝑠(ς))

# Constants
# ---------

# Exact CODATA2022 value for Ru
const Ru = uconvert(u"kJ/kmol/K", MolarGasConstant)
export Ru

# Legacy CODATA1986 value - baked into NASA-9 coefficients
const RuCODATA1986 = 8.314510u"kJ/kmol/K"
export RuCODATA1986

# Utilities
# ---------

# Precision of
precof(x::Real) = typeof(x)
precof(x::Quantity{ℙ}) where ℙ = ℙ

# Precision Composition Simplification
⊚(p::Type{ℙ}, f::Function) where {ℙ <: FLOAT} = precof(f(300)) == ℙ ? f : p ∘ f

# Chained Precision Composition Simplification
⊚(
    p::Type{ℙ},
    c::ComposedFunction{<:Union{Type{ℚ}, typeof(float)}}
) where {ℙ <: FLOAT, ℚ <: FLOAT} = ⊚(p, c.inner)

# Auxiliary methods
function subscript(x::Int)
    asSub(c::Char) = Char(Int(c) - Int('0') + Int('₀'))
    return map(asSub, "$(x)")
end

pDeco(::Type{Float16}) = subscript(16)
pDeco(::Type{Float32}) = subscript(32)
pDeco(::Type{Float64}) = subscript(64)

# Numerical integrator
# --------------------

HITYP = Union{Integer, Rational, Float32, Float64}
LOTYP = Union{Integer, Rational, Float16}
LOLIM = Union{ℍ, Quantity{ℍ}} where {ℍ <: LOTYP}
HILIM = Union{ℍ, Quantity{ℍ}} where {ℍ <: HITYP}

function ∫(𝑔::𝔽, a::HILIM, b::HILIM) where {𝔽 <: Function}
    ℙ = typeof(promote(a, b, one(Float32))[1])
    return quadgk(𝑔, a, b, rtol = eps(ℙ) * 2 << 6)[1]
end

function ∫(𝑔::𝔽, a::LOLIM, b::LOLIM) where {𝔽 <: Function}
    # isapprox(a, b, atol = max(eps(a), eps(b))) && return zero(Float16) * unit(𝑔(a) * a)
    a != b || return zero(Float16) * unit(𝑔(a) * a)
    sa, sb, ss = b > a ? Float16.((a, b, 1)) : Float16.((b, a, -1))
    n = min(Int(trunc((sb - sa) / eps(sb))), 256)
    x = range(sa, step = (sb - sa) / n, length = n + 1) |> collect
    y = map(𝑔, x)
    return integrate(x, y, Trapezoidal()) * ss
end
