# Type aliasing
# -------------

# IEEE-754 normalized floating point types of half, single, and double precision
FLOAT = Base.IEEEFloat

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

# Thermodynamic unit conversion/stripping
kSI(x::Real) = x
kSI(x::PRES) = uconvert(u"kPa", x).val
kSI(x::TEMP) = uconvert(u"K", x).val
kSI(x::MOLW) = uconvert(u"kg/kmol", x).val

function kSI(x::MASS)
    return if x isa VOLU
        uconvert(u"m^3/kg", x).val
    elseif x isa ENER
        uconvert(u"kJ/kg", x).val
    elseif x isa ENTR
        uconvert(u"kJ/kg/K", x).val
    elseif x isa DENS
        uconvert(u"kg/m^3", x).val
    end
end

function kSI(x::MOLR)
    return if x isa VOLU
        uconvert(u"m^3/kmol", x).val
    elseif x isa ENER
        uconvert(u"kJ/kmol", x).val
    elseif x isa ENTR
        uconvert(u"kJ/kmol/K", x).val
    elseif x isa DENS
        uconvert(u"kmol/m^3", x).val
    end
end

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
    a ≈ b && return zero(Float16)
    sa, sb, ss = b > a ? Float16.((a, b, 1)) : Float16.((b, a, -1))
    n = min(Int(trunc((sb - sa) / eps(sb))), 256)
    x = range(sa, step = (sb - sa) / n, length = n + 1) |> collect
    y = map(𝑔, x)
    return integrate(x, y, Trapezoidal()) * ss
end

## function ∫(𝑔::𝔽, a::Quantity{𝔸}, b::Quantity{𝔹}) where {𝔽 <: Function, 𝔸 <: LOTYP, 𝔹 <: LOTYP}
##     return ∫(𝑔.f┆R, a.val, b.val) * unit(a) * unit(𝑔(a))
## end
