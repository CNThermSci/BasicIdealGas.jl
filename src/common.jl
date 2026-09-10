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
⊚(p::Type{ℙ}, f::Function) where {ℙ <: FLOAT} = f(ℙ(300u"K")) isa ℙ ? f : p ∘ f

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

HILIM = Union{
    Quantity{Float32}, Quantity{Float64}, Quantity{Integer}, Quantity{Rational},
    Float32, Float64, Integer, Rational,
}

LOLIM = Union{
    Quantity{Float16}, Quantity{Integer}, Quantity{Rational},
    Float16, Integer, Rational,
}

function ∫(𝑔::Function, a::HILIM, b::HILIM)
    ℙ = typeof(promote(a, b, one(Float32))[1])
    return quadgk(𝑔, a, b, rtol = eps(ℙ) * 2 << 6)[1]
end

function ∫(𝑔::Function, a::LOLIM, b::LOLIM)
    a32, b32 = Float32.((a, b))
    n = max(Int(ceil((b32 - a32) / 0.25f0)), 32)
    x32 = range(a32, step = (b32 - a32) / n, length = n + 1) |> collect
    y32 = map(𝑔, x32)
    return Float16(integrate(x32, y32, Trapezoidal()))
end
