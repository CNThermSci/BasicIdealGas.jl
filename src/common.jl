# Type aliasing
# -------------

# IEEE-754 normalized floating point types of half, single, and double precision
FLOAT = Base.IEEEFloat

# Thermodynamic state function Quantity type alias
PRES = Quantity{ℙ, dimension(u"kPa")} where {ℙ <: Real}
TEMP = Quantity{ℙ, dimension(u"K")} where {ℙ <: Real}
MOLW = Quantity{ℙ, dimension(u"kg/kmol")} where {ℙ <: Real}
VOLU = Union{
    Quantity{ℙ, dimension(u"m^3/kg")},
    Quantity{ℙ, dimension(u"m^3/kmol")},
} where {ℙ <: Real}
ENER = Union{
    Quantity{ℙ, dimension(u"kJ/kg")},
    Quantity{ℙ, dimension(u"kJ/kmol")},
} where {ℙ <: Real}
ENTR = Union{
    Quantity{ℙ, dimension(u"kJ/kg/K")},
    Quantity{ℙ, dimension(u"kJ/kmol/K")},
} where {ℙ <: Real}
DENS = Union{
    Quantity{ℙ, dimension(u"kg/m^3")},
    Quantity{ℙ, dimension(u"kmol/m^3")},
} where {ℙ <: Real}

# Termodynamic base Unions
MASS = Union{
    Quantity{ℙ, dimension(u"m^3/kg")},
    Quantity{ℙ, dimension(u"kJ/kg")},
    Quantity{ℙ, dimension(u"kJ/kg/K")},
    Quantity{ℙ, dimension(u"kg/m^3")},
}
MOLR = Union{
    Quantity{ℙ, dimension(u"m^3/kmol")},
    Quantity{ℙ, dimension(u"kJ/kmol")},
    Quantity{ℙ, dimension(u"kJ/kmol/K")},
    Quantity{ℙ, dimension(u"kmol/m^3")},
}

# Constants
# ---------

universal_R = 8.31447

# Utilities
# ---------

# Precision Composition Simplification
⊚(p::Type{ℙ}, f::Function) where {ℙ <: FLOAT} = f(1) isa ℙ ? f : p ∘ f

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

function ∫(
        𝑔::Function,
        a::Union{Float32, Float64, Integer, Rational},
        b::Union{Float32, Float64, Integer, Rational},
    )
    ℙ = typeof(promote(a, b)[1])
    return quadgk(𝑔, a, b, rtol = eps(ℙ) * 2 << 6)[1]
end

function ∫(
        𝑔::Function,
        a::Union{Float16, Integer, Rational},
        b::Union{Float16, Integer, Rational},
    )
    a32, b32 = Float32.((a, b))
    n = max(Int(ceil((b32 - a32) / 0.25f0)), 32)
    x32 = range(a32, step = (b32 - a32) / n, length = n + 1) |> collect
    y32 = map(𝑔, x32)
    return Float16(integrate(x32, y32, Trapezoidal()))
end
