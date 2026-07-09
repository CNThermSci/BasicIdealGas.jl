# Type aliasing
# -------------

# IEEE-754 normalized floating point types of half, single, and double precision
FLOAT = Base.IEEEFloat

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

function ∫(𝑔::Function, a::ℙ, b::ℙ) where {ℙ <: Union{Float32, Float64}}
    return quadgk(𝑔, a, b, rtol = eps(ℙ) * 2 << 6)[1]
end

function ∫(𝑔::Function, a::Float16, b::Float16)
    n = max(Int(ceil((b - a) / Float16(0.25))), 32)
    x = range(a, step = (b - a) / n, length = n + 1) |> collect
    y = map(𝑔, x)
    return integrate(x, y, Trapezoidal())
end
