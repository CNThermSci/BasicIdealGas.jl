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
⊚(p::Type{ℙ}, c::ComposedFunction{Type{ℚ}}) where {ℙ <: FLOAT, ℚ <: FLOAT} = ⊚(p, c.inner)

# Auxiliary methods
function subscript(x::Int)
    asSub(c::Char) = Char(Int(c) - Int('0') + Int('₀'))
    return map(asSub, "$(x)")
end

pDeco(::Type{Float16}) = subscript(16)
pDeco(::Type{Float32}) = subscript(32)
pDeco(::Type{Float64}) = subscript(64)
