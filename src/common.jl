# Type aliasing
# -------------

# IEEE-754 normalized floating point types of half, single, and double precision
FLOAT = Base.IEEEFloat

# Utilities
# ---------

# Precision Composition Simplification
⊚(p::Type{ℙ}, f::Function) where {ℙ <: FLOAT} = f(1) isa ℙ ? f : p ∘ f

# Chained Precision Composition Simplification
⊚(p::Type{ℙ}, c::ComposedFunction{Type{ℚ}}) where {ℙ <: FLOAT, ℚ <: FLOAT} = ⊚(p, c.inner)
