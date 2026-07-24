# proppair.jl - Property Pair

# Structure (type) definition
# ---------------------------

struct PropPair{ℙ <: FLOAT}
    𝑃::ℙ                # kPa
    𝑇::ℙ                # K
    # Internal, validating constructors
    function PropPair(
            P::ℙ,
            T::ℙ,
        ) where {ℙ <: FLOAT}
        @assert(P > zero(ℙ), "P = $(P) ⩽ 0")
        @assert(T > zero(ℙ), "T = $(T) ⩽ 0")
        return new{ℙ}(P, T)
    end
end
