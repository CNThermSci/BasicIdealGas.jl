# idealState.jl - Ideal Gas States

# Structure (type) definition
# ---------------------------

struct IdealState{ℙ <: FLOAT}
    𝐺::IdealGas{ℙ}
    𝑃::ℙ
    𝑇::ℙ
    # Internal, validating constructors
    function IdealState(
            G::IdealGas{ℙ},
            P::ℙ,
            T::ℙ,
        ) where {ℙ <: FLOAT}
        new{ℙ}(G, P, T)
    end
end

# Base.getproperty
# ----------------

import Base: getproperty, propertynames

function Base.getproperty(obj::IdealState, s::Symbol)
    if s in (:𝐺, :𝑃, :𝑇)
        return getfield(obj, s)
    elseif s == :gas
        return getfield(obj, :𝐺)
    elseif s == :P
        return getfield(obj, :𝑃) * u"kPa"
    elseif s == :T
        return getfield(obj, :𝑇) * u"K"
    end
end

Base.propertynames(::IdealState) = (:gas, :P, :T)

# Export
# ------

export IdealState

# "ﬆ" is U+FB06
