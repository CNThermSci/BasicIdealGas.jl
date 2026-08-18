# interact.jl - Heat-work interactions

# Structure (type) definition
# ---------------------------

struct Interact{ℙ <: FLOAT}
    𝑞::ℙ            # specific heat interaction, kJ/kg
    𝑤::ℙ            # specific work interaction, kJ/kg
    # Internal constructors (no validation needed beyond normal julia)
    function Interact(
            q::ℙ,
            w::ℙ,
        ) where {ℙ <: FLOAT}
        return new{ℙ}(q, w)
    end
end

# External constructors
# ---------------------

# Set type conversion / 1 indirection
Interact{ℙ}(
    q::Real,
    w::Real,
) where {ℙ} = Interact(ℙ.((q, w))...)

# Conversions
# -----------

import Base: convert, promote_rule
import Base: Float16, Float32, Float64
import Base: getproperty, propertynames
export Interact

convert(::Type{Interact{ℙ}}, ξ::Interact{ℙ}) where {ℙ <: FLOAT} = ξ

function convert(::Type{Interact{ℙ}}, ξ::Interact{ℚ}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return Interact{ℙ}(ξ.𝑞, ξ.𝑤)
end

Float16(ξ::Interact) = convert(Interact{Float16}, ξ)
Float32(ξ::Interact) = convert(Interact{Float32}, ξ)
Float64(ξ::Interact) = convert(Interact{Float64}, ξ)

# Promotions
# ----------

function promote_rule(::Type{Interact{ℙ}}, ::Type{Interact{ℚ}}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return Interact{promote_type(ℙ, ℚ)}
end
function (ξ::Interact{ℙ})(
        ;
        q::Union{Missing, Real, ENER} = missing,
        w::Union{Missing, Real, ENER} = missing,
    ) where {ℙ}
    return if count(x -> !isa(x, Missing), (q, w)) == 0
        pairs((q = ξ.𝑞, w = ξ.𝑤))
    else
        Interact{ℙ}(
            q isa Missing ? ξ.𝑞 : q,
            w isa Missing ? ξ.𝑤 : w,
        )
    end
end
function Base.getproperty(ξ::Interact, sy::Symbol)
    # Raw fields
    if sy in fieldnames(Interact)
        return getfield(ξ, sy)
    end
    # User-facing state function accessors (with units)
    if sy == :q
        return getfield(ξ, :𝑞) * u"kJ/kg"
    elseif sy == :w
        return getfield(ξ, :𝑤) * u"kJ/kg"
    end
end
Base.propertynames(ξ::Interact) = (fieldnames(Interact)..., :q, :w)

