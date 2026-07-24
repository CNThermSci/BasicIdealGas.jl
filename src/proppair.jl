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

# External constructors
# ---------------------

# Set type conversion / 1 indirection
PropPair{ℙ}(
    P::Real,
    T::Real,
) where {ℙ} = PropPair(ℙ.((P, T))...)

# Ideal gas model type conversion / 2 indirections
function PropPair(
        P::Real,
        T::Real,
    )
    ℙ = promote_type(typeof.((P, T))...)
    ℙ = ℙ <: FLOAT ? ℙ : Float64
    return PropPair{ℙ}(P, T)
end

# Set type with unit conversion and stripping / 2 indirections
function PropPair{ℙ}(
        P::PRES,
        T::TEMP,
    ) where {ℙ <: FLOAT}
    return PropPair{ℙ}(kSI.((P, T))...)
end

# Heat model type with unit conversion and stripping / 3 indirections
function PropPair(
        P::PRES{𝔸},
        T::TEMP{𝔹},
    ) where {𝔸 <: Real, 𝔹 <: Real}
    ℙ = promote_type(𝔸, 𝔹)
    ℙ = ℙ <: FLOAT ? ℙ : Float64
    return PropPair{ℙ}(P, T)
end

# Conversions
# -----------

import Base: convert

convert(::Type{PropPair{ℙ}}, ξ::PropPair{ℙ}) where {ℙ <: FLOAT} = ξ

function convert(::Type{PropPair{ℙ}}, ξ::PropPair{ℚ}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return PropPair{ℙ}(ξ.𝑃, ξ.𝑇)
end

import Base: Float16, Float32, Float64

Float16(ξ::PropPair) = convert(PropPair{Float16}, ξ)
Float32(ξ::PropPair) = convert(PropPair{Float32}, ξ)
Float64(ξ::PropPair) = convert(PropPair{Float64}, ξ)

# Promotions
# ----------

import Base: promote_rule

function promote_rule(::Type{PropPair{ℙ}}, ::Type{PropPair{ℚ}}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return PropPair{promote_type(ℙ, ℚ)}
end

# Functor
# -------

function (ξ::PropPair{ℙ})(
        ;
        P::Union{Missing, Real, PRES} = missing,
        T::Union{Missing, Real, TEMP} = missing,
    ) where {ℙ}
    return if count(x -> !isa(x, Missing), (P, T)) == 0
        # pars variant
        pairs(
            (
                P = ξ.P,
                T = ξ.T,
            )
        )
    else
        # copy-edit variant
        PropPair{ℙ}(
            P isa Missing ? ξ.𝑃 : P,
            T isa Missing ? ξ.𝑇 : T,
        )
    end
end

# Export
# ------

export PropPair

# Base.getproperty
# ----------------

import Base: getproperty, propertynames

function Base.getproperty(ξ::PropPair, sy::Symbol)
    # Raw fields
    if sy in fieldnames(PropPair)
        return getfield(ξ, sy)
    end
    # User-facing state function accessors (with units)
    if sy == :P
        return getfield(ξ, :𝑃) * u"kPa"
    elseif sy == :T
        return getfield(ξ, :𝑇) * u"K"
    end
end

Base.propertynames(ξ::PropPair) = (
    fieldnames(PropPair)...,
    :P, :T,
)

# User-facing functions
# ---------------------

function Base.show(io::IO, ::MIME"text/plain", st::PropPair{ℙ}) where {ℙ <: FLOAT}
    return print(
        io,
        "@$(pDeco(ℙ))($(@sprintf("%.*g", 5, st.𝑃)) kPa, $(@sprintf("%.*g", 5, st.𝑇)) K)"
    )
end
