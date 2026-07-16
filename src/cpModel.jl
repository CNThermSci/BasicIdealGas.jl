# cpModel.jl - General Specific Heat Model

# Structure (type) definition
# ---------------------------

struct SpecificHeat{ℙ <: FLOAT}
    ID::Symbol      # Model ID, as in :cubic, etc...
    𝑓::Function     # K -> kJ/kmol/K
    𝑀::ℙ            # kg/kmol
    Tmin::ℙ         # K
    Tref::ℙ         # K
    Tmax::ℙ         # K
    uref::ℙ         # kJ/kmol
    sref::ℙ         # kJ/kmol⋅K
    𝑅::ℙ            # kJ/kmol⋅K
    # Internal constructors
    # Validating
    SpecificHeat(
        ID::Symbol,
        𝑓::Function,
        𝑀::ℙ,
        Tmin::ℙ,
        Tref::ℙ,
        Tmax::ℙ,
        uref::ℙ,
        sref::ℙ,
        𝑅::ℙ = ℙ(universal_R),
        B::Symbol = :MO
    ) where {ℙ <: FLOAT} = begin
        @assert(ID != Symbol(""), "Error: Empty model ID")
        @assert(𝑀 > zero(ℙ), "Error: M <= 0")
        @assert(zero(ℙ) <= Tmin <= Tref < Tmax, "Error: Temperature values")
        @assert(𝑅 > zero(ℙ), "Error: 𝑅 <= 0")
        @assert(B in (:MA, :MO), "Error: B should be either :MA or :MO")
        return B == :MA ? (
                new{ℙ}(ID, ℙ ⊚ T -> 𝑓(T) * 𝑀, 𝑀, Tmin, Tref, Tmax, uref * 𝑀, sref * 𝑀, 𝑅)
            ) : (
                new{ℙ}(ID, ℙ ⊚ 𝑓, 𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅)
            )
    end
end

# External constructors
# ---------------------

# Set type conversion / 1 indirection
function SpecificHeat{ℙ}(
        ID::Symbol,
        𝑓::Function,
        𝑀::Real,
        Tmin::Real,
        Tref::Real,
        Tmax::Real,
        uref::Real,
        sref::Real,
        𝑅::Real = ℙ(universal_R),
        B::Symbol = :MO,
    ) where {ℙ <: FLOAT}
    return SpecificHeat(ID, ℙ ⊚ 𝑓, ℙ.((𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅))..., B)
end

# Promotion type conversion / 2 indirections
function SpecificHeat(
        ID::Symbol,
        𝑓::Function,
        𝑀::Real,
        Tmin::Real,
        Tref::Real,
        Tmax::Real,
        uref::Real,
        sref::Real,
        𝑅::Real = universal_R,
        B::Symbol = :MO,
    )
    ℙ = promote_type(typeof.((𝑀, Tmin, Tref, Tmax, uref, sref))...) # Default 𝑅 left out
    ℙ = ℙ <: FLOAT ? ℙ : Float64
    return SpecificHeat{ℙ}(ID, 𝑓, 𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅, B)
end

# Set type with unit conversion and stripping / 2 indirections
function SpecificHeat{ℙ}(
        ID::Symbol,
        𝑓::Function,
        𝑀::Union{Real, MWGT},
        Tmin::Union{Real, TEMP},
        Tref::Union{Real, TEMP},
        Tmax::Union{Real, TEMP},
        uref::ENER,
        sref::ENTR
        𝑅::Union{Real, ENTR} = universal_R,
    ) where {ℙ <: FLOAT}
    _𝑀 = 𝑀 isa Quantity ? uconvert(u"kg/kmol", 𝑀).val : 𝑀
    _uMO = uref isa Quantity{<:Real, dimension(u"kJ/kmol")} ? (
            uconvert(u"kJ/kmol", uref).val
        ) : (
            uconvert(u"kJ/kg", uref).val * _M
        )
    _sMO = sref isa Quantity{<:Real, dimension(u"kJ/kmol/K")} ? (
            uconvert(u"kJ/kmol/K", sref).val
        ) : (
            uconvert(u"kJ/kg/K", sref).val * _M
        )
    return SpecificHeat{ℙ}(
        ID,
        𝑓,
        _𝑀,
        Tmin isa Quantity ? uconvert(u"K", Tmin).val : Tmin,
        Tref isa Quantity ? uconvert(u"K", Tref).val : Tref,
        Tmax isa Quantity ? uconvert(u"K", Tmax).val : Tmax,
        _uMO,
        _sMO,
        𝑅 isa Quantity ? uconvert(u"kJ/kmol/K", 𝑅).val : 𝑅,
        :MO,
    )
end

# Promotion type with unit conversion and stripping / 3 indirections
function SpecificHeat(
        ID::Symbol,
        𝑓::Function,
        𝑀::Union{𝕄, Quantity{<:𝕄, dimension(u"kg/kmol")}},
        Tmin::Union{𝕀, Quantity{<:𝕀, dimension(u"K")}},
        Tref::Union{𝔸, Quantity{<:𝔸, dimension(u"K")}},
        Tmax::Union{𝔼, Quantity{<:𝔼, dimension(u"K")}},
        uref::Union{
            Quantity{<:𝕌, dimension(u"kJ/kmol")},
            Quantity{<:𝕌, dimension(u"kJ/kg")},
        },
        sref::Union{
            Quantity{<:𝕊, dimension(u"kJ/kmol/K")},
            Quantity{<:𝕊, dimension(u"kJ/kg/K")},
        },
        𝑅::Union{Real, Quantity{<:Real, dimension(u"kJ/kmol/K")}} = universal_R,
    ) where {𝕄 <: Real, 𝕀 <: Real, 𝔸 <: Real, 𝔼 <: Real, 𝕌 <: Real, 𝕊 <: Real}
    ℙ = promote_type(𝕄, 𝕀, 𝔸, 𝔼, 𝕌, 𝕊) # Default R left out
    ℙ = ℙ <: FLOAT ? ℙ : Float64
    return SpecificHeat{ℙ}(ID, 𝑓, 𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅)
end

# Conversions
# -----------

import Base: convert

convert(::Type{SpecificHeat{ℙ}}, ξ::SpecificHeat{ℙ}) where {ℙ <: FLOAT} = ξ

function convert(::Type{SpecificHeat{ℙ}}, ξ::SpecificHeat{ℚ}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return SpecificHeat{ℙ}(
        ξ.ID, ξ.𝑓, ξ.𝑀, ξ.Tmin, ξ.Tref, ξ.Tmax, ξ.uref, ξ.sref, ξ.𝑅
    )
end

import Base: Float16, Float32, Float64

Float16(ξ::SpecificHeat) = convert(SpecificHeat{Float16}, ξ)
Float32(ξ::SpecificHeat) = convert(SpecificHeat{Float32}, ξ)
Float64(ξ::SpecificHeat) = convert(SpecificHeat{Float64}, ξ)

# Promotions
# ----------

import Base: promote_rule

function promote_rule(
        ::Type{SpecificHeat{ℙ}},
        ::Type{SpecificHeat{ℚ}}
    ) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return SpecificHeat{promote_type(ℙ, ℚ)}
end

# Export
# ------

export SpecificHeat

# Show
# ----

function Base.show(io::IO, S::SpecificHeat{ℙ}) where {ℙ <: FLOAT}
    rng = "[$(@sprintf("%.*g K", 5, S.Tmin)) $(@sprintf("%.*g K", 5, S.Tmax))]"
    return print(io, "$(S.ID) cp$(pDeco(ℙ))(T) $(rng)")
end

# SpecificHeat Helper functions
# -----------------------------

∫┆T(C::SpecificHeat, T::Real) = ∫(T -> C.𝑓(T) / T, C.Tref, T)

# User-facing functions
# ---------------------

𝗯(C::SpecificHeat, T::Real) = @assert(C.Tmin <= T <= C.Tmax, "T out of bounds")

import Base: cp

cp┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = (𝗯(C, T); C.𝑓(T) / C.𝑅)
cv┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = cp┆R(C, T) - one(ℙ)
ga(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = (𝗯(C, T); x = C.𝑓(T); x / (x - C.𝑅))

export cp┆R, cv┆R, ga

function R(C::SpecificHeat, B::Symbol = :MA)
    @assert B in (:MA, :MO)
    return B == :MO ? C.𝑅 : C.𝑅 / C.𝑀
end
cp(C::SpecificHeat{ℙ}, T::Real, B::Symbol = :MA) where {ℙ <: FLOAT} = cp┆R(C, T) * R(C, B)
cv(C::SpecificHeat{ℙ}, T::Real, B::Symbol = :MA) where {ℙ <: FLOAT} = cv┆R(C, T) * R(C, B)

export R, cp, cv

∫cp┆R(C::SpecificHeat{ℙ}, T::ℙ) where {ℙ <: FLOAT} = (𝗯(C, T); ∫(C.𝑓, C.Tref, T) / C.𝑅)
∫cp┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = ∫cp┆R(C, ℙ(T))
∫cv┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = ∫cp┆R(C, T) - ℙ(T) + C.Tref
u┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = ∫cv┆R(C, T) + C.uref / C.𝑅
h┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = u┆R(C, T) + ℙ(T)

export ∫cp┆R, ∫cv┆R, u┆R, h┆R

u(C::SpecificHeat{ℙ}, T::Real, B::Symbol = :MA) where {ℙ <: FLOAT} = u┆R(C, T) * R(C, B)
h(C::SpecificHeat{ℙ}, T::Real, B::Symbol = :MA) where {ℙ <: FLOAT} = h┆R(C, T) * R(C, B)

export u, h

∫cp┆RT(C::SpecificHeat, T::Real) = (𝗯(C, T); ∫┆T(C, T) / C.𝑅)
s0┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = ∫cp┆RT(C, T) + C.sref / C.𝑅

export ∫cp┆RT, s0┆R

s0(C::SpecificHeat{ℙ}, T::Real, B::Symbol = :MA) where {ℙ <: FLOAT} = s0┆R(C, T) * R(C, B)
Pr(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = exp(∫cp┆RT(C, T))
vr(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = ℙ(T) / Pr(C, T)

export s0, Pr, vr

# Base.getproperty
# ----------------

import Base: getproperty, propertynames

function Base.getproperty(ξ::SpecificHeat, sy::Symbol)
    # Raw fields
    if sy in fieldnames(SpecificHeat) return getfield(ξ, sy) end
    # Convenience accessors/transformers
    if sy in (:f, :mod, :modMO)
        return getfield(ξ, :𝑓)
    elseif sy in (:fMA, :modMA)
        return T -> getfield(ξ, :𝑓)(T) / getfield(ξ, :𝑀)
    end
    # Porcelain accessors (with units)
    if sy == :M
        return getfield(ξ, :𝑀) * u"kg/kmol"
    elseif sy in (:R, :RMA)
        return R(ξ, :MA) * u"kJ/kg/K"
    elseif sy in (:RU, :RMO)
        return getfield(ξ, :𝑅) * u"kJ/kmol/K"
    end
    # Pretty print
    if sy == :view
        xmin, xmax = getfield(ξ, :Tmin), getfield(ξ, :Tmax)
        x = range(xmin, stop = xmax, length = 33)
        y = map(T -> cp(ξ, T, :MA), x)
        plt = lineplot(
            x, y, xlabel = "T [K]", ylabel = "cp (T)", name = "⠤⠤⠤⠤ [kJ/kg·K]",
            xlim = (xmin, xmax), width = 32, height = 6,
            border = :ascii, color = :white, compact_labels = true,
        )
        print(join([repr(ξ), string(plt)], "\n"))
    end
end

Base.propertynames(::SpecificHeat) = (
    :ID, :𝑓, :𝑀, :Tmin, :Tmax, :Tref, :uref, :sref, :𝑅,
    :f, :mod, :modMO, :fMA, :modMA, :M, :R, :RMA, :RU, :RMO,
    :view,
)
