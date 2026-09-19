# idealEoS.jl - Ideal Gas Equation of State

# Structure (type) definition
# ---------------------------

struct IdealGas{ℙ <: FLOAT}
    form::String            # formula
    name::String            # name
    hmod::SpecificHeat{ℙ}   # heat model
    𝑃ref::Quantity{ℙ, dimension(u"kPa"), typeof(u"kPa")}
    function IdealGas(
            FORM::AbstractString,
            NAME::AbstractString,
            HMOD::SpecificHeat{ℙ},
            PREF::Quantity{ℙ, dimension(u"kPa"), typeof(u"kPa")} = one(ℙ) * u"kPa"
        ) where {ℙ <: FLOAT}
        @assert(length(FORM) > 0, "Error: Empty formula")
        @assert(length(NAME) > 0, "Error: Empty name")
        @assert(PREF > 0u"kPa", "Error: Pref <= 0 kPa")
        return new{ℙ}(String(FORM), String(NAME), HMOD, PREF)
    end
end

# External constructors
# ---------------------

# Set precision conversion / 1 indirection
function IdealGas{ℙ}(
        FORM::AbstractString,
        NAME::AbstractString,
        HMOD::SpecificHeat,
        PREF::Union{Real, PRES} = one(ℙ) * u"kPa",
    ) where {ℙ}
    Pref = PREF isa PRES ? uconvert(u"kPa", PREF) : PREF * u"kPa"
    return IdealGas(FORM, NAME, ℙ.((HMOD, Pref))...)
end

# Heat model type conversion / 2 indirections
function IdealGas(
        FORM::AbstractString,
        NAME::AbstractString,
        HMOD::SpecificHeat{ℙ},
        PREF::Union{Real, PRES} = one(ℙ) * u"kPa",
    ) where {ℙ}
    return IdealGas{ℙ}(FORM, NAME, HMOD, PREF)
end

# Conversions
# -----------

import Base: convert

convert(::Type{IdealGas{ℙ}}, ξ::IdealGas{ℙ}) where {ℙ <: FLOAT} = ξ

function convert(::Type{IdealGas{ℙ}}, ξ::IdealGas{ℚ}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return IdealGas{ℙ}(ξ.form, ξ.name, ξ.hmod, ξ.𝑃ref)
end

import Base: Float16, Float32, Float64

Float16(ξ::IdealGas) = convert(IdealGas{Float16}, ξ)
Float32(ξ::IdealGas) = convert(IdealGas{Float32}, ξ)
Float64(ξ::IdealGas) = convert(IdealGas{Float64}, ξ)

# Promotions
# ----------

import Base: promote_rule

function promote_rule(::Type{IdealGas{ℙ}}, ::Type{IdealGas{ℚ}}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return IdealGas{promote_type(ℙ, ℚ)}
end

# Export
# ------

export IdealGas

# User-facing functions
# ---------------------

function Base.show(io::IO, ::MIME"text/plain", ξ::IdealGas{ℙ}) where {ℙ <: FLOAT}
    return print(
        io,
        "$(ξ.form) gas, ",
        repr(MIME"text/plain"(), ξ.hmod),
    )
end

for FUNC in (:R,)
    @eval begin
        $FUNC(ξ::IdealGas, B::Symbol = :MA) = $FUNC(ξ.hmod, B)
    end
end

R(ξ::IdealGas, 𝑣::VOLU) = 𝑣 isa MASS ? R(ξ, :MA) : R(ξ, :MO)

for FUNC in (:cp┆R, :cv┆R, :ga, :∫cp┆R, :∫cv┆R, :u┆R, :h┆R, :∫cp┆RT, :s0┆R, :Pr, :vr)
    @eval begin
        $FUNC(ξ::IdealGas, 𝑇::TEMP) = $FUNC(ξ.hmod, 𝑇)
        $FUNC(ξ::IdealGas, 𝑇::Real) = $FUNC(ξ.hmod, 𝑇)
    end
end

for FUNC in (:cp, :cv, :u, :h, :s0)
    @eval begin
        $FUNC(ξ::IdealGas, 𝑇::TEMP, B::Symbol = :MA) = $FUNC(ξ.hmod, 𝑇, B)
        $FUNC(ξ::IdealGas, 𝑇::Real, B::Symbol = :MA) = $FUNC(ξ.hmod, 𝑇, B)
    end
end

# Internal Positional P, T, V, ρ, s functions
# -------------------------------------------

# Pressure
function _P(ξ::IdealGas{ℙ}, 𝑇::TEMP, 𝑣::VOLU) where {ℙ}
    return uconvert(u"kPa", R(ξ, 𝑣) * ℙ(𝑇 / 𝑣))
end
_P(ξ::IdealGas{ℙ}, 𝑃::PRES) where {ℙ} = uconvert(u"kPa", ℙ(𝑃))
_P(ξ::IdealGas{ℙ}, P::Real) where {ℙ} = ℙ(P) * u"kPa"

# Temperature
function _T(ξ::IdealGas{ℙ}, 𝑃::PRES, 𝑣::VOLU) where {ℙ}
    return uconvert(u"K", ℙ(𝑃 * 𝑣) / R(ξ, 𝑣))
end
_T(ξ::IdealGas{ℙ}, 𝑇::TEMP) where {ℙ} = uconvert(u"K", ℙ(T))
_T(ξ::IdealGas{ℙ}, T::Real) where {ℙ} = ℙ(T) * u"K"

# Specific volume
function _v(ξ::IdealGas{ℙ}, 𝑃::PRES, 𝑇::TEMP, B::Symbol = :MA) where {ℙ}
    UNIT = B == :MA ? u"m^3/kg" : u"m^3/kmol"
    return uconvert(UNIT, R(ξ, B) * ℙ(𝑇 / 𝑃))
end
function _v(ξ::IdealGas{ℙ}, 𝑣::VOLU) where {ℙ}
    UNIT = 𝑣 isa MASS ? u"m^3/kg" : u"m^3/kmol"
    return uconvert(UNIT, ℙ(𝑣))
end
function _v(ξ::IdealGas{ℙ}, 𝑣::VOLU, B::Symbol) where {ℙ}
    iUNIT = 𝑣 isa MASS ? u"m^3/kg" : u"m^3/kmol"
    oUNIT = B == :MA ? u"m^3/kg" : u"m^3/kmol"
    return if iUNIT == oUNIT
        # No base change / same dims / units can still differ
        uconvert(oUNIT, ℙ(𝑣))
    else
        # Base change to B
        uconvert(oUNIT, B == :MA ? 𝑣 / ξ.𝑀 : 𝑣 * ξ.𝑀)
    end
end
function _v(ξ::IdealGas{ℙ}, v::Real, B::Symbol = :MA) where {ℙ}
    return ℙ(v) * (B == :MA ? u"m^3/kg" : u"m^3/kmol")
end

# Density
_ρ(ξ::IdealGas, 𝑃::PRES, 𝑇::TEMP, B::Symbol = :MA) = inv(_v(ξ, 𝑃, 𝑇, B))

# Specific entropy
function _s(ξ::IdealGas{ℙ}, 𝑃::PRES, 𝑇::TEMP, B::Symbol = :MA) where {ℙ}
    return s0(ξ, 𝑇, B) - R(ξ, B) * log(ℙ(𝑃) / ξ.𝑃ref)
end

# PTV helper
# ----------

PTv = (
    ξ::IdealGas;
    P::Union{Real, PRES, Missing} = missing,
    T::Union{Real, TEMP, Missing} = missing,
    v::Union{Real, VOLU, Missing} = missing,
    B::Symbol = :MA,
) -> begin
    miss = [ i[1] for i in [(:P, P), (:T, T), (:v, v)] if ismissing(i[2]) ]
    length(miss) <= 1 ||
        throw(ArgumentError(@sprintf("Underspecified state: missing (%s)", join(miss, ", "))))
    return if ismissing(P)
        𝑇 = _T(ξ, T)
        𝑣 = _v(ξ, v, B)
        _P(ξ, 𝑇, 𝑣), 𝑇, 𝑣
    elseif ismissing(T)
        𝑃 = _P(ξ, P)
        𝑣 = _v(ξ, v, B)
        𝑃, _T(ξ, 𝑃, 𝑣), 𝑣
    else
        𝑃 = _P(ξ, P)
        𝑇 = _T(ξ, T)
        𝑃, 𝑇, _v(ξ, 𝑃, 𝑇, B)
    end
end

kws = (
    ξ::IdealGas;
    P::Union{Real, PRES, Missing} = missing,
    T::Union{Real, TEMP, Missing} = missing,
    v::Union{Real, VOLU, Missing} = missing,
    B::Symbol = :MA,
) -> begin
    𝑃, 𝑇, 𝑣 = PTv(ξ, P = P, T = T, v = v, B = B)
    _s(ξ, 𝑃, 𝑇, B)
end

kwa = (
    ξ::IdealGas;
    P::Union{Real, PRES, Missing} = missing,
    T::Union{Real, TEMP, Missing} = missing,
    v::Union{Real, VOLU, Missing} = missing,
    B::Symbol = :MA,
) -> begin
    𝑃, 𝑇, 𝑣 = PTv(ξ, P = P, T = T, v = v, B = B)
    𝑠 = _s(ξ, 𝑃, 𝑇, B)
    𝑢 = u(ξ, 𝑇, B)
    𝑢 - 𝑇 * 𝑠
end

kwg = (
    ξ::IdealGas;
    P::Union{Real, PRES, Missing} = missing,
    T::Union{Real, TEMP, Missing} = missing,
    v::Union{Real, VOLU, Missing} = missing,
    B::Symbol = :MA,
) -> begin
    𝑃, 𝑇, 𝑣 = PTv(ξ, P = P, T = T, v = v, B = B)
    𝑠 = _s(ξ, 𝑃, 𝑇, B)
    ℎ = h(ξ, 𝑇, B)
    ℎ - 𝑇 * 𝑠
end

kwβ = (
    ξ::IdealGas;
    P::Union{Real, PRES, Missing} = missing,
    T::Union{Real, TEMP, Missing} = missing,
    v::Union{Real, VOLU, Missing} = missing,
    B::Symbol = :MA,
) -> begin
    if !ismissing(T)
        return inv(_T(ξ, T))
    end
    𝑃, 𝑇, 𝑣 = PTv(ξ, P = P, T = T, v = v, B = B)
    inv(𝑇)
end

kwκT = (
    ξ::IdealGas;
    P::Union{Real, PRES, Missing} = missing,
    T::Union{Real, TEMP, Missing} = missing,
    v::Union{Real, VOLU, Missing} = missing,
    B::Symbol = :MA,
) -> begin
    if !ismissing(P)
        return inv(_P(ξ, P))
    end
    𝑃, 𝑇, 𝑣 = PTv(ξ, P = P, T = T, v = v, B = B)
    inv(𝑃)
end

kwκs = (
    ξ::IdealGas;
    P::Union{Real, PRES, Missing} = missing,
    T::Union{Real, TEMP, Missing} = missing,
    v::Union{Real, VOLU, Missing} = missing,
    B::Symbol = :MA,
) -> begin
    𝑃, 𝑇, 𝑣 = PTv(ξ, P = P, T = T, v = v, B = B)
    inv(𝑃 * ga(ξ, 𝑇))
end

# Base.getproperty
# ----------------

import Base: getproperty, propertynames

function Base.getproperty(ξ::IdealGas, sy::Symbol)
    # Raw fields
    if sy in fieldnames(IdealGas)
        return getfield(ξ, sy)
    end
    # Convenience raw field aliases
    if sy == :Pref
        return getfield(ξ, :𝑃ref)
    end
    # Short-circuit SpecificHeat model accessors
    𝐶 = getfield(ξ, :hmod)
    if sy in propertynames(𝐶)
        if sy ∉ (props_T(𝐶)..., props_T_B(𝐶)...)
            return getproperty(𝐶, sy)
        elseif sy ∈ props_T(𝐶)
            # This allows an 𝑓(T) be calc'd from 𝑓(T(P, T, v, B))
            # Makes sense only at the IdealGas level (can't fallback to SpecificHeat directly)
            return (; P = missing, T = missing, v = missing, B = :MA) -> begin
                getproperty(𝐶, sy)(PTv(ξ; P = P, T = T, v = v, B = B)[2])
            end
        elseif sy ∈ props_T_B(𝐶)
            # This allows an 𝑓(T, B) be calc'd from 𝑓(T(P, T, v, B), B)
            # Makes sense only at the IdealGas level (can't fallback to SpecificHeat directly)
            return (; P = missing, T = missing, v = missing, B = :MA) -> begin
                getproperty(𝐶, sy)(PTv(ξ; P = P, T = T, v = v, B = B)[2], B)
            end
        end
    end
    # OOP-style covenience functions (formerly exported ones)
    return if sy == :P
        (; P = missing, T = missing, v = missing, B = :MA) -> PTv(ξ; P = P, T = T, v = v, B = B)[1]
    elseif sy == :T
        (; P = missing, T = missing, v = missing, B = :MA) -> PTv(ξ; P = P, T = T, v = v, B = B)[2]
    elseif sy == :v
        (; P = missing, T = missing, v = missing, B = :MA) -> PTv(ξ; P = P, T = T, v = v, B = B)[3]
    elseif sy == :vMA
        (; P = missing, T = missing, v = missing) -> PTv(ξ; P = P, T = T, v = v, B = :MA)[3]
    elseif sy == :vMO
        (; P = missing, T = missing, v = missing) -> PTv(ξ; P = P, T = T, v = v, B = :MO)[3]
    elseif sy == :ρ
        (; P = missing, T = missing, v = missing, B = :MA) -> inv(PTv(ξ; P = P, T = T, v = v, B = B)[3])
    elseif sy == :ρMA
        (; P = missing, T = missing, v = missing) -> inv(PTv(ξ; P = P, T = T, v = v, B = :MA)[3])
    elseif sy == :ρMO
        (; P = missing, T = missing, v = missing) -> inv(PTv(ξ; P = P, T = T, v = v, B = :MO)[3])
    elseif sy == :s
        (; P = missing, T = missing, v = missing, B = :MA) -> kws(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :sMA
        (; P = missing, T = missing, v = missing) -> kws(ξ; P = P, T = T, v = v, B = :MA)
    elseif sy == :sMO
        (; P = missing, T = missing, v = missing) -> kws(ξ; P = P, T = T, v = v, B = :MO)
    elseif sy == :a
        (; P = missing, T = missing, v = missing, B = :MA) -> kwa(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :aMA
        (; P = missing, T = missing, v = missing) -> kwa(ξ; P = P, T = T, v = v, B = :MA)
    elseif sy == :aMO
        (; P = missing, T = missing, v = missing) -> kwa(ξ; P = P, T = T, v = v, B = :MO)
    elseif sy == :g
        (; P = missing, T = missing, v = missing, B = :MA) -> kwg(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :gMA
        (; P = missing, T = missing, v = missing) -> kwg(ξ; P = P, T = T, v = v, B = :MA)
    elseif sy == :gMO
        (; P = missing, T = missing, v = missing) -> kwg(ξ; P = P, T = T, v = v, B = :MO)
    elseif sy in (:β, :beta)
        (; P = missing, T = missing, v = missing, B = :MA) -> kwβ(ξ; P = P, T = T, v = v, B = B)
    elseif sy in (:κT, :kappaT)
        (; P = missing, T = missing, v = missing, B = :MA) -> kwκT(ξ; P = P, T = T, v = v, B = B)
    elseif sy in (:κs, :kappas)
        (; P = missing, T = missing, v = missing, B = :MA) -> kwκs(ξ; P = P, T = T, v = v, B = B)
    end
end

Base.propertynames(ξ::IdealGas) = (
    :form, :name, :hmod, :𝑃ref,
    :Pref,
    propertynames(getfield(ξ, :hmod))...,
    :P, :T, :v, :ρ, :s,
)
