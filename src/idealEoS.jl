# idealEoS.jl - Ideal Gas Equation of State

# Structure (type) definition
# ---------------------------

struct IdealGas{ℙ <: FLOAT} <: ThermModel{ℙ}
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

# Internal, positional, dispatched, no default args P, T, v, ρ functions
# ----------------------------------------------------------------------

P(ξ::IdealGas{ℙ}, 𝑇::TEMP, 𝑣::VOLU) where {ℙ} = PP(R(ξ.hmod, 𝑣) * ℙ(𝑇 / 𝑣))
T(ξ::IdealGas{ℙ}, 𝑃::PRES, 𝑣::VOLU) where {ℙ} = TT(ℙ(𝑃 * 𝑣) / R(ξ.hmod, 𝑣))
v(ξ::IdealGas{ℙ}, 𝑃::PRES, 𝑇::TEMP, B::Symbol) where {ℙ} = vv(R(ξ.hmod, B) * ℙ(𝑇 / 𝑃))
ρ(ξ::IdealGas, 𝑃::PRES, 𝑇::TEMP, B::Symbol) inv(v(ξ, 𝑃, 𝑇, B))
s(ξ::IdealGas{ℙ}, 𝑃::PRES, 𝑇::TEMP, B::Symbol) where {ℙ} = ss(s0(ξ.hmod, 𝑇, B) - R(ξ.hmod, B) * log(ℙ(𝑃) / ξ.𝑃ref))

# Internal, keyworded, default base property functions
# ----------------------------------------------------

# Pressure
function P(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        kw...,
    )
    return ismissing(P) ? BasicIdealGas.P(ξ, TT(ξ, T), vv(ξ, v)) : PP(ξ, P)
end

# Temperature
function T(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        kw...,
    )
    return ismissing(T) ? BasicIdealGas.T(ξ, PP(ξ, P), vv(ξ, v)) : TT(ξ, T)
end

# Specific volume
function v(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        B::Union{Symbol, Missing} = missing,
        kw...,
    )
    return if ismissing(v)
        BasicIdealGas.v(ξ, PP(ξ, P), TT(ξ, T), ismissing(B) ? :MA : B)
    else
        𝑣 = vv(ξ, v)
        if ismissing(B) || (𝑣 isa MASS && B == :MA) || (𝑣 isa MOLR && B == :MO)
            𝑣
        else
            𝑣 isa MASS ? 𝑣 * ξ.M : 𝑣 / ξ.M
        end
    end
end

# Density
function ρ(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        B::Union{Symbol, Missing} = missing,
        kw...,
    )
    return inv(BasicIdealGas.v(ξ; P = P, T = T, v = v, B = B, kw...))
end

# Helper PT function
function PT(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        kw...,
    )
    miss = [ i[1] for i in [(:P, P), (:T, T), (:v, v)] if ismissing(i[2]) ]
    length(miss) <= 1 ||
        throw(ArgumentError(@sprintf("Underspecified state: missing (%s)", join(miss, ", "))))
    return if ismissing(P)
        𝑇 = TT(ξ, T)
        𝑣 = vv(ξ, v)
        BasicIdealGas.P(ξ, 𝑇, 𝑣), 𝑇
    elseif ismissing(T)
        𝑃 = PP(ξ, P)
        𝑣 = vv(ξ, v)
        𝑃, BasicIdealGas.T(ξ, 𝑃, 𝑣)
    else
        PP(ξ, P), TT(ξ, T)
    end
end

# Specific entropy
function s(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        B::Union{Symbol, Missing} = missing,
        kw...,
    )
    return s(ξ, PT(ξ, P = P, T = T, v = v)..., ismissing(B) ? :MA : B)
end

# Specific Helmholtz energy
function a(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        B::Union{Symbol, Missing} = missing,
        kw...,
    )
    𝑃, 𝑇 = PT(ξ, P = P, T = T, v = v)
    𝑢 = u(ξ.hmod, T, ismissing(B) ? :MA : B)
    𝑠 = s(ξ, 𝑃, 𝑇, ismissing(B) ? :MA : B)
    return 𝑢 - 𝑇 * 𝑠
end

# Specific Gibbs energy
function g(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        B::Union{Symbol, Missing} = missing,
        kw...,
    )
    𝑃, 𝑇 = PT(ξ, P = P, T = T, v = v)
    ℎ = h(ξ.hmod, T, ismissing(B) ? :MA : B)
    𝑠 = s(ξ, 𝑃, 𝑇, ismissing(B) ? :MA : B)
    return ℎ - 𝑇 * 𝑠
end

# Coefficient of volume expansion, β
function β(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        kw...,
    )
    return inv(T(ξ, P = P, T = T, v = v))
end

# Isothermal compressibility, κT
function κT(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        kw...,
    )
    return inv(P(ξ, P = P, T = T, v = v))
end

# TODO: below

__κs = (
    ξ::IdealGas;
    P::Union{Real, PRES, Missing} = missing,
    T::Union{Real, TEMP, Missing} = missing,
    v::Union{Real, VOLU, Missing} = missing,
    B::Symbol = :MA,
) -> begin
    𝑃, 𝑇, 𝑣 = PTv(ξ, P = P, T = T, v = v, B = B)
    inv(𝑃 * ga(ξ, 𝑇))
end

# Base.getproperty - user-facing, oop-style
# -----------------------------------------

import Base: getproperty, propertynames

function Base.getproperty(ξ::IdealGas{ℙ}, sy::Symbol) where {ℙ}
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
            return (; P = missing, T = missing, v = missing, B = missing) -> begin
                getproperty(𝐶, sy)(_T(ξ; P = P, T = T, v = v, B = B))
            end
        elseif sy ∈ props_T_B(𝐶)
            # This allows an 𝑓(T, B) be calc'd from 𝑓(T(P, T, v, B), B)
            # Makes sense only at the IdealGas level (can't fallback to SpecificHeat directly)
            return (; P = missing, T = missing, v = missing, B = missing) -> begin
                getproperty(𝐶, sy)(_T(ξ; P = P, T = T, v = v), ismissing(B) ? :MA : B)
            end
        end
    end
    # OOP-style covenience functions (formerly exported ones)
    return if sy == :P
        (; P = missing, T = missing, v = missing, B = missing) -> _P(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :T
        (; P = missing, T = missing, v = missing, B = missing) -> _T(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :v
        (; P = missing, T = missing, v = missing, B = missing) -> _v(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :vMA
        (; P = missing, T = missing, v = missing) -> _v(ξ; P = P, T = T, v = v, B = :MA)
    elseif sy == :vMO
        (; P = missing, T = missing, v = missing) -> _v(ξ; P = P, T = T, v = v, B = :MO)
    elseif sy == :ρ
        (; P = missing, T = missing, v = missing, B = :MA) -> inv(PTv(ξ; P = P, T = T, v = v, B = B)[3])
    elseif sy == :ρMA
        (; P = missing, T = missing, v = missing) -> inv(PTv(ξ; P = P, T = T, v = v, B = :MA)[3])
    elseif sy == :ρMO
        (; P = missing, T = missing, v = missing) -> inv(PTv(ξ; P = P, T = T, v = v, B = :MO)[3])
    elseif sy == :uMA
        (; P = missing, T = missing, v = missing, B = missing) -> begin
            getproperty(𝐶, :u)(ξ, _T(ξ; P = P, T = T, v = v, B = B), ismissing(B) ? :MA : B)
        end
    elseif sy == :uMO
        (; P = missing, T = missing, v = missing) -> begin
            getproperty(𝐶, :u)(ξ, _T(ξ; P = P, T = T, v = v), B = :MO)
        end
    elseif sy == :s
        (; P = missing, T = missing, v = missing, B = :MA) -> __s(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :sMA
        (; P = missing, T = missing, v = missing) -> __s(ξ; P = P, T = T, v = v, B = :MA)
    elseif sy == :sMO
        (; P = missing, T = missing, v = missing) -> __s(ξ; P = P, T = T, v = v, B = :MO)
    elseif sy == :a
        (; P = missing, T = missing, v = missing, B = :MA) -> __a(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :aMA
        (; P = missing, T = missing, v = missing) -> __a(ξ; P = P, T = T, v = v, B = :MA)
    elseif sy == :aMO
        (; P = missing, T = missing, v = missing) -> __a(ξ; P = P, T = T, v = v, B = :MO)
    elseif sy == :g
        (; P = missing, T = missing, v = missing, B = :MA) -> __g(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :gMA
        (; P = missing, T = missing, v = missing) -> __g(ξ; P = P, T = T, v = v, B = :MA)
    elseif sy == :gMO
        (; P = missing, T = missing, v = missing) -> __g(ξ; P = P, T = T, v = v, B = :MO)
    elseif sy in (:β, :beta)
        (; P = missing, T = missing, v = missing, B = :MA) -> __β(ξ; P = P, T = T, v = v, B = B)
    elseif sy in (:κT, :kappaT)
        (; P = missing, T = missing, v = missing, B = :MA) -> __κT(ξ; P = P, T = T, v = v, B = B)
    elseif sy in (:κs, :kappas)
        (; P = missing, T = missing, v = missing, B = :MA) -> __κs(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :k
        (; P = missing, T = missing, v = missing, B = :MA) -> getproperty(𝐶, :ga)(_T(ξ, T))
    elseif sy == :c
        (; P = missing, T = missing, v = missing, B = :MA) -> begin
            𝑇 = _T(ξ, T)
            γ = getproperty(𝐶, :ga)(𝑇)
            𝑅 = getproperty(𝐶, :RMA)
            uconvert(u"m/s", √(γ * 𝑅 * 𝑇))
        end
    elseif sy in (:μJT, :muJT)
        (; P = missing, T = missing, v = missing, B = :MA) -> zero(ℙ) * u"K/kPa"
    elseif sy in (:μs, :mus)
        (; P = missing, T = missing, v = missing, B = :MA) -> begin
            𝑃, 𝑇, 𝑣 = PTv(ξ, P = P, T = T, v = v, B = B)
            uconvert(u"K/kPa", 𝑣 / getproperty(𝐶, :cp)(𝑇))
        end
    end
end

Base.propertynames(ξ::IdealGas) = (
    :form, :name, :hmod, :𝑃ref,
    :Pref,
    propertynames(getfield(ξ, :hmod))...,
    :P, :T, :v, :vMA, :vMO, :ρ, :ρMA, :ρMO, :s, :sMA, :sMO,
    :a, :aMA, :aMO, :g, :gMA, :gMO, :β, :beta, :κT, :kappaT,
    :κs, :kappas, :k, :c, :μJT, :muJT, :μs, :mus,
)
