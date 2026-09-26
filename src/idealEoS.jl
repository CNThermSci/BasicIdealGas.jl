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
ρ(ξ::IdealGas, 𝑃::PRES, 𝑇::TEMP, B::Symbol) = inv(v(ξ, 𝑃, 𝑇, B))
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

# Helper PT, PTv functions
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

function PTv(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        B::Union{Symbol, Missing} = missing,
        kw...,
    )
    miss = [ i[1] for i in [(:P, P), (:T, T), (:v, v)] if ismissing(i[2]) ]
    length(miss) <= 1 ||
        throw(ArgumentError(@sprintf("Underspecified state: missing (%s)", join(miss, ", "))))
    return if ismissing(P)
        𝑇 = TT(ξ, T)
        𝑣 = vv(ξ, v)
        BasicIdealGas.P(ξ, 𝑇, 𝑣), 𝑇, 𝑣
    elseif ismissing(T)
        𝑃 = PP(ξ, P)
        𝑣 = vv(ξ, v)
        𝑃, BasicIdealGas.T(ξ, 𝑃, 𝑣), 𝑣
    else
        𝑃 = PP(ξ, P)
        𝑇 = TT(ξ, T)
        return 𝑃, 𝑇, BasicIdealGas.v(ξ, 𝑃, 𝑇, ismissing(B) ? :MA : B)
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
    ℎ = h(ξ.hmod, 𝑇, ismissing(B) ? :MA : B)
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
    return inv(BasicIdealGas.T(ξ, P = P, T = T, v = v))
end

# Isothermal compressibility, κT
function κT(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        kw...,
    )
    return inv(BasicIdealGas.P(ξ, P = P, T = T, v = v))
end

# Isentropic compressibility, κs
function κs(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        kw...,
    )
    𝑃, 𝑇 = PT(ξ, P = P, T = T, v = v)
    return inv(𝑃 * γ(ξ.hmod, 𝑇))
end

# Isentropic expansion exponent, k
function k(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        kw...,
    )
    𝑇 = T(ξ, P = P, T = T, v = v)
    return γ(ξ.hmod, 𝑇)
end

# Adiabatic speed of sound, cs
function cs(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        kw...,
    )
    𝑇 = BasicIdealGas.T(ξ, P = P, T = T, v = v)
    return uconvert(u"m/s", √(γ(ξ.hmod, 𝑇) * R(ξ.hmod, :MA) * 𝑇))
end

# Joule-Thomson coefficient, μJT
μJT(ξ::IdealGas{ℙ}; kw...) where {ℙ} = zero(ℙ) * u"K/kPa"

# Isentropic expansion coefficient, μs
function μs(
        ξ::IdealGas;
        P::Union{PRES, Real, Missing} = missing,
        T::Union{TEMP, Real, Missing} = missing,
        v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
        kw...,
    )
    𝑃, 𝑇, 𝑣 = PTv(ξ, P = P, T = T, v = v, B = :MA)
    return uconvert(u"K/kPa", 𝑣 / cp(ξ.hmod, 𝑇, :MA))
end

# Base.getproperty - user-facing, oop-style
# -----------------------------------------

fields(ξ::IdealGas) = (:form, :name, :hmod, :Pref)
props_UNB(ξ::IdealGas) = (:P, :T, :β, :κT, :κs, :k, :cs, :μJT, :μs)
props_BAS(ξ::IdealGas) = (:v, :ρ, :u, :h, :s, :a, :g)
props_CMP(ξ::IdealGas) = ([Symbol(string(i) * string(j)) for i in props_BAS(ξ) for j in (:MA, :MO)]...,)
props(ξ::IdealGas) = (props_UNB(ξ)..., props_BAS(ξ)..., props_CMP(ξ)...)

function Base.getproperty(ξ::IdealGas{ℙ}, sy::Symbol) where {ℙ}
    # Convenience raw field accessors
    if sy == :form
        return getfield(ξ, :form)
    elseif sy == :name
        return getfield(ξ, :name)
    elseif sy in (:hmod, :heat)
        return getfield(ξ, :hmod)
    elseif sy in (:𝑃ref, :Pref)
        return getfield(ξ, :𝑃ref)
    end
    # Heat model object
    𝐶 = getfield(ξ, :hmod)
    # Pretty-print
    if sy == :view
        return 𝐶.view
    end
    # IdealGas properties
    if sy in props_UNB(ξ)
        return eval(sy)
    elseif sy in props_BAS(ξ)
        return eval(sy)
    elseif sy in props_CMP(ξ)
        fn = Symbol(string(sy)[1:(end - 2)])
        BA = Symbol(last(string(sy), 2))
        return (
            ξ::IdealGas;
            P::Union{PRES, Real, Missing} = missing,
            T::Union{TEMP, Real, Missing} = missing,
            v::Union{VOLU, Tuple{Real, Symbol}, Missing} = missing,
            kw...,
        ) -> eval(fn)(ξ, ismissing(T) ? 𝑇 : T, BA)
    end
    # SpecificHeat model property fallbacks
    if sy in props(𝐶)
        return getproperty(𝐶, sy)
    end
end

Base.propertynames(ξ::IdealGas) = tuple(
    Set([fields(ξ)..., props(ξ)..., fields(ξ.hmod)..., props(ξ.hmod)..., :view])...,
)
