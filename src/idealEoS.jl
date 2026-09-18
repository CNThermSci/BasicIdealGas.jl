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

# Internal positional EoS functions
function _P(ξ::IdealGas{ℙ}, 𝑇::TEMP, 𝑣::VOLU, B::Symbol = :MA) where {ℙ}
    return uconvert(u"kPa", R(ξ, B) * ℙ(𝑇 / 𝑣))
end

function _T(ξ::IdealGas{ℙ}, 𝑃::PRES, 𝑣::VOLU, B::Symbol = :MA) where {ℙ}
    return uconvert(u"K", ℙ(𝑃 * 𝑣) / R(ξ, B))
end

function _v(ξ::IdealGas{ℙ}, 𝑃::PRES, 𝑇::TEMP, B::Symbol = :MA) where {ℙ}
    UNIT = B == :MA ? u"m^3/kg" : u"m^3/kmol"
    return uconvert(UNIT, R(ξ, B) * ℙ(𝑇 / 𝑃))
end

function _ρ(ξ::IdealGas{ℙ}, 𝑃::PRES, 𝑇::TEMP, B::Symbol = :MA) where {ℙ}
    return inv(_v(ξ, 𝑃, 𝑇, B))
end

# Internal positional entropy function
function _s(ξ::IdealGas{ℙ}, 𝑃::PRES, 𝑇::TEMP, B::Symbol = :MA) where {ℙ}
    return s0(ξ, 𝑇, B) - R(ξ, B) * log(ℙ(𝑃) / ξ.𝑃ref)
end

# Internal keyworded generic function
kwP = (
    ξ::IdealGas;
    P::Union{Real, PRES, Missing} = missing,
    T::Union{Real, TEMP, Missing} = missing,
    v::Union{Real, VOLU, Missing} = missing,
    B::Symbol = :MA,
) -> begin
    if !ismissing(P)
        return P isa PRES ? uconvert(u"kPa", P) : P * u"kPa"
    end
    count(x -> ismissing(x), (T, v)) > 0 && throw(
        ArgumentError("Unspecified state: (T = $(T), v = $(v))")
    )
    𝑇 = T isa TEMP ? T : T * u"K"
    𝑣 = v isa VOLU ? v : v * (B == :MA ? u"m^3/kg" : u"m^3/kmol")
    _P(ξ, 𝑇, 𝑣, B)
end

kwT = (
    ξ::IdealGas;
    P::Union{Real, PRES, Missing} = missing,
    T::Union{Real, TEMP, Missing} = missing,
    v::Union{Real, VOLU, Missing} = missing,
    B::Symbol = :MA,
) -> begin
    if !ismissing(T)
        return T isa TEMP ? uconvert(u"K", T) : T * u"K"
    end
    count(x -> ismissing(x), (P, v)) > 0 && throw(
        ArgumentError("Unspecified state: (P = $(P), v = $(v))")
    )
    𝑃 = P isa PRES ? P : P * u"kPa"
    𝑣 = v isa VOLU ? v : v * (B == :MA ? u"m^3/kg" : u"m^3/kmol")
    _T(ξ, 𝑃, 𝑣, B)
end

kwv = (
    ξ::IdealGas;
    P::Union{Real, PRES, Missing} = missing,
    T::Union{Real, TEMP, Missing} = missing,
    v::Union{Real, VOLU, Missing} = missing,
    B::Symbol = :MA,
) -> begin
    if !ismissing(v)
        UNIT = B == :MA ? u"m^3/kg" : u"m^3/kmol"
        return v isa VOLU ? uconvert(UNIT, v) : v * UNIT
    end
    count(x -> ismissing(x), (P, T)) > 0 && throw(
        ArgumentError("Unspecified state: (P = $(P), T = $(T))")
    )
    𝑃 = P isa PRES ? P : P * u"kPa"
    𝑇 = T isa TEMP ? T : T * u"K"
    _v(ξ, 𝑃, 𝑇, B)
end

kwρ = (
    ξ::IdealGas;
    P::Union{Real, PRES, Missing} = missing,
    T::Union{Real, TEMP, Missing} = missing,
    v::Union{Real, VOLU, Missing} = missing,
    B::Symbol = :MA,
) -> inv(kwv(ξ; P = P, T = T, B = B))

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
    if sy in propertynames(getfield(ξ, :hmod))
        return getproperty(getfield(ξ, :hmod), sy)
    end
    # OOP-style covenience functions (formerly exported ones)
    if sy == :P
        return (; P = missing, T = missing, v = missing, B = :MA) -> kwP(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :T
        return (; P = missing, T = missing, v = missing, B = :MA) -> kwT(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :v
        return (; P = missing, T = missing, v = missing, B = :MA) -> kwv(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :ρ
        return (; P = missing, T = missing, v = missing, B = :MA) -> kwρ(ξ; P = P, T = T, v = v, B = B)
    elseif sy == :s
        return (
            ;
            P::Union{Real, PRES, Missing} = missing,
            T::Union{Real, TEMP, Missing} = missing,
            v::Union{Real, VOLU, Missing} = missing,
            B::Symbol = :MA,
        ) -> begin
            @assert(
                count(x -> ismissing(x), (P, T, v)) == 1,
                "exactly two P-T-v state functions must be specified!"
            )
            return if ismissing(P)
                𝑇 = T isa TEMP ? T : T * u"K"
                𝑣 = v isa VOLU ? v : v * (B == :MA ? u"m^3/kg" : u"m^3/kmol")
                _s(ξ, _P(ξ, 𝑇, 𝑣, B), 𝑇, B)
            elseif ismissing(T)
                𝑃 = P isa PRES ? P : P * u"kPa"
                𝑣 = v isa VOLU ? v : v * (B == :MA ? u"m^3/kg" : u"m^3/kmol")
                _s(ξ, 𝑃, _T(ξ, 𝑃, 𝑣, B), B)
            else
                𝑇 = T isa TEMP ? T : T * u"K"
                𝑃 = P isa PRES ? P : P * u"kPa"
                _s(ξ, 𝑃, 𝑇, B)
            end
        end
    end
end

Base.propertynames(ξ::IdealGas) = (
    :form, :name, :hmod, :𝑃ref,
    :Pref,
    propertynames(getfield(ξ, :hmod))...,
    :P, :T, :v, :ρ, :s,
)
