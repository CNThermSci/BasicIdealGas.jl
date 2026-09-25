# cpModel.jl - General Specific Heat Model

# Structure (type) definition
# ---------------------------

struct SpecificHeat{ℙ <: FLOAT} <: ThermModel{ℙ}
    ID::Symbol      # Model ID, as in :cubic, etc...
    f┆R::Function   # Unitless function cp(T)/R: ℙ -> ℙ
    𝑀::Quantity{ℙ, dimension(u"kg/kmol"), typeof(u"kg/kmol")}
    𝑇min::Quantity{ℙ, dimension(u"K"), typeof(u"K")}
    𝑇ref::Quantity{ℙ, dimension(u"K"), typeof(u"K")}
    𝑇max::Quantity{ℙ, dimension(u"K"), typeof(u"K")}
    𝑢ref::Quantity{ℙ, dimension(u"kJ/kmol"), typeof(u"kJ/kmol")}
    𝑠ref::Quantity{ℙ, dimension(u"kJ/kmol/K"), typeof(u"kJ/kmol/K")}
    𝑅::Quantity{ℙ, dimension(u"kJ/kmol/K"), typeof(u"kJ/kmol/K")}
    # Internal, validating constructor
    SpecificHeat(
        ID::Symbol,
        f┆R::Function,
        𝑀::Quantity{ℙ, dimension(u"kg/kmol"), typeof(u"kg/kmol")},
        𝑇min::Quantity{ℙ, dimension(u"K"), typeof(u"K")},
        𝑇ref::Quantity{ℙ, dimension(u"K"), typeof(u"K")},
        𝑇max::Quantity{ℙ, dimension(u"K"), typeof(u"K")},
        𝑢ref::Quantity{ℙ, dimension(u"kJ/kmol"), typeof(u"kJ/kmol")},
        𝑠ref::Quantity{ℙ, dimension(u"kJ/kmol/K"), typeof(u"kJ/kmol/K")},
        𝑅::Quantity{ℙ, dimension(u"kJ/kmol/K"), typeof(u"kJ/kmol/K")} = ℙ(Ru),
    ) where {ℙ <: FLOAT} = begin
        @assert(ID != Symbol(""), "Error: Empty model ID")
        @assert(𝑀 > zero(ℙ) * u"kg/kmol", "Error: M <= 0 kg/kmol")
        @assert(zero(ℙ) * u"K" <= 𝑇min <= 𝑇ref < 𝑇max, "Error: Temperature values")
        @assert(𝑅 > zero(ℙ) * u"kJ/kmol/K", "Error: 𝑅 <= 0 kJ/kmol/K")
        wf┆R = (T::Quantity{𝔽, dimension(u"K")} where {𝔽 <: Real}) -> (ℙ ⊚ f┆R)(ℙ(uconvert(u"K", T).val))
        return new{ℙ}(ID, wf┆R, 𝑀, 𝑇min, 𝑇ref, 𝑇max, 𝑢ref, 𝑠ref, 𝑅)
    end
end

# External constructors
# ---------------------

# Set precision conversion / 1 indirection
function SpecificHeat{ℙ}(
        ID::Symbol,
        f┆R::Function,
        𝑀::Union{Real, MOLW},
        𝑇min::Union{Real, TEMP},
        𝑇ref::Union{Real, TEMP},
        𝑇max::Union{Real, TEMP},
        𝑢ref::ENER,
        𝑠ref::ENTR,
        𝑅::ENTR = ℙ(Ru),
    ) where {ℙ <: FLOAT}
    M = MM(𝑀)
    Tmin = TT(𝑇min)
    Tref = TT(𝑇ref)
    Tmax = TT(𝑇max)
    uref = ee(𝑢ref)
    sref = ss(𝑠ref)
    R = ss(𝑅)
    return SpecificHeat(ID, f┆R, ℙ.((M, Tmin, Tref, Tmax, uref, sref, 𝑅))...)
end

# Promotion type conversion / 2 indirections
function SpecificHeat(
        ID::Symbol,
        f┆R::Function,
        𝑀::Union{Real, MOLW},
        𝑇min::Union{Real, TEMP},
        𝑇ref::Union{Real, TEMP},
        𝑇max::Union{Real, TEMP},
        𝑢ref::ENER,
        𝑠ref::ENTR,
        𝑅::ENTR = Ru,
    )
    ℙ = promote_type(precof.((𝑀, 𝑇min, 𝑇ref, 𝑇max, 𝑢ref, 𝑠ref))...) # Default 𝑅 left out
    ℙ = ℙ <: FLOAT ? ℙ : Float64
    return SpecificHeat{ℙ}(ID, f┆R, 𝑀, 𝑇min, 𝑇ref, 𝑇max, 𝑢ref, 𝑠ref, 𝑅)
end

# Conversions
# -----------

import Base: convert

convert(::Type{SpecificHeat{ℙ}}, ξ::SpecificHeat{ℙ}) where {ℙ <: FLOAT} = ξ

function convert(::Type{SpecificHeat{ℙ}}, ξ::SpecificHeat{ℚ}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return if ξ.f┆R.f┆R isa ComposedFunction
        SpecificHeat{ℙ}(ξ.ID, ξ.f┆R.f┆R.inner, ξ.𝑀, ξ.𝑇min, ξ.𝑇ref, ξ.𝑇max, ξ.𝑢ref, ξ.𝑠ref, ξ.𝑅)
    else
        SpecificHeat{ℙ}(ξ.ID, ξ.f┆R.f┆R, ξ.𝑀, ξ.𝑇min, ξ.𝑇ref, ξ.𝑇max, ξ.𝑢ref, ξ.𝑠ref, ξ.𝑅)
    end
end

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

pretty(ξ::SpecificHeat{ℙ}) where {ℙ <: FLOAT} = "$(ξ.ID) cp$(pDeco(ℙ))(T)"

function Base.show(io::IO, ::MIME"text/plain", ξ::SpecificHeat{ℙ}) where {ℙ <: FLOAT}
    return print(io, pretty(ξ))
end

# Internal Property Calculations: positional, dispatched, no default args
# -----------------------------------------------------------------------

import Base: cp

# Temperature args types:
#   𝑇::TEMP                 => dispatched fallback
#   T::Real                 => add units and fallback
#   θ::Union{Real, TEMP}    => fallback

# Ancillary: bounds
𝗯(ξ::SpecificHeat{ℙ}, 𝑇::TEMP) where {ℙ <: FLOAT} = begin
    T = ℙ(𝑇)
    msg = "T = $(@sprintf("%.*g K", 5, T.val)) out of bounds"
    @assert(ξ.𝑇min <= T <= ξ.𝑇max, msg)
end
𝗯(ξ::SpecificHeat, T::Real) = 𝗯(ξ, TT(T))

# Primitives
cp┆R(ξ::SpecificHeat{ℙ}, 𝑇::TEMP) where {ℙ <: FLOAT} = (𝗯(ξ, 𝑇); ξ.f┆R(ℙ(𝑇)))
cp┆R(ξ::SpecificHeat, T::Real) = cp┆R(ξ, TT(T))
cv┆R(ξ::SpecificHeat{ℙ}, θ::Union{Real, TEMP}) where {ℙ <: FLOAT} = cp┆R(ξ, θ) - one(ℙ)

# Base-selected R
function R(ξ::SpecificHeat, B::Symbol)
    @assert B in (:MA, :MO)
    return B == :MO ? ξ.𝑅 : ξ.𝑅 / ξ.𝑀
end
R(ξ::SpecificHeat, B::MASS) = R(ξ, :MA)
R(ξ::SpecificHeat, B::MOLR) = R(ξ, :MO)

# Derived, base-independent properties
γ(ξ::SpecificHeat, θ::Union{Real, TEMP}) = cp┆R(ξ, θ) / cv┆R(ξ, θ)
∫cp┆R(ξ::SpecificHeat{ℙ}, 𝑇::TEMP) where {ℙ <: FLOAT} = (𝗯(ξ, 𝑇); ∫(ξ.f┆R, ξ.𝑇ref, ℙ(𝑇)))
∫cp┆R(ξ::SpecificHeat, T::Real) = ∫cp┆R(ξ, TT(T))
∫cv┆R(ξ::SpecificHeat{ℙ}, 𝑇::TEMP) where {ℙ <: FLOAT} = ∫cp┆R(ξ, 𝑇) - ℙ(𝑇) + ξ.𝑇ref
∫cv┆R(ξ::SpecificHeat, T::Real) = ∫cv┆R(ξ, TT(T))
u┆R(ξ::SpecificHeat, θ::Union{Real, TEMP}) = ∫cv┆R(ξ, θ) + ξ.𝑢ref / ξ.𝑅
h┆R(ξ::SpecificHeat{ℙ}, θ::Union{Real, TEMP}) where {ℙ <: FLOAT} = u┆R(ξ, θ) + ℙ(θ)
∫cp┆RT(ξ::SpecificHeat{ℙ}, 𝑇::TEMP) where {ℙ <: FLOAT} = (𝗯(ξ, 𝑇); ∫(T -> ξ.f┆R(T) / T, ξ.𝑇ref, ℙ(𝑇)))
∫cp┆RT(ξ::SpecificHeat, T::Real) = ∫cp┆RT(ξ, TT(T))
s0┆R(ξ::SpecificHeat, θ::Union{Real, TEMP}) = ∫cp┆RT(ξ, θ) + ξ.𝑠ref / ξ.𝑅
Pr(ξ::SpecificHeat, θ::Union{Real, TEMP}) = exp(∫cp┆RT(ξ, θ))
vr(ξ::SpecificHeat{ℙ}, θ::Union{Real, TEMP}) where {ℙ <: FLOAT} = ℙ(θ) / Pr(ξ, θ)

# Derived, based properties
cp(ξ::SpecificHeat, θ::Union{Real, TEMP}, B::Union{Symbol, MASS, MOLR}) = cp┆R(ξ, θ) * R(ξ, B)
cv(ξ::SpecificHeat, θ::Union{Real, TEMP}, B::Union{Symbol, MASS, MOLR}) = cv┆R(ξ, θ) * R(ξ, B)
u(ξ::SpecificHeat, θ::Union{Real, TEMP}, B::Union{Symbol, MASS, MOLR}) = u┆R(ξ, θ) * R(ξ, B)
h(ξ::SpecificHeat, θ::Union{Real, TEMP}, B::Union{Symbol, MASS, MOLR}) = h┆R(ξ, θ) * R(ξ, B)
s0(ξ::SpecificHeat, θ::Union{Real, TEMP}, B::Union{Symbol, MASS, MOLR}) = s0┆R(ξ, θ) * R(ξ, B)

# Base.getproperty - user-facing, oop-style
# -----------------------------------------

fields(ξ::SpecificHeat) = (:ID, :f, :M, :R, :RMO, :RMA, :Tmin, :Tref, :Tmax, :uref, :sref)
props_T(ξ::SpecificHeat) = (:γ, :Pr, :vr)
props_T_B(ξ::SpecificHeat) = (:cp, :cv, :u, :h, :s0)

import Base: getproperty, propertynames

function Base.getproperty(ξ::SpecificHeat, sy::Symbol)
    # Convenience accessors/transformers
    if sy == :ID
        return getfield(ξ, :ID)
    elseif sy in (:f┆R, :f)
        return getfield(ξ, :f┆R)
    elseif sy in (:𝑀, :M)
        return getfield(ξ, :𝑀)
    elseif sy in (:𝑅, :R, :RMO)
        return getfield(ξ, :𝑅)
    elseif sy == :RMA
        return R(ξ, :MA)
    elseif sy in (:𝑇min, :Tmin)
        return getfield(ξ, :𝑇min)
    elseif sy in (:𝑇ref, :Tref)
        return getfield(ξ, :𝑇ref)
    elseif sy in (:𝑇max, :Tmax)
        return getfield(ξ, :𝑇max)
    elseif sy in (:𝑢ref, :uref)
        return getfield(ξ, :𝑢ref)
    elseif sy in (:𝑠ref, :sref)
        return getfield(ξ, :𝑠ref)
    end
    # Pretty print
    if sy == :view
        xmin, xmax = getfield(ξ, :𝑇min), getfield(ξ, :𝑇max)
        x = range(xmin, stop = xmax, length = 33)
        y = map(T -> cp(ξ, T, :MA), x)
        plt = lineplot(
            x, y, xlabel = "T", ylabel = "cp (T)",
            xlim = (xmin, xmax), width = 32, height = 6,
            border = :ascii, color = :white, compact_labels = true,
        )
        return println(join([pretty(ξ), string(plt)], "\n"))
    end
    # OOP-style covenience functions (formerly exported ones)
    if sy in props_T(ξ)
        return (
            𝑇::Union{Real, TEMP, Missing} = missing;
            T::Union{Real, TEMP, Missing} = missing,
            kw...,
        ) -> eval(sy)(ξ, ismissing(T) ? 𝑇 : T)
    elseif sy in props_T_B(ξ)
        return (
            𝑇::Union{Real, TEMP, Missing} = missing,
            𝐵::Symbol = :MA;
            T::Union{Real, TEMP, Missing} = missing,
            B::Union{Symbol, Missing} = missing,
            kw...,
        ) -> eval(sy)(ξ, ismissing(T) ? 𝑇 : T, ismissing(B) ? 𝐵 : B)
    elseif sy in [ Symbol(string(i) * string(j)) for i in props_T_B(ξ) for j in (:MA, :MO) ]
        fn = Symbol(string(sy)[1:end-2])
        BA = Symbol(last(string(sy), 2))
        return (
            𝑇::Union{Real, TEMP, Missing} = missing;
            T::Union{Real, TEMP, Missing} = missing,
            kw...,
        ) -> eval(fn)(ξ, ismissing(T) ? 𝑇 : T, BA)
    end
end

Base.propertynames(ξ::SpecificHeat) = (fields(ξ)..., props_T(ξ)..., props_T_B(ξ)..., :view)
