# cpModel.jl - General Specific Heat Model

# Structure (type) definition
# ---------------------------

struct SpecificHeat{ℙ <: FLOAT}
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
        @assert(B in (:MA, :MO), "Error: B should be either :MA or :MO")
        wf┆R = (T::Quantity{𝔽, dimension(u"K")} where {𝔽 <: Real}) -> f┆R(ℙ(uconvert(u"K", T).val))
        return new{ℙ}(ID, ℙ ⊚ wf┆R, 𝑀, 𝑇min, 𝑇ref, 𝑇max, 𝑢ref, 𝑠ref, 𝑅)
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
    M = 𝑀 isa MOLW ? uconvert(u"kg/kmol", 𝑀) : 𝑀 * u"kg/kmol"
    Tmin = 𝑇min isa TEMP ? uconvert(u"K", 𝑇min) : 𝑇min * u"K"
    Tref = 𝑇ref isa TEMP ? uconvert(u"K", 𝑇ref) : 𝑇ref * u"K"
    Tmax = 𝑇max isa TEMP ? uconvert(u"K", 𝑇max) : 𝑇max * u"K"
    uref = 𝑢ref isa MASS ? uconvert(u"kJ/kmol", 𝑢ref * 𝑀) : uconvert(u"kJ/kmol", 𝑢ref)
    sref = 𝑠ref isa MASS ? uconvert(u"kJ/kmol/K", 𝑠ref * 𝑀) : uconvert(u"kJ/kmol/K", 𝑠ref)
    R = 𝑅 isa MASS ? uconvert(u"kJ/kmol/K", 𝑅 * 𝑀) : uconvert(u"kJ/kmol/K", 𝑅)
    return SpecificHeat(ID, f┆R, ℙ.((𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅))..., B)
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
    ) where {ℙ <: FLOAT}
    ℙ = promote_type(precof.((𝑀, 𝑇min, 𝑇ref, 𝑇max, 𝑢ref, 𝑠ref))...) # Default 𝑅 left out
    ℙ = ℙ <: FLOAT ? ℙ : Float64
    return SpecificHeat{ℙ}(ID, 𝑓, 𝑀, 𝑇min, 𝑇ref, 𝑇max, 𝑢ref, 𝑠ref, 𝑅)
end

# Conversions
# -----------

import Base: convert

convert(::Type{SpecificHeat{ℙ}}, ξ::SpecificHeat{ℙ}) where {ℙ <: FLOAT} = ξ

function convert(::Type{SpecificHeat{ℙ}}, ξ::SpecificHeat{ℚ}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return SpecificHeat{ℙ}(
        ξ.ID, ξ.f┆R, ξ.𝑀, ξ.𝑇min, ξ.𝑇ref, ξ.𝑇max, ξ.𝑢ref, ξ.𝑠ref, ξ.𝑅
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

function Base.show(io::IO, ::MIME"text/plain", ξ::SpecificHeat{ℙ}) where {ℙ <: FLOAT}
    # rng = "[$(@sprintf("%.*g K", 5, ξ.Tmin)) $(@sprintf("%.*g K", 5, ξ.Tmax))]"
    return print(
        io,
        "$(ξ.ID) cp$(pDeco(ℙ))(T)"
    )
end

# SpecificHeat Helper functions
# -----------------------------

∫┆T(ξ::SpecificHeat, 𝑇::TEMP) = ∫(T -> ξ.f┆R(T) / T, ξ.𝑇ref, 𝑇)
∫┆T(ξ::SpecificHeat, 𝑇::Real) = ∫┆T(ξ, 𝑇 * u"K")

# User-facing functions
# ---------------------

𝗯(ξ::SpecificHeat, 𝑇::TEMP) = begin
    msg = "T = $(@sprintf("%.*g K", 5, 𝑇.val)) out of bounds"
    @assert(ξ.𝑇min <= 𝑇 <= ξ.𝑇max, msg)
end
𝗯(ξ::SpecificHeat, 𝑇::Real) = 𝗯(ξ, 𝑇 * u"K")

import Base: cp

cp┆R(ξ::SpecificHeat{ℙ}, 𝑇::TEMP) where {ℙ <: FLOAT} = (𝗯(ξ, 𝑇); ξ.f┆R(𝑇))
cp┆R(ξ::SpecificHeat{ℙ}, 𝑇::Real) where {ℙ <: FLOAT} = cp┆R(ξ, 𝑇 * u"K")
cv┆R(ξ::SpecificHeat{ℙ}, 𝑇) where {ℙ <: FLOAT} = cp┆R(ξ, 𝑇) - one(ℙ)
ga(ξ::SpecificHeat{ℙ}, 𝑇) where {ℙ <: FLOAT} = (x = cp┆R(ξ, 𝑇); x / (x - one(ℙ)))

function R(ξ::SpecificHeat, B::Symbol = :MA)
    @assert B in (:MA, :MO)
    return B == :MO ? ξ.𝑅 : ξ.𝑅 / ξ.𝑀
end

cp(ξ::SpecificHeat{ℙ}, 𝑇, B = :MA) where {ℙ <: FLOAT} = cp┆R(ξ, 𝑇) * R(ξ, B)
cv(ξ::SpecificHeat{ℙ}, 𝑇, B = :MA) where {ℙ <: FLOAT} = cv┆R(ξ, 𝑇) * R(ξ, B)

∫cp┆R(ξ::SpecificHeat{ℙ}, 𝑇::ℙ) where {ℙ <: FLOAT} = (𝗯(ξ, 𝑇); ∫(ξ.𝑓, ξ.𝑇ref, 𝑇) / ξ.𝑅)
∫cp┆R(ξ::SpecificHeat{ℙ}, 𝑇::Real) where {ℙ <: FLOAT} = ∫cp┆R(ξ, ℙ(𝑇))

∫cv┆R(ξ::SpecificHeat{ℙ}, 𝑇::Real) where {ℙ <: FLOAT} = ∫cp┆R(ξ, 𝑇) - ℙ(𝑇) + ξ.𝑇ref

u┆R(ξ::SpecificHeat{ℙ}, 𝑇::Real) where {ℙ <: FLOAT} = ∫cv┆R(ξ, 𝑇) + ξ.𝑢ref / ξ.𝑅
h┆R(ξ::SpecificHeat{ℙ}, 𝑇::Real) where {ℙ <: FLOAT} = u┆R(ξ, 𝑇) + ℙ(𝑇)
u(ξ::SpecificHeat{ℙ}, 𝑇::Real, B::Symbol = :MA) where {ℙ <: FLOAT} = u┆R(ξ, 𝑇) * R(ξ, B)
h(ξ::SpecificHeat{ℙ}, 𝑇::Real, B::Symbol = :MA) where {ℙ <: FLOAT} = h┆R(ξ, 𝑇) * R(ξ, B)
∫cp┆RT(ξ::SpecificHeat, 𝑇::Real) = (𝗯(ξ, 𝑇); ∫┆T(ξ, 𝑇) / ξ.𝑅)
s0┆R(ξ::SpecificHeat{ℙ}, 𝑇::Real) where {ℙ <: FLOAT} = ∫cp┆RT(ξ, 𝑇) + ξ.𝑠ref / ξ.𝑅
s0(ξ::SpecificHeat{ℙ}, 𝑇::Real, B::Symbol = :MA) where {ℙ <: FLOAT} = s0┆R(ξ, 𝑇) * R(ξ, B)
Pr(ξ::SpecificHeat{ℙ}, 𝑇::Real) where {ℙ <: FLOAT} = exp(∫cp┆RT(ξ, 𝑇))
vr(ξ::SpecificHeat{ℙ}, 𝑇::Real) where {ℙ <: FLOAT} = ℙ(𝑇) / Pr(ξ, 𝑇)

# Base.getproperty
# ----------------

import Base: getproperty, propertynames

function Base.getproperty(ξ::SpecificHeat, sy::Symbol)
    # Raw fields
    if sy in fieldnames(SpecificHeat)
        return getfield(ξ, sy)
    end
    # Convenience accessors/transformers
    if sy == :f
        return getfield(ξ, :𝑓)
    elseif sy == :fMA
        return T -> getfield(ξ, :𝑓)(T) / getfield(ξ, :𝑀)
    end
    # Porcelain accessors (with units)
    if sy == :M
        return getfield(ξ, :𝑀) * u"kg/kmol"
    elseif sy == :R
        return getfield(ξ, :𝑅) * u"kJ/kmol/K"
    elseif sy == :RMA
        return R(ξ, :MA) * u"kJ/kg/K"
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
        return print(join([repr(ξ), string(plt)], "\n"))
    end
    # OOP-style covenience functions (formerly exported ones)
    oop_style_funcs_1 = (
        :cp┆R, :cv┆R, :ga, :R, :∫cp┆R, :∫cv┆R,
        :u┆R, :h┆R, :∫cp┆RT, :s0┆R, :Pr, :vr,
    )
    oop_style_funcs_2 = (
        :cp, :cv, :u, :h, :s0,
    )
    if sy in oop_style_funcs_1
        return (T::Real,) -> eval(sy)(ξ, T)
    elseif sy in oop_style_funcs_2
        return (T::Real, B::Symbol = :MA) -> eval(sy)(ξ, T, B)
    end
end

Base.propertynames(::SpecificHeat) = (
    :ID, :𝑓, :𝑀, :Tmin, :Tmax, :Tref, :uref, :sref, :𝑅,
    :f, :fMA, :M, :R, :RMA, :view,
    :cp┆R, :cv┆R, :ga, :R, :∫cp┆R, :∫cv┆R,
    :u┆R, :h┆R, :∫cp┆RT, :s0┆R, :Pr, :vr,
    :cp, :cv, :u, :h, :s0,
)
