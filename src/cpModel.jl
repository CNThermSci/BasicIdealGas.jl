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
        𝑀::MOLW{𝔸},
        𝑇min::TEMP{𝔹},
        𝑇ref::TEMP{ℂ},
        𝑇max::TEMP{𝔻},
        𝑢ref::ENER{𝔼},
        𝑠ref::ENTR{𝔽},
        𝑅::ENTR{𝔾} = 𝔾(Ru),
    ) where {ℙ <: FLOAT, 𝔸 <: Real, 𝔹 <: Real, ℂ <: Real, 𝔻 <: Real, 𝔼 <: Real, 𝔽 <: Real, 𝔾 <: Real}
    return SpecificHeat(ID, f┆R, ℙ.((𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅))..., B)
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
        𝑅::Real = Ru,
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
        𝑀::Union{Real, MOLW},
        Tmin::Union{Real, TEMP},
        Tref::Union{Real, TEMP},
        Tmax::Union{Real, TEMP},
        uref::ENER,
        sref::ENTR,
        𝑅::Union{Real, ENTR} = Ru,
    ) where {ℙ <: FLOAT}
    uref = uref isa MASS ? kSI(uref) * kSI(𝑀) : kSI(uref)
    sref = sref isa MASS ? kSI(sref) * kSI(𝑀) : kSI(sref)
    𝑅 = 𝑅 isa MASS ? kSI(𝑅) * kSI(𝑀) : kSI(𝑅)
    return SpecificHeat{ℙ}(ID, 𝑓, kSI.((𝑀, Tmin, Tref, Tmax))..., uref, sref, 𝑅, :MO)
end

# Promotion type with unit conversion and stripping / 3 indirections
function SpecificHeat(
        ID::Symbol,
        𝑓::Function,
        𝑀::Union{𝕄, MOLW{𝕄}},
        Tmin::Union{𝕀, TEMP{𝕀}},
        Tref::Union{𝔼, TEMP{𝔼}},
        Tmax::Union{𝔸, TEMP{𝔸}},
        uref::ENER{𝕌},
        sref::ENTR{𝕊},
        𝑅::Union{Real, ENTR} = Ru,
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

function Base.show(io::IO, ::MIME"text/plain", ξ::SpecificHeat{ℙ}) where {ℙ <: FLOAT}
    # rng = "[$(@sprintf("%.*g K", 5, ξ.Tmin)) $(@sprintf("%.*g K", 5, ξ.Tmax))]"
    return print(
        io,
        "$(ξ.ID) cp$(pDeco(ℙ))(T)"
    )
end

# SpecificHeat Helper functions
# -----------------------------

∫┆T(C::SpecificHeat, T::Real) = ∫(T -> C.𝑓(T) / T, C.Tref, T)

# User-facing functions
# ---------------------

𝗯(C::SpecificHeat, T::Real) = begin
    msg = "T = $(@sprintf("%.*g K", 5, T)) out of bounds"
    @assert(C.Tmin <= T <= C.Tmax, msg)
end

import Base: cp

cp┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = (𝗯(C, T); C.𝑓(T) / C.𝑅)
cv┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = cp┆R(C, T) - one(ℙ)
ga(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = (𝗯(C, T); x = C.𝑓(T); x / (x - C.𝑅))

function R(C::SpecificHeat, B::Symbol = :MA)
    @assert B in (:MA, :MO)
    return B == :MO ? C.𝑅 : C.𝑅 / C.𝑀
end

cp(C::SpecificHeat{ℙ}, T::Real, B::Symbol = :MA) where {ℙ <: FLOAT} = cp┆R(C, T) * R(C, B)
cv(C::SpecificHeat{ℙ}, T::Real, B::Symbol = :MA) where {ℙ <: FLOAT} = cv┆R(C, T) * R(C, B)
∫cp┆R(C::SpecificHeat{ℙ}, T::ℙ) where {ℙ <: FLOAT} = (𝗯(C, T); ∫(C.𝑓, C.Tref, T) / C.𝑅)
∫cp┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = ∫cp┆R(C, ℙ(T))
∫cv┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = ∫cp┆R(C, T) - ℙ(T) + C.Tref
u┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = ∫cv┆R(C, T) + C.uref / C.𝑅
h┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = u┆R(C, T) + ℙ(T)
u(C::SpecificHeat{ℙ}, T::Real, B::Symbol = :MA) where {ℙ <: FLOAT} = u┆R(C, T) * R(C, B)
h(C::SpecificHeat{ℙ}, T::Real, B::Symbol = :MA) where {ℙ <: FLOAT} = h┆R(C, T) * R(C, B)
∫cp┆RT(C::SpecificHeat, T::Real) = (𝗯(C, T); ∫┆T(C, T) / C.𝑅)
s0┆R(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = ∫cp┆RT(C, T) + C.sref / C.𝑅
s0(C::SpecificHeat{ℙ}, T::Real, B::Symbol = :MA) where {ℙ <: FLOAT} = s0┆R(C, T) * R(C, B)
Pr(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = exp(∫cp┆RT(C, T))
vr(C::SpecificHeat{ℙ}, T::Real) where {ℙ <: FLOAT} = ℙ(T) / Pr(C, T)

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
