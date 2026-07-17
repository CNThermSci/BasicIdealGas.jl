# idealState.jl - Ideal Gas States

# Structure (type) definition
# ---------------------------

struct IdealState{ℙ <: FLOAT}
    𝐺::IdealGas{ℙ}
    𝑃::ℙ                # kPa
    𝑇::ℙ                # K
    # Internal, validating constructors
    function IdealState(
            G::IdealGas{ℙ},
            P::ℙ,
            T::ℙ,
        ) where {ℙ <: FLOAT}
        @assert(G.hmod.Tmin <= T <= G.hmod.Tmax, "T = $(T) out of range for $(G.hmod)")
        return new{ℙ}(G, P, T)
    end
end

# External constructors
# ---------------------

# Set type conversion / 1 indirection
IdealState{ℙ}(
    G::IdealGas,
    P::Real,
    T::Real,
) where {ℙ} = IdealState(ℙ.((G, P, T))...)

# Ideal gas model type conversion / 2 indirections
IdealState(
    G::IdealGas{ℙ},
    P::Real,
    T::Real,
) where {ℙ} = IdealState{ℙ}(G, P, T)

# Set type with unit conversion and stripping / 2 indirections
function IdealState{ℙ}(
        G::IdealGas,
        P::Quantity{<:Real, dimension(u"kPa")},
        T::Quantity{<:Real, dimension(u"K")},
    ) where {ℙ <: FLOAT}
    return IdealState{ℙ}(
        G,
        P isa Quantity ? uconvert(u"kPa", P).val : P,
        T isa Quantity ? uconvert(u"K", T).val : T,
    )
end

# Heat model type with unit conversion and stripping / 3 indirections
IdealState(
    G::IdealGas{ℙ},
    P::Quantity{<:Real, dimension(u"kPa")},
    T::Quantity{<:Real, dimension(u"K")},
) where {ℙ <: FLOAT} = IdealState{ℙ}(G, P, T)

# Conversions
# -----------

import Base: convert

convert(::Type{IdealState{ℙ}}, ξ::IdealState{ℙ}) where {ℙ <: FLOAT} = ξ

function convert(::Type{IdealState{ℙ}}, ξ::IdealState{ℚ}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return IdealState{ℙ}(ξ.𝐺, ξ.𝑃, ξ.𝑇)
end

import Base: Float16, Float32, Float64

Float16(ξ::IdealState) = convert(IdealState{Float16}, ξ)
Float32(ξ::IdealState) = convert(IdealState{Float32}, ξ)
Float64(ξ::IdealState) = convert(IdealState{Float64}, ξ)

# Promotions
# ----------

import Base: promote_rule

function promote_rule(::Type{IdealState{ℙ}}, ::Type{IdealState{ℚ}}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return IdealState{promote_type(ℙ, ℚ)}
end

# Functor
# -------

function (ξ::IdealState{ℙ})(
    ;
    P::Union{Missing, Real, Quantity{<:Real, dimension(u"kPa")}} = missing,
    T::Union{Missing, Real, Quantity{<:Real, dimension(u"K")}} = missing,
) where {ℙ}
    return if count(x -> !isa(x, Missing), (P, T)) == 0
        # named tuple variant
        (
            M = ξ.M,
            RMO = ξ.RMO,
            R = ξ.R,
            P = ξ.P,
            T = ξ.T,
            v = ξ.v,
            vMO = ξ.vMO,
            u = ξ.u,
            uMO = ξ.uMO,
            h = ξ.h,
            hMO = ξ.hMO,
            s = ξ.s,
            sMO = ξ.sMO,
        )
    else
        # copy-edit variant
        IdealState{ℙ}(
            ξ.𝐺,
            P isa Missing ? ξ.𝑃 : P isa Quantity ? uconvert(u"kPa", P).val : P,
            T isa Missing ? ξ.𝑇 : T isa Quantity ? uconvert(u"K", T).val : T,
        )
    end
end

# Export
# ------

export IdealState

# Base.getproperty
# ----------------

import Base: getproperty, propertynames

function Base.getproperty(ξ::IdealState, sy::Symbol)
    # Raw fields
    if sy in fieldnames(IdealState)
        return getfield(ξ, sy)
    elseif sy in fieldnames(IdealGas)
        return getfield(getfield(ξ, :𝐺), sy)
    elseif sy in fieldnames(SpecificHeat)
        return getfield(getfield(getfield(ξ, :𝐺), :hmod), sy)
    end
    # User-facing state function accessors (with units)
    GAS, P, T = map(sy -> getfield(ξ, sy), (:𝐺, :𝑃, :𝑇))
    MOD = getfield(GAS, :hmod)
    if sy == :gas
        return GAS
    elseif sy == :mod
        return MOD
    elseif sy == :P
        return getfield(ξ, :𝑃) * u"kPa"
    elseif sy == :T
        return getfield(ξ, :𝑇) * u"K"
    elseif sy in (:γ, :ga)
        return ga(MOD, T)
    elseif sy == :v
        return _v(GAS, P, T, :MA) * u"m^3/kg"
    elseif sy == :vMO
        return _v(GAS, P, T, :MO) * u"m^3/kmol"
    elseif sy == :ρ
        return _ρ(GAS, P, T, :MA) * u"kg/m^3"
    elseif sy == :ρMO
        return _ρ(GAS, P, T, :MO) * u"kmol/m^3"
    elseif sy == :cp
        return cp(MOD, T, :MA) * u"kJ/kg/K"
    elseif sy == :cpMO
        return cp(MOD, T, :MO) * u"kJ/kmol/K"
    elseif sy == :cv
        return cv(MOD, T, :MA) * u"kJ/kg/K"
    elseif sy == :cvMO
        return cv(MOD, T, :MO) * u"kJ/kmol/K"
    elseif sy == :u
        return u(MOD, T, :MA) * u"kJ/kg"
    elseif sy == :uMO
        return u(MOD, T, :MO) * u"kJ/kmol"
    elseif sy == :h
        return h(MOD, T, :MA) * u"kJ/kg"
    elseif sy == :hMO
        return h(MOD, T, :MO) * u"kJ/kmol"
    elseif sy == :s0
        return s0(MOD, T, :MA) * u"kJ/kg/K"
    elseif sy == :s0MO
        return s0(MOD, T, :MO) * u"kJ/kmol/K"
    elseif sy == :s
        return _s(GAS, P, T, :MA) * u"kJ/kg/K"
    elseif sy == :sMO
        return _s(GAS, P, T, :MO) * u"kJ/kmol/K"
    elseif sy == :Pr
        return Pr(MOD, T)
    elseif sy == :vr
        return vr(MOD, T)
    end
end

Base.propertynames(ξ::IdealState) = (
    fieldnames(IdealState)...,
    fieldnames(IdealGas)...,
    fieldnames(SpecificHeat)...,
    :gas, :mod, :P, :T, :γ, :ga, :v, :vMO, :ρ, :ρMO, :cp, :cpMO,
    :cv, :cvMO, :u, :uMO, :h, :hMO, :s0, :s0MO, :s, :sMO, :Pr, :vr,
)

# User-facing functions
# ---------------------

function Base.show(io::IO, ::MIME"text/plain", st::IdealState{ℙ}) where {ℙ <: FLOAT}
    return print(
        io,
        repr(MIME"text/plain"(), st.gas),
        " @($(@sprintf("%.*g", 5, st.𝑃)) kPa, $(@sprintf("%.*g", 5, st.𝑇)) K)"
    )
end
