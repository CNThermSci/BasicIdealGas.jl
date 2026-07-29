# idealProcs.jl - Ideal Gas Processes

# Ancillary type definitions
# --------------------------

struct Interact{ℙ <: FLOAT}
    𝑞::ℙ            # specific heat interaction, kJ/kg
    𝑤::ℙ            # specific work interaction, kJ/kg
    function Interact(
            q::ℙ,
            w::ℙ,
        ) where {ℙ <: FLOAT}
        return new{ℙ}(q, w)
    end
end

Interact{ℙ}(
    q::Real,
    w::Real,
) where {ℙ} = Interact(ℙ.((q, w))...)

import Base: convert, promote_rule
import Base: Float16, Float32, Float64
import Base: getproperty, propertynames
export Interact

convert(::Type{Interact{ℙ}}, ξ::Interact{ℙ}) where {ℙ <: FLOAT} = ξ
function convert(::Type{Interact{ℙ}}, ξ::Interact{ℚ}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return Interact{ℙ}(ξ.𝑞, ξ.𝑤)
end
Float16(ξ::Interact) = convert(Interact{Float16}, ξ)
Float32(ξ::Interact) = convert(Interact{Float32}, ξ)
Float64(ξ::Interact) = convert(Interact{Float64}, ξ)
function promote_rule(::Type{Interact{ℙ}}, ::Type{Interact{ℚ}}) where {ℙ <: FLOAT, ℚ <: FLOAT}
    return Interact{promote_type(ℙ, ℚ)}
end
function (ξ::Interact{ℙ})(
        ;
        q::Union{Missing, Real, ENER} = missing,
        w::Union{Missing, Real, ENER} = missing,
    ) where {ℙ}
    return if count(x -> !isa(x, Missing), (q, w)) == 0
        pairs((q = ξ.𝑞, w = ξ.𝑤))
    else
        Interact{ℙ}(
            q isa Missing ? ξ.𝑞 : q,
            w isa Missing ? ξ.𝑤 : w,
        )
    end
end
function Base.getproperty(ξ::Interact, sy::Symbol)
    # Raw fields
    if sy in fieldnames(Interact)
        return getfield(ξ, sy)
    end
    # User-facing state function accessors (with units)
    if sy == :q
        return getfield(ξ, :𝑞) * u"kJ/kg"
    elseif sy == :w
        return getfield(ξ, :𝑤) * u"kJ/kg"
    end
end
Base.propertynames(ξ::Interact) = (fieldnames(Interact)..., :q, :w)

# Structure (type) definition
# ---------------------------

struct IdealProc{ℙ <: FLOAT}
    𝐺::IdealGas{ℙ}              # gas
    𝑖::PropPair{ℙ}              # process initial state
    𝑓::PropPair{ℙ}              # process final state
    intr::Interact{ℙ}           # process interactions
    path::Vector{PropPair{ℙ}}   # process path
    proc::Symbol                # process type
    # Internal, validating constructors
    function IdealProc(
            G::IdealGas{ℙ},
            i::PropPair{ℙ},
            f::PropPair{ℙ},
            intr::Interact{ℙ},
            path::Vector{PropPair{ℙ}},
            proc::Symbol,
        ) where {ℙ <: FLOAT}
        @assert(G.hmod.Tmin <= i.𝑇 <= G.hmod.Tmax, "T = $(i.𝑇) out of range for $(G.hmod)")
        @assert(G.hmod.Tmin <= f.𝑇 <= G.hmod.Tmax, "T = $(f.𝑇) out of range for $(G.hmod)")
        @assert(proc in (:P, :T, :v, :u, :h, :s, :poly, :other), "Invalid process type: '$(proc)'")
        return new{ℙ}(G, i, f, intr, path, proc)
    end
end
# Isobaric processes
# ------------------

isoP_T(ξ::IdealState, T::Union{Missing, Real, TEMP}) = ξ(T = T)

isoP_v(ξ::IdealState, v::Real, B::Symbol) = ξ(T = _T(ξ.gas, ξ.𝑃, v, B))
isoP_v(ξ::IdealState, v::VOLU) = isoP_v(ξ, kSI(v), v isa MOLR ? :MO : :MA)

isoP_u(ξ::IdealState, u::Real, B::Symbol) =
    ξ(
    T = find_zero(
        B == :MO ? T -> kSI(ξ(T = T).uMO) - u : T -> kSI(ξ(T = T).u) - u,
        (ξ.Tmin, ξ.Tmax), Bisection()
    )
)
isoP_u(ξ::IdealState, u::ENER) = isoP_u(ξ, kSI(u), u isa MOLR ? :MO : :MA)

isoP_h(ξ::IdealState, h::Real, B::Symbol) =
    ξ(
    T = find_zero(
        B == :MO ? T -> kSI(ξ(T = T).hMO) - h : T -> kSI(ξ(T = T).h) - h,
        (ξ.Tmin, ξ.Tmax), Bisection()
    )
)
isoP_h(ξ::IdealState, h::ENER) = isoP_h(ξ, kSI(h), h isa MOLR ? :MO : :MA)

isoP_s(ξ::IdealState, s::Real, B::Symbol) =
    ξ(
    T = find_zero(
        B == :MO ? T -> kSI(ξ(T = T).sMO) - s : T -> kSI(ξ(T = T).s) - s,
        (ξ.Tmin, ξ.Tmax), Bisection()
    )
)
isoP_s(ξ::IdealState, s::ENTR) = isoP_s(ξ, kSI(s), s isa MOLR ? :MO : :MA)

function isoP(
        ξ::IdealState;
        T::Union{Missing, ℙ, TEMP{ℙ}} where {ℙ <: Real} = missing,
        v::Union{Missing, Tuple{ℙ, Symbol}, VOLU{ℙ}} where {ℙ <: Real} = missing,
        u::Union{Missing, Tuple{ℙ, Symbol}, ENER{ℙ}} where {ℙ <: Real} = missing,
        h::Union{Missing, Tuple{ℙ, Symbol}, ENER{ℙ}} where {ℙ <: Real} = missing,
        s::Union{Missing, Tuple{ℙ, Symbol}, ENTR{ℙ}} where {ℙ <: Real} = missing,
    )
    @assert(
        count(x -> !isa(x, Missing), (T, v, u, h, s)) == 1,
        "exactly one end-state function must be specified!"
    )
    return if !ismissing(T)
        isoP_T(ξ, T)
    elseif !ismissing(v)
        v isa Tuple ? isoP_v(ξ, v...) : isoP_v(ξ, v)
    elseif !ismissing(u)
        u isa Tuple ? isoP_u(ξ, u...) : isoP_u(ξ, u)
    elseif !ismissing(h)
        h isa Tuple ? isoP_h(ξ, h...) : isoP_h(ξ, h)
    elseif !ismissing(s)
        s isa Tuple ? isoP_s(ξ, s...) : isoP_s(ξ, s)
    end
end

export isoP

# Isothermal processes
# --------------------

isoT_P(ξ::IdealState, P::Union{Missing, Real, PRES}) = ξ(P = P)

isoT_v(ξ::IdealState, v::Real, B::Symbol) = ξ(P = _P(ξ.gas, ξ.𝑇, v, B))
isoT_v(ξ::IdealState, v::VOLU) = isoT_v(ξ, kSI(v), v isa MOLR ? :MO : :MA)

isoT_s(ξ::IdealState, s::Real, B::Symbol) =
    ξ(P = ξ.𝑃 * exp(B == :MO ? (kSI(ξ.sMO) - s) / kSI(ξ.R) : (kSI(ξ.s) - s) / kSI(ξ.RMA)))
isoT_s(ξ::IdealState, s::ENTR) = isoT_s(ξ, kSI(s), s isa MOLR ? :MO : :MA)

function isoT(
        ξ::IdealState;
        P::Union{Missing, ℙ, PRES{ℙ}} where {ℙ <: Real} = missing,
        v::Union{Missing, Tuple{ℙ, Symbol}, VOLU{ℙ}} where {ℙ <: Real} = missing,
        s::Union{Missing, Tuple{ℙ, Symbol}, ENTR{ℙ}} where {ℙ <: Real} = missing,
    )
    @assert(
        count(x -> !isa(x, Missing), (P, v, s)) == 1,
        "exactly one end-state function must be specified!"
    )
    return if !ismissing(P)
        isoT_P(ξ, P)
    elseif !ismissing(v)
        v isa Tuple ? isoT_v(ξ, v...) : isoT_v(ξ, v)
    elseif !ismissing(s)
        s isa Tuple ? isoT_s(ξ, s...) : isoT_s(ξ, s)
    end
end

isou = isoT
isoh = isoT

export isoT, isou, isoh

# Isochoric processes
# -------------------

isov_P(ξ::IdealState, P::Union{Missing, Real, PRES}) = begin
    v1 = _v(ξ.𝐺, ξ.𝑃, ξ.𝑇, :MO)
    ξ(P = P, T = _T(ξ.𝐺, P, v1, :MO))
end

isov_T(ξ::IdealState, T::Union{Missing, Real, TEMP}) = begin
    v1 = _v(ξ.𝐺, ξ.𝑃, ξ.𝑇, :MO)
    ξ(P = _P(ξ.𝐺, T, v1, :MO), T = T)
end

isov_u(ξ::IdealState, u::Real, B::Symbol) = begin
    v1 = _v(ξ.𝐺, ξ.𝑃, ξ.𝑇, :MO)
    T2 = isoP_u(ξ, u, B).𝑇
    ξ(P = _P(ξ.𝐺, T2, v1, :MO), T = T2)
end
isov_u(ξ::IdealState, u::ENER) = isov_u(ξ, kSI(u), u isa MOLR ? :MO : :MA)

isov_h(ξ::IdealState, h::Real, B::Symbol) = begin
    v1 = _v(ξ.𝐺, ξ.𝑃, ξ.𝑇, :MO)
    T2 = isoP_h(ξ, h, B).𝑇
    ξ(P = _P(ξ.𝐺, T2, v1, :MO), T = T2)
end
isov_h(ξ::IdealState, h::ENER) = isov_h(ξ, kSI(h), h isa MOLR ? :MO : :MA)

isov_s(ξ::IdealState, s::Real, B::Symbol) = begin
    v1 = _v(ξ.𝐺, ξ.𝑃, ξ.𝑇, :MO)
    T2 = find_zero(
        B == :MO ? T -> kSI(isov_T(ξ, T).sMO) - s : T -> kSI(isov_T(ξ, T).s) - s,
        (ξ.Tmin, ξ.Tmax), Bisection()
    )
    ξ(P = _P(ξ.𝐺, T2, v1, :MO), T = T2)
end
isov_s(ξ::IdealState, s::ENTR) = isov_s(ξ, kSI(s), s isa MOLR ? :MO : :MA)

function isov(
        ξ::IdealState;
        P::Union{Missing, ℙ, PRES{ℙ}} where {ℙ <: Real} = missing,
        T::Union{Missing, ℙ, TEMP{ℙ}} where {ℙ <: Real} = missing,
        u::Union{Missing, Tuple{ℙ, Symbol}, ENER{ℙ}} where {ℙ <: Real} = missing,
        h::Union{Missing, Tuple{ℙ, Symbol}, ENER{ℙ}} where {ℙ <: Real} = missing,
        s::Union{Missing, Tuple{ℙ, Symbol}, ENTR{ℙ}} where {ℙ <: Real} = missing,
    )
    @assert(
        count(x -> !isa(x, Missing), (P, T, u, h, s)) == 1,
        "exactly one end-state function must be specified!"
    )
    return if !ismissing(P)
        isov_P(ξ, P)
    elseif !ismissing(T)
        isov_T(ξ, T)
    elseif !ismissing(u)
        u isa Tuple ? isov_u(ξ, u...) : isov_u(ξ, u)
    elseif !ismissing(h)
        h isa Tuple ? isov_h(ξ, h...) : isov_h(ξ, h)
    elseif !ismissing(s)
        s isa Tuple ? isov_s(ξ, s...) : isov_s(ξ, s)
    end
end

export isov

# Isentropic processes
# --------------------

isos_P(ξ::IdealState, P::Union{Missing, Real, PRES}) = begin
    s1 = _s(ξ.𝐺, ξ.𝑃, ξ.𝑇, :MO)
    ξ(P = P, T = find_zero(T -> kSI(ξ(P = P, T = T).sMO) - s1, (ξ.Tmin, ξ.Tmax), Bisection()))
end

isos_T(ξ::IdealState, T::Union{Missing, Real, TEMP}) = begin
    s1 = _s(ξ.𝐺, ξ.𝑃, ξ.𝑇, :MO)
    ξ(P = find_zero(P -> kSI(ξ(P = P, T = T).sMO) - s1, (1.0e-9, 1.0e+9), Bisection()), T = T)
end

isos_v(ξ::IdealState, v::Real, B::Symbol) = begin
    s1 = _s(ξ.𝐺, ξ.𝑃, ξ.𝑇, :MO)
    T2 = find_zero(T -> kSI(isoT_v(ξ(T = T), v, B).sMO) - s1, (ξ.Tmin, ξ.Tmax), Bisection())
    isoT_v(ξ(T = T2), v, B)
end
isos_v(ξ::IdealState, v::VOLU) = isos_v(ξ, kSI(v), v isa MOLR ? :MO : :MA)

function isos(
        ξ::IdealState;
        P::Union{Missing, ℙ, PRES{ℙ}} where {ℙ <: Real} = missing,
        T::Union{Missing, ℙ, TEMP{ℙ}} where {ℙ <: Real} = missing,
        v::Union{Missing, Tuple{ℙ, Symbol}, VOLU{ℙ}} where {ℙ <: Real} = missing,
    )
    @assert(
        count(x -> !isa(x, Missing), (P, T, v)) == 1,
        "exactly one end-state function must be specified!"
    )
    return if !ismissing(P)
        isos_P(ξ, P)
    elseif !ismissing(T)
        isos_T(ξ, T)
    elseif !ismissing(v)
        v isa Tuple ? isos_v(ξ, v...) : isos_v(ξ, v)
    end
end

export isos
