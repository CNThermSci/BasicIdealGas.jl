# idealProcs.jl - Ideal Gas Processes

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

export isoT

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

