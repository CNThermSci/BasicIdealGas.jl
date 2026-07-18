# idealProcs.jl - Ideal Gas Processes

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
        v::Union{Missing, Tuple{ℚ, Symbol}, VOLU{ℚ}} where {ℚ <: Real} = missing,
        s::Union{Missing, Tuple{ℝ, Symbol}, ENTR{ℝ}} where {ℝ <: Real} = missing,
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

# Isobaric processes
# ------------------

# TODO: refactor as with isoT

isoP_T(ξ::IdealState, T::Union{Missing, Real, TEMP}) = ξ(T = T)

isoP_v(ξ::IdealState, v::Real, B::Symbol) = ξ(T = _T(ξ.gas, ξ.𝑃, v, B))
isoP_v(ξ::IdealState, v::VOLU) = isoP_v(ξ, kSI(v), v isa MOLR ? :MO : :MA)

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
        P::Union{Missing, ℙ, PRES{ℙ}} where {ℙ <: Real} = missing,
        v::Union{Missing, Tuple{ℚ, Symbol}, VOLU{ℚ}} where {ℚ <: Real} = missing,
        s::Union{Missing, Tuple{ℝ, Symbol}, ENTR{ℝ}} where {ℝ <: Real} = missing,
    )
    @assert(
        count(x -> !isa(x, Missing), (P, v, s)) == 1,
        "exactly one end-state function must be specified!"
    )
    return if !ismissing(P)
        isoP_T(ξ, T)
    elseif !ismissing(v)
        v isa Tuple ? isoP_v(ξ, v...) : isoP_v(ξ, v)
    elseif !ismissing(s)
        s isa Tuple ? isoP_s(ξ, s...) : isoP_s(ξ, s)
    end
end

export isoP
