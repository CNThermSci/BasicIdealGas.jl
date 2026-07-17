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
        P::Union{Missing, ℙ, PRES{ℙ}} where ℙ<:Real = missing,
        v::Union{Missing, Tuple{ℚ, Symbol}, VOLU{ℚ}} where ℚ<:Real = missing,
        s::Union{Missing, Tuple{ℝ, Symbol}, ENTR{ℝ}} where ℝ<:Real = missing,
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

function isoP(
        FR::IdealState{ℙ};
        T::Union{
            Missing,
            Quantity{<:Real, dimension(u"K")},
        } = missing,
        v::Union{
            Missing,
            Quantity{<:Real, dimension(u"m^3/kg")},
            Quantity{<:Real, dimension(u"m^3/kmol")},
        } = missing,
        s::Union{
            Missing,
            Quantity{<:Real, dimension(u"kJ/kg/K")},
            Quantity{<:Real, dimension(u"kJ/kmol/K")},
        } = missing,
    ) where {ℙ}
    @assert(
        count(x -> !isa(x, Missing), (T, v, s)) == 1,
        "exactly one end-state function must be specified!"
    )
    return if !ismissing(T)
        FR(T = uconvert(u"K", ℙ(T)).val)
    elseif !ismissing(v)
        if dimension(v) == dimension(u"m^3/kg")
            FR(T = _T(FR.gas, FR.𝑃, uconvert(u"m^3/kg", ℙ(v)).val, :MA))
        else
            FR(T = _T(FR.gas, FR.𝑃, uconvert(u"m^3/kmol", ℙ(v)).val, :MO))
        end
    elseif !ismissing(s)
        if dimension(s) == dimension(u"kJ/kg/K")
            FR(T = find_zero(T -> FR(T = T).s - ℙ(s), (FR.Tmin, FR.Tmax), Bisection()))
        else
            FR(T = find_zero(T -> FR(T = T).sMO - ℙ(s), (FR.Tmin, FR.Tmax), Bisection()))
        end
    end
end

export isoP
