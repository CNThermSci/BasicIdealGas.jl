# idealProcs.jl - Ideal Gas Processes

isoT_P(
    ξ::IdealState,
    P::Union{Missing, Real, Quantity{<:Real, dimension(u"kPa")}},
) = ξ(P = P)

isoT_v(ξ::IdealState, v::Real, B::Symbol) = ξ(P = _P(ξ.gas, ξ.𝑇, v, B))
isoT_v(ξ::IdealState, v::Quantity{<:Real, dimension(u"m^3/kg")}) =
    isoT_v(ξ, uconvert(u"m^3/kg", v).val, :MA)
isoT_v(ξ::IdealState, v::Quantity{<:Real, dimension(u"m^3/kmol")}) =
    isoT_v(ξ, uconvert(u"m^3/kmol", v).val, :MO)

isoT_s(ξ::IdealState{ℙ}, s::Quantity{<:Real, dimension(u"kJ/kg/K")}) where {ℙ} =
    ξ(P = ξ.P * exp((ξ.s - ℙ(s)) / ξ.R))
isoT_s(ξ::IdealState{ℙ}, s::Quantity{<:Real, dimension(u"kJ/kmol/K")}) where {ℙ} =
    ξ(P = ξ.P * exp((ξ.sMO - ℙ(s)) / ξ.RMO))
isoT_s(ξ::IdealState, s::Real, B::Symbol) = 
    B == :MA ? isoT_s(ξ, s * u"kJ/kg/K") : isoT_s(ξ, s * u"kJ/kmol/K")

function isoT(
        FR::IdealState;
        P::Union{
            Missing,
            Quantity{<:Real, dimension(u"kPa")},
            Real,
        } = missing,
        v::Union{
            Missing,
            Quantity{<:Real, dimension(u"m^3/kg")},
            Quantity{<:Real, dimension(u"m^3/kmol")},
            Tuple{<:Real, Symbol},
        } = missing,
        s::Union{
            Missing,
            Quantity{<:Real, dimension(u"kJ/kg/K")},
            Quantity{<:Real, dimension(u"kJ/kmol/K")},
            Tuple{<:Real, Symbol},
        } = missing,
    )
    @assert(
        count(x -> !isa(x, Missing), (P, v, s)) == 1,
        "exactly one end-state function must be specified!"
    )
    return if !ismissing(P)
        isoT_P(ξ, P)
    elseif !ismissing(v)
        isoT_v(ξ, v)
    elseif !ismissing(s)
        isoT_s(ξ, s)
    end
end

export isoT

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
