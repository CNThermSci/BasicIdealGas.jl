# idealProcs.jl - Ideal Gas Processes

function isoT(
        FR::IdealState{ℙ};
        P::Union{
            Missing,
            Quantity{<:Real, dimension(u"kPa")},
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
        count(x -> !isa(x, Missing), (P, v, s)) == 1,
        "exactly one end-state function must be specified!"
    )
    return if !ismissing(P)
        FR(P = uconvert(u"kPa", ℙ(P)).val)
    elseif !ismissing(v)
        if dimension(v) == dimension(u"m^3/kg")
            FR(P = _P(FR.gas, FR.𝑇, uconvert(u"m^3/kg", ℙ(v)).val, :MA))
        else
            FR(P = _P(FR.gas, FR.𝑇, uconvert(u"m^3/kmol", ℙ(v)).val, :MO))
        end
    elseif !ismissing(s)
        if dimension(s) == dimension(u"kJ/kg/K")
            FR(P = FR.P * exp((FR.s - ℙ(s)) / FR.R))
        else
            FR(P = FR.P * exp((FR.sMO - ℙ(s)) / FR.RMO))
        end
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
            # FR(T = )
        else
        end
    end
end

export isoP
