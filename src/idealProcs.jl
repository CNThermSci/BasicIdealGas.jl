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
        count(x -> !isa(x, Missing), (P, v)) == 1,
        "exactly one P-T-v state function must be specified!"
    )
    return if !ismissing(P)
        IdealState{ℙ}(FR.gas, uconvert(u"kPa", P).val, FR.𝑇)
    elseif !ismissing(v)
        if dimension(v) == dimension(u"m^3/kg")
            IdealState{ℙ}(
                FR.gas,
                _P(FR.gas, FR.𝑇, uconvert(u"m^3/kg", v).val, :MA),
                FR.𝑇
            )
        else
            IdealState{ℙ}(
                FR.gas,
                _P(FR.gas, FR.𝑇, uconvert(u"m^3/kmol", v).val, :MO),
                FR.𝑇
            )
        end
    end
end

export isoT
