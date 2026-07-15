# idealProcs.jl - Ideal Gas Processes

function isoT(
        FR::IdealState{ℙ};
        P::Union{Missing, Real} = missing,
        v::Union{Missing, Real} = missing,
    ) where {ℙ}
    @assert(
        count(x -> isa(x, Real), (P, v)) == 1,
        "exactly one P-T-v state function must be specified!"
    )
    return if !ismissing(P)
        IdealState{ℙ}(FR.gas, P, FR.𝑇)
    elseif !ismissing(v)
        IdealState{ℙ}(FR.gas, _P(FR.gas, FR.𝑇, v, :MA), FR.𝑇)
    end
end

export isoT
