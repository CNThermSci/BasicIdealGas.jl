using InteractiveUtils

# Collects and returns Union members into a DataType[], allowing loop through Union types
function union2vec(theU::Union)
    ret = DataType[]
    while theU isa Union
        push!(ret, theU.a)
        theU = theU.b
    end
    push!(ret, theU)
    return ret
end

@testset "cpModel.test.jl: inner constructor return types                           " begin
    for ℙ in union2vec(Base.IEEEFloat)
        ID, B = :CO2_CB7_test, :MO
        CP = T -> ℙ(22.26 + 5.891e-2 * T - 3.501e-5 * T^2 + 7.469e-9 * T^3)
        Tmin, Tmax, Tref  = 273, 1800, 298
        uref, sref, M, RU = 6885, 213.685, 44.01, ℙ(8.31447)
        pars = ℙ.((M, Tmin, Tmax, Tref, uref, sref))
        @test SpecificHeat(ID, CP, pars..., B, RU) isa SpecificHeat{ℙ}
    end
end

