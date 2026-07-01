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
        for 𝔽 in (ℙ, float)
            ID, B = :CO2_CB7_test, :MO
            FN = T -> 𝔽(22.26 + 5.891e-2 * T - 3.501e-5 * T^2 + 7.469e-9 * T^3)
            Tmin, Tref, Tmax  = 273, 298, 1800
            uref, sref, M, RU = 6885, 213.685, 44.01, ℙ(8.31447)
            pars = ℙ.((M, Tmin, Tref, Tmax, uref, sref, RU))
            @test SpecificHeat(ID, FN, pars..., B) isa SpecificHeat{ℙ}
        end
    end
end

# Adjacent swaps of a tuple
adjswp(t::Tuple) = [
    ntuple(i -> i == j ? t[j+1] : i == j+1 ? t[j] : t[i], length(t)) for j in 1:length(t)-1
]

@testset "cpModel.test.jl: inner constructor validations                            " begin
    for ℙ in union2vec(Base.IEEEFloat)
        ID, B = :CO2_CB7_test, :MO
        FN = T -> 22.26 + 5.891e-2 * T - 3.501e-5 * T^2 + 7.469e-9 * T^3
        Tneg, Tmin, Tref, Tmax = -1, 273, 298, 1800
        uref, sref, M, RU = 6885, 213.685, 44.01, 8.31447
        pars = ℙ.((M, Tmin, Tref, Tmax, uref, sref, RU))
        @test_throws "Error: Empty ID" SpecificHeat(Symbol(""), FN, pars..., B)
        pars = ℙ.((0, Tmin, Tref, Tmax, uref, sref, RU))
        @test_throws "Error: M <= 0" SpecificHeat(ID, FN, pars..., B)
        pars = ℙ.((-M, Tmin, Tref, Tmax, uref, sref, RU))
        @test_throws "Error: M <= 0" SpecificHeat(ID, FN, pars..., B)
        for temp in adjswp(ℙ.((Tneg, Tmin, Tref, Tmax)))
            pars = (ID, FN, ℙ.((M, temp[2:end]..., uref, sref, RU))..., B)
            @test_throws "Error: Temperature values" SpecificHeat(pars...)
        end
        pars = ℙ.((M, Tmin, Tref, Tmax, uref, sref, 0))
        @test_throws "Error: RU <= 0" SpecificHeat(ID, FN, pars..., B)
        pars = ℙ.((M, Tmin, Tref, Tmax, uref, sref, -RU))
        @test_throws "Error: RU <= 0" SpecificHeat(ID, FN, pars..., B)
        pars = (ID, FN, ℙ.((M, Tmin, Tref, Tmax, uref, sref, RU))...)
        for b in (:ma, :mo, :other, Symbol(""))
            @test_throws "Error: B should be either :MA or :MO" SpecificHeat(pars..., b)
        end
    end
end

@testset "cpModel.test.jl: outer constructor return types                           " begin
    ID, B = :CO2_CB7_test, :MO
    FN = T -> 22.26 + 5.891e-2 * T - 3.501e-5 * T^2 + 7.469e-9 * T^3
    Tmin, Tref, Tmax = 273, 298, 1800
    uref, sref, M, RU = 6885, 213685//1000, BigFloat("44.01"), π
    # Set type conversion / 1 indirection
    for ℙ in union2vec(Base.IEEEFloat)
        pars = (ID, FN, M, Tmin, Tref, Tmax, uref, sref, RU, B)
        @test SpecificHeat{ℙ}(pars...) isa SpecificHeat{ℙ}
    end
    # Set type with unit conversion and stripping / 2 indirections
    for ℙ in union2vec(Base.IEEEFloat)
        pars = (ID, FN, M*u"kg/kmol", Tmin, Tref, Tmax, uref*u"kJ/kmol", sref*u"kJ/kmol/K")
        @test SpecificHeat{ℙ}(pars..., RU) isa SpecificHeat{ℙ}
    end
    uref, sref, M, RU = 6885, 213685//1000, 4401//100, 8.31447
    # Promotion type conversion / 2 indirections
    for ℙ in union2vec(Base.IEEEFloat)
        pars = (ID, FN, M, Tmin, Tref, Tmax, uref, sref, ℙ(RU), B)
        @test SpecificHeat(pars...) isa SpecificHeat{ℙ}
    end
    # Promotion type with unit conversion and stripping / 3 indirections
    for ℙ in union2vec(Base.IEEEFloat)
        pars = (ID, FN, M*u"kg/kmol", Tmin, Tref, Tmax, uref*u"kJ/kmol", sref*u"kJ/kmol/K")
        @test SpecificHeat(pars..., ℙ(RU)) isa SpecificHeat{ℙ}
    end
end

