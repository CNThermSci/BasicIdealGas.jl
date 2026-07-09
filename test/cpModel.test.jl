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
            ID, B = :cubic, :MO
            𝑓 = T -> 𝔽(22.26 + 5.891e-2 * T - 3.501e-5 * T^2 + 7.469e-9 * T^3)
            Tmin, Tref, Tmax = 273, 298, 1800
            uref, sref, 𝑀, 𝑅 = 6885, 213.685, 44.01, ℙ(8.31447)
            pars = ℙ.((𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅))
            @test SpecificHeat(ID, 𝑓, pars..., B) isa SpecificHeat{ℙ}
        end
    end
end

# Adjacent swaps of a tuple
adjswp(t::Tuple) = [
    ntuple(i -> i == j ? t[j + 1] : i == j + 1 ? t[j] : t[i], length(t)) for j in 1:(length(t) - 1)
]

@testset "cpModel.test.jl: inner constructor validations                            " begin
    for ℙ in union2vec(Base.IEEEFloat)
        ID, B = :cubic, :MO
        𝑓 = T -> 22.26 + 5.891e-2 * T - 3.501e-5 * T^2 + 7.469e-9 * T^3
        Tneg, Tmin, Tref, Tmax = -1, 273, 298, 1800
        uref, sref, 𝑀, 𝑅 = 6885, 213.685, 44.01, 8.31447
        pars = ℙ.((0, Tmin, Tref, Tmax, uref, sref, 𝑅))
        @test_throws "Error: Empty model ID" SpecificHeat(Symbol(""), 𝑓, pars..., B)
        @test_throws "Error: M <= 0" SpecificHeat(ID, 𝑓, pars..., B)
        pars = ℙ.((-𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅))
        @test_throws "Error: M <= 0" SpecificHeat(ID, 𝑓, pars..., B)
        for temp in adjswp(ℙ.((Tneg, Tmin, Tref, Tmax)))
            pars = (ID, 𝑓, ℙ.((𝑀, temp[2:end]..., uref, sref, 𝑅))..., B)
            @test_throws "Error: Temperature values" SpecificHeat(pars...)
        end
        pars = ℙ.((𝑀, Tmin, Tref, Tmax, uref, sref, 0))
        @test_throws "Error: 𝑅 <= 0" SpecificHeat(ID, 𝑓, pars..., B)
        pars = ℙ.((𝑀, Tmin, Tref, Tmax, uref, sref, -𝑅))
        @test_throws "Error: 𝑅 <= 0" SpecificHeat(ID, 𝑓, pars..., B)
        pars = (ID, 𝑓, ℙ.((𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅))...)
        for b in (:ma, :mo, :other, Symbol(""))
            @test_throws "Error: B should be either :MA or :MO" SpecificHeat(pars..., b)
        end
    end
end

@testset "cpModel.test.jl: outer constructor return types                           " begin
    ID, B = :cubic, :MO
    𝑓 = T -> 22.26 + 5.891e-2 * T - 3.501e-5 * T^2 + 7.469e-9 * T^3
    Tmin, Tref, Tmax = 273, 298, 1800
    uref, sref, 𝑀, 𝑅 = 6885, 213685 // 1000, BigFloat("44.01"), π
    # Set type conversion / 1 indirection
    for ℙ in union2vec(Base.IEEEFloat)
        pars = (ID, 𝑓, 𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅, B)
        @test SpecificHeat{ℙ}(pars...) isa SpecificHeat{ℙ}
    end
    # Set type with unit conversion and stripping / 2 indirections
    for ℙ in union2vec(Base.IEEEFloat)
        pars = (ID, 𝑓, 𝑀 * u"kg/kmol", Tmin, Tref, Tmax, uref * u"kJ/kmol", sref * u"kJ/kmol/K")
        @test SpecificHeat{ℙ}(pars..., 𝑅) isa SpecificHeat{ℙ}
    end
    uref, sref, 𝑀, 𝑅 = 6885, 213685 // 1000, 4401 // 100, 8.31447
    # Promotion type conversion / 2 indirections
    for ℙ in union2vec(Base.IEEEFloat)
        pars = (ID, 𝑓, ℙ(𝑀), Tmin, Tref, Tmax, uref, sref, 𝑅, B)
        @test SpecificHeat(pars...) isa SpecificHeat{ℙ}
    end
    # Promotion type with unit conversion and stripping / 3 indirections
    for ℙ in union2vec(Base.IEEEFloat)
        pars = (ID, 𝑓, ℙ(𝑀 * u"kg/kmol"), Tmin, Tref, Tmax, uref * u"kJ/kmol", sref * u"kJ/kmol/K")
        @test SpecificHeat(pars..., 𝑅) isa SpecificHeat{ℙ}
    end
end

@testset "cpModel.test.jl: constructor's optional arguments                         " begin
    ID = :cubic
    𝑓 = T -> 22.26 + 5.891e-2 * T - 3.501e-5 * T^2 + 7.469e-9 * T^3
    Tmin, Tref, Tmax = 273, 298, 1800
    uref, sref, 𝑀, 𝑅 = 6885, 213685 // 1000, 4401 // 100, BasicIdealGas.universal_R
    # Promotion type conversion / 2 indirections
    for ℙ in union2vec(Base.IEEEFloat)
        pars = (ID, 𝑓, ℙ(𝑀), Tmin, Tref, Tmax, uref, sref, 𝑅)
        MOLR = SpecificHeat(pars..., :MO)
        MASS = SpecificHeat(pars..., :MA)
        AUTO = SpecificHeat(pars...)
        @test MOLR != MASS
        @test MOLR == AUTO
        @test MASS != AUTO
        auto = SpecificHeat(pars[1:(end - 1)]...)
        @test auto == AUTO
    end
    # Set type conversion / 1 indirection
    for ℙ in union2vec(Base.IEEEFloat)
        pars = (ID, 𝑓, 𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅)
        MOLR = SpecificHeat{ℙ}(pars..., :MO)
        MASS = SpecificHeat{ℙ}(pars..., :MA)
        AUTO = SpecificHeat{ℙ}(pars...)
        @test MOLR != MASS
        @test MOLR == AUTO
        @test MASS != AUTO
        auto = SpecificHeat{ℙ}(pars[1:(end - 1)]...)
        @test auto == AUTO
    end
    # Internal constructor / no indirection
    for ℙ in union2vec(Base.IEEEFloat)
        pars = (ID, 𝑓, ℙ.((𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅))...)
        MOLR = SpecificHeat(pars..., :MO)
        MASS = SpecificHeat(pars..., :MA)
        AUTO = SpecificHeat(pars...)
        @test MOLR != MASS
        @test MOLR == AUTO
        @test MASS != AUTO
        auto = SpecificHeat(pars[1:(end - 1)]...)
        @test auto == AUTO
    end
end

@testset "cpModel.test.jl: type conversions                                         " begin
    ID, B = :cubic, :MO
    𝑓 = T -> 22.26 + 5.891e-2 * T - 3.501e-5 * T^2 + 7.469e-9 * T^3
    Tmin, Tref, Tmax = 273.0, 298.0, 1800.0
    uref, sref, 𝑀, 𝑅 = 6885.0, 213.685, 44.01, 8.31447
    pars = (ID, 𝑓, 𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅, B)
    SH = Dict(
        Float16 => SpecificHeat{Float16}(pars...),
        Float32 => SpecificHeat{Float32}(pars...),
        Float64 => SpecificHeat{Float64}(pars...),
    )
    # Conversions
    for orig in union2vec(Base.IEEEFloat)
        for dest in union2vec(Base.IEEEFloat)
            # Lossless 𝑓 conversion
            @test SH[orig].𝑓 === orig(SH[dest]).𝑓
            # Type conversion
            @test typeof(SH[orig]) === typeof(orig(SH[dest]))
        end
    end
    # Lossless narrowings
    @test SH[Float64] === Float64(SH[Float64])
    @test SH[Float32] === Float32(SH[Float64])
    @test SH[Float32] === Float32(SH[Float32])
    @test SH[Float16] === Float16(SH[Float64])
    @test SH[Float16] === Float16(SH[Float32])
    @test SH[Float16] === Float16(SH[Float16])
end

@testset "cpModel.test.jl: type promotions                                          " begin
    ID, B = :cubic, :MO
    𝑓 = T -> 22.26 + 5.891e-2 * T - 3.501e-5 * T^2 + 7.469e-9 * T^3
    Tmin, Tref, Tmax = 273.0, 298.0, 1800.0
    uref, sref, 𝑀, 𝑅 = 6885.0, 213.685, 44.01, 8.31447
    pars = (ID, 𝑓, 𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅, B)
    SH = Dict(
        Float16 => SpecificHeat{Float16}(pars...),
        Float32 => SpecificHeat{Float32}(pars...),
        Float64 => SpecificHeat{Float64}(pars...),
    )
    # Promotions
    for orig in union2vec(Base.IEEEFloat)
        for dest in union2vec(Base.IEEEFloat)
            @test all(
                [
                    i isa SpecificHeat{promote_type(orig, dest)}
                        for i in promote(SH[orig], SH[dest])
                ]
            )
        end
    end
end

@testset "cpModel.test.jl: user-facing functions: bound temperature intervals       " begin
    # Bounds checks
    for FUNC in (
            :bounds, :cp┆R, :cv┆R, :ga,
            :cp, :cv,
            :∫cp┆R, :∫cv┆R, :u┆R, :h┆R,
            :u, :h,
            :∫cp┆RT, :s0┆R,
            :s0, :Pr, :vr,
        )
        @test_throws "T out of bounds" eval(
            quote
                bounds = BasicIdealGas.bounds
                ID, 𝑓 = :const, T -> 22.26
                Tmin, Tref, Tmax = 273, 298, 1800
                uref, sref, 𝑀 = 6885, 213.685, 44.01
                C = SpecificHeat(ID, 𝑓, 𝑀, Tmin, Tref, Tmax, uref, sref)
                $FUNC(C, prevfloat(C.Tmin))
            end
        )
        @test_throws "T out of bounds" eval(
            quote
                bounds = BasicIdealGas.bounds
                ID, 𝑓 = :const, T -> 22.26
                Tmin, Tref, Tmax = 273, 298, 1800
                uref, sref, 𝑀 = 6885, 213.685, 44.01
                C = SpecificHeat(ID, 𝑓, 𝑀, Tmin, Tref, Tmax, uref, sref)
                $FUNC(C, nextfloat(C.Tmax))
            end
        )
    end
end
