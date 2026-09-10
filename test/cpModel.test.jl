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

cp_R = Dict(
    :const => T -> 5//2,
    :cubic => T -> begin
        t = T / 1000
        r = 8314 // 1000
        ((2226//100) + (5891//100)*t -(3501//100)*t^2 +(7469//1000)*t^3) / r
    end,
)

function cu_cp_R(ℙ::Type{T} where {T <: Base.IEEEFloat})
    return T -> begin
        t = ℙ(T / 1000)
        r = ℙ(8314 // 1000)
        c = map(ℙ, [(2226//100) , (5891//100), -(3501//100), (7469//1000)])
        sum([ c[i+1]*t^i for i in 0:3 ]) / r
    end
end

@testset "cpModel.test.jl: inner constructor return types                           " begin
    for ℙ in union2vec(Base.IEEEFloat)
        for 𝔽 in (union2vec(Base.IEEEFloat)..., float)
            ID = :cubic
            f┆R = 𝔽 ∘ cp_R[:cubic]
            Tmin, Tref, Tmax = 273u"K", 298u"K", 1800u"K"
            uref, sref, 𝑀, 𝑅 = 6885u"kJ/kmol", 213.685u"kJ/kmol/K", 44.01u"kg/kmol", Ru
            pars = ℙ.((𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅))
            @test SpecificHeat(ID, f┆R, pars...) isa SpecificHeat{ℙ}
            f┆R = 𝔽 ∘ cu_cp_R(ℙ)
            # Constructor function smart composition simplifies f┆R away into a Function:
            @test !(SpecificHeat(ID, f┆R, pars...).f┆R isa ComposedFunction)
            # Inner f┆R not a float-returning function exception
            g┆R = 𝔽 ∘ cp_R[:const]
            @test SpecificHeat(ID, g┆R, pars...).f┆R.f┆R isa ComposedFunction
        end
    end
end

# Adjacent swaps of a tuple
adjswp(t::Tuple) = [
    ntuple(i -> i == j ? t[j + 1] : i == j + 1 ? t[j] : t[i], length(t)) for j in 1:(length(t) - 1)
]

@testset "cpModel.test.jl: inner constructor validations                            " begin
    for ℙ in union2vec(Base.IEEEFloat)
        ID = :cubic
        f┆R = cp_R[:cubic]
        Tneg, Tmin, Tref, Tmax = -300u"K", 273u"K", 298u"K", 1800u"K"
        uref, sref, 𝑀, 𝑅 = 6885u"kJ/kmol", 213.685u"kJ/kmol/K", 44.01u"kg/kmol", Ru
        pars = ℙ.((0, Tmin, Tref, Tmax, uref, sref, 𝑅))
        @test_throws "Error: Empty model ID" SpecificHeat(Symbol(""), f┆R, pars...)
        @test_throws "Error: M <= 0" SpecificHeat(ID, f┆R, pars...)
        pars = ℙ.((-𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅))
        @test_throws "Error: M <= 0" SpecificHeat(ID, f┆R, pars...)
        for temp in adjswp(ℙ.((Tneg, Tmin, Tref, Tmax)))
            pars = (ID, f┆R, ℙ.((𝑀, temp[2:end]..., uref, sref, 𝑅))...)
            @test_throws "Error: Temperature values" SpecificHeat(pars...)
        end
        pars = ℙ.((𝑀, Tmin, Tref, Tmax, uref, sref, 0Ru))
        @test_throws "Error: 𝑅 <= 0" SpecificHeat(ID, f┆R, pars...)
        pars = ℙ.((𝑀, Tmin, Tref, Tmax, uref, sref, -𝑅))
        @test_throws "Error: 𝑅 <= 0" SpecificHeat(ID, f┆R, pars...)
    end
end

@testset "cpModel.test.jl: outer constructor return types                           " begin
    ID = :cubic
    f┆R = cp_R[:cubic]
    Tmin, Tref, Tmax = 273, 298, 1800
    uref, sref = 6885u"kJ/kmol", (213685 // 1000)u"kJ/kmol/K"
    𝑀, 𝑅 = BigFloat("44.01"), π * u"kJ/kmol/K"
    # Set type conversion
    for ℙ in union2vec(Base.IEEEFloat)
        pars = (ID, f┆R, 𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅)
        @test SpecificHeat{ℙ}(pars...) isa SpecificHeat{ℙ}
        pars = (ID, f┆R, 𝑀 * u"kg/kmol", Tmin * u"K", Tref, Tmax, uref, sref)
        @test SpecificHeat{ℙ}(pars..., 𝑅) isa SpecificHeat{ℙ}
    end
    𝑀, 𝑅 = 4401//100, 8.31447 * u"kJ/kmol/K"
    # Promotion type conversion
    for ℙ in union2vec(Base.IEEEFloat)
        pars = (ID, f┆R, ℙ(𝑀), Tmin, Tref, Tmax, uref, sref, 𝑅)
        @test SpecificHeat(pars...) isa SpecificHeat{ℙ}
        pars = (ID, f┆R, ℙ(𝑀 * u"kg/kmol"), Tmin * u"K", Tref, Tmax, uref, sref, 𝑅)
        @test SpecificHeat(pars...) isa SpecificHeat{ℙ}
    end
end

@testset "cpModel.test.jl: type conversions                                         " begin
    ID = :cubic
    f┆R = cp_R[:cubic]
    Tmin, Tref, Tmax = 273.0, 298.0, 1800.0
    uref, sref, 𝑀, 𝑅 = 6885u"kJ/kmol", 213.685u"kJ/kmol/K", 44.01u"kg/kmol", Ru
    pars = (ID, f┆R, 𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅)
    SH = Dict(
        Float16 => SpecificHeat{Float16}(pars...),
        Float32 => SpecificHeat{Float32}(pars...),
        Float64 => SpecificHeat{Float64}(pars...),
    )
    # Conversions
    for orig in union2vec(Base.IEEEFloat)
        for dest in union2vec(Base.IEEEFloat)
            # Lossless f┆R conversion
            @test SH[orig].f┆R === orig(SH[dest]).f┆R
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
    ID = :cubic
    f┆R = cp_R[:cubic]
    Tmin, Tref, Tmax = 273.0, 298.0, 1800.0
    uref, sref, 𝑀, 𝑅 = 6885u"kJ/kmol", 213.685u"kJ/kmol/K", 44.01u"kg/kmol", Ru
    pars = (ID, f┆R, 𝑀, Tmin, Tref, Tmax, uref, sref, 𝑅)
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
    bounds = BasicIdealGas.𝗯
    ID, f┆R = :const, T -> 22.26
    Tmin, Tref, Tmax = 273, 298, 1800
    uref, sref, 𝑀 = 6885u"kJ/kmol", 213.685u"kJ/kmol/K", 44.01u"kg/kmol"
    C = SpecificHeat(ID, f┆R, 𝑀, Tmin, Tref, Tmax, uref, sref)
    @test_throws AssertionError bounds(C, prevfloat(C.Tmin))
    @test_throws AssertionError bounds(C, nextfloat(C.Tmax))
end

@testset "cpModel.test.jl: user-facing functions: thermodynamic consistencies       " begin
    cp┆R = BasicIdealGas.cp┆R
    cv┆R = BasicIdealGas.cv┆R
    ga = BasicIdealGas.ga
    R = BasicIdealGas.R
    cp = BasicIdealGas.cp
    cv = BasicIdealGas.cv
    ∫cp┆R = BasicIdealGas.∫cp┆R
    ∫cv┆R = BasicIdealGas.∫cv┆R
    u┆R = BasicIdealGas.u┆R
    h┆R = BasicIdealGas.h┆R
    u = BasicIdealGas.u
    h = BasicIdealGas.h
    ∫cp┆RT = BasicIdealGas.∫cp┆RT
    s0┆R = BasicIdealGas.s0┆R
    s0 = BasicIdealGas.s0
    Pr = BasicIdealGas.Pr
    vr = BasicIdealGas.vr
    # Float16 are tested but may overflow depending on model function form and argument type
    for ℙ in [Float32, Float64]
        Tmin, Tref, Tmax = 273u"K", 298u"K", 1800u"K"
        uref, sref, 𝑀 = 6885u"kJ/kmol", 213.685u"kJ/kmol/K", 44.01u"kg/kmol"
        𝑅 = Ru
        C = SpecificHeat{ℙ}(:cubic, cp_R[:cubic], 𝑀, Tmin, Tref, Tmax, uref, sref)
        G = SpecificHeat{ℙ}(:const, cp_R[:const], 𝑀, Tmin, Tref, Tmax, uref, sref)
        for T in (Tmin, (Tmin + Tmax) / 2, Tmax)
            @test C.f┆R(T) isa ℙ
            @test cp┆R(C, T) ≈ C.f┆R(T)
            @test cv┆R(C, T) ≈ C.f┆R(T) - one(ℙ)
            @test ga(C, T) ≈ cp┆R(C, T) / cv┆R(C, T) ≈ cp(C, T) / cv(C, T)
            @test R(C, :MO) == C.𝑅
            @test R(C, :MA) ≈ C.𝑅 / C.𝑀
            for B in (:MA, :MO)
                @test cp(C, T, B) ≈ cp┆R(C, T) * R(C, B)
                @test cv(C, T, B) ≈ cv┆R(C, T) * R(C, B)
                @test cp(C, T, B) ≈ cv(C, T) + R(C, B)
                @test ga(C, T) ≈ cp(C, T, B) / cv(C, T, B)
            end
            @test ∫cp┆R(G, T) ≈ (5 // 2) * (ℙ(T) - C.Tref)
            @test ∫cv┆R(G, T) ≈ (3 // 2) * (ℙ(T) - C.Tref)
            for H in (C, G)
                @test ∫cv┆R(H, T) ≈ ∫cp┆R(H, T) - (ℙ(T) - C.Tref)
                @test u┆R(H, T) ≈ ∫cv┆R(H, T) + C.uref / C.𝑅
                @test h┆R(H, T) ≈ u┆R(H, T) + ℙ(T)
            end
            for B in (:MA, :MO)
                @test u(C, T) ≈ u┆R(C, T) * R(C)
                @test h(C, T) ≈ h┆R(C, T) * R(C)
                @test h(C, T) ≈ u(C, T) + R(C) * ℙ(T)
            end
            @test ∫cp┆RT(G, T) ≈ (5 // 2) * log(ℙ(T) / C.Tref)
            @test s0┆R(G, T) ≈ ∫cp┆RT(G, T) + C.sref / C.𝑅
            for B in (:MA, :MO)
                @test s0(C, T) ≈ s0┆R(C, T) * R(C)
            end
            @test Pr(C, C.Tref) ≈ one(ℙ)
            @test Pr(C, T) ≈ exp(∫cp┆RT(C, T))
            @test vr(C, T) * Pr(C, T) ≈ ℙ(T)
        end
    end
end
