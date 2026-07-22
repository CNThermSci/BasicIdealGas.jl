using QuadGK
using NumericalIntegration

import Base: ComposedFunction
import BasicIdealGas: FLOAT, universal_R, subscript, pDeco, ⊚, ∫

@testset "common.test.jl: type aliases and constants                                " begin
    @test FLOAT === Base.IEEEFloat
    @test universal_R == 8.31447
    @test universal_R isa Float64
end

@testset "common.test.jl: subscript utilities (subscript, pDeco)                    " begin
    @test subscript(0) == "₀"
    @test subscript(16) == "₁₆"
    @test subscript(123) == "₁₂₃"
    @test pDeco(Float16) == "₁₆"
    @test pDeco(Float32) == "₃₂"
    @test pDeco(Float64) == "₆₄"
end

@testset "common.test.jl: precision composition simplification ⊚                    " begin
    # Named functions to allow === comparison
    f64(x) = Float64(x)
    f32(x) = Float32(x)
    f16(x) = Float16(x)

    # Basic method: returns f unchanged when f(1) already has precision ℙ
    @test ⊚(Float64, f64) === f64
    @test ⊚(Float32, f32) === f32
    @test ⊚(Float16, f16) === f16

    # Basic method: otherwise returns a ComposedFunction p ∘ f
    let comp = ⊚(Float64, f32)
        @test comp isa ComposedFunction
        @test comp(1) === Float64(1.0)
    end

    # Chained method: strips outer float-type / float wrapper
    let c1 = Float32 ∘ f64
        @test ⊚(Float64, c1) === f64
    end
    let c2 = Float64 ∘ (Float32 ∘ sin)
        @test ⊚(Float64, c2) === sin
    end
    let c3 = float ∘ sin
        res = ⊚(Float32, c3)
        @test res isa ComposedFunction
        @test res(1) === Float32(sin(1))
    end
end

@testset "common.test.jl: numerical integrator ∫                                    " begin
    # First method (Float64 / Float32 / mixed)
    @test ∫(identity, 0.0, 1.0) ≈ 0.5
    @test typeof(∫(identity, 0.0, 1.0)) === Float64

    @test ∫(x -> x^2, 0.0, 2.0) ≈ 8 / 3
    @test typeof(∫(x -> x^2, 0.0, 2.0)) === Float64

    @test ∫(identity, 0.0f0, 1.0f0) ≈ 0.5f0
    @test typeof(∫(identity, 0.0f0, 1.0f0)) === Float32

    # Mixed Integer / Float64 promotes to Float64
    @test ∫(identity, 0, 1.0) ≈ 0.5
    @test typeof(∫(identity, 0, 1.0)) === Float64

    # Rational with Float64 (first method only)
    @test ∫(identity, 0 // 1, 1.0) ≈ 0.5
    @test typeof(∫(identity, 0 // 1, 1.0)) === Float64

    # Second method (Float16 bounds, Integer / Rational mixed with Float16)
    @test ∫(identity, Float16(0), Float16(1)) ≈ Float16(0.5)
    @test typeof(∫(identity, Float16(0), Float16(1))) === Float16

    @test ∫(identity, 0, Float16(1)) ≈ Float16(0.5)
    @test typeof(∫(identity, 0, Float16(1))) === Float16

    @test ∫(identity, 0 // 1, Float16(1)) ≈ Float16(0.5)
    @test typeof(∫(identity, 0 // 1, Float16(1))) === Float16
end
