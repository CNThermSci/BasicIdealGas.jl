# BasicIdealGas

Basic ideal gas models in engineering thermodynamics.

## Description

`BasicIdealGas.jl` is a package developed in the context of undergraduate mechanical engineering course on internal combustion engine simulation at the equilibrium thermodynamics level, also known as 0-D models. It provides types for basic ideal gas functionality from a hierarchy of `Type`s:

- `SpecificHeat{ℙ <: Base.IEEEFloat}`: A generic ideal gas specific heat model parameterized by
  the precision `ℙ <: Base.IEEEFloat`. `SpecificHeat` methods are the ultimate fallback for
  ideal gas calculations that are solely dependent on temperature, as well as for storing and
  retrieving gas constants, such as the molecular weight and gas constant;

- `IdealGas{ℙ <: Base.IEEEFloat}`: A precision-parametric type for basic ideal gas EoS and entropy
  calculations. `IdealGas{ℙ}` objects include a `SpecificHeat{ℙ}` member data. `IdealGas`
  introduces Equation of State calculations atop of the ones covered by the `SpecificHeat` data
  member, including the ideal gas $P$-$T$-$v$ behavior, as well as entropy, $s:s(P, T)$, ones.

- `IdealState{ℙ <: Base.IEEEFloat}`: A precision-parametric type for an ideal gas at a
  determined `(P, T)` state. Since here the state is known, `IdealState` object instances are
  able to return ideal gas properties (in the Thermodynamic sense) through properties (in the
  julia language sense).

## Common Design Choices

- All data fields are stored as plain `ℙ <: Base.IEEEFloat` types;

- It is _assumed_ that values are in kSI system, i.e., energy in $kJ$, temperatures in $K$,
  pressure in $kPa$, specific volumes in $m³/kmol$ (molar base, `:MO`) or $m³/kg$ (mass base, `:MA`);

- Instantiated `SpecificHeat` object fields store base-dependent data in the molar base;

- `SpecificHeat` functions that return based amounts default to the mass base, `:MA`, due to
  the envisioned engineering applications; however, based amount calculating functions accept a
  `Symbol`ic base argument---either `:MA`, or `:MO`.

- Constructors accept any unambiguous combination of `Real` and `Quantity{<:Real}` arguments;

- User-facing outputs are frequently accessed through fields and properties (in the julia
  langauge sense);

## Examples

### Example 1 – `SpecificHeat`

*Instantiation:*

```julia
julia> using BasicIdealGas

julia> C = SpecificHeat(:cubic, T -> 22.26 + 5.891e-2 * T - 3.501e-5 * T^2 + 7.469e-9 * T^3, 44.01, 273, 298, 1800, 6885, 213.685)
cubic cp₆₄(T) [273.0 1800.0]

julia> dump(C)
SpecificHeat{Float64}
  ID: Symbol cubic
  𝑓: #2 (function of type var"#2#3")
  𝑀: Float64 44.01
  Tmin: Float64 273.0
  Tref: Float64 298.0
  Tmax: Float64 1800.0
  uref: Float64 6885.0
  sref: Float64 213.685
  𝑅: Float64 8.31447

julia> C.f(300)
36.983763
```

It is worth noting that (i) each specific heat model may have it's own gas constant—this is so
due to legacy databases such as NASA Glenn coefficients employing universal gas constants of
slighlty different precision than today's accepted value; (ii) although the `𝑀`, `Tmin`, etc.
values are stored in plain `Base.IEEEFloat`s, as shown, the model function, i.e., the `𝑓` field
has no such return type information, even though it's return type is checked, and does return
consistently typed values, in the above case, a `Float64`.

*Precision conversion:*

```julia
julia> Float32(C)
cubic cp₃₂(T) [273.0 1800.0]

julia> dump(Float32(C))
SpecificHeat{Float32}
  ID: Symbol cubic
  𝑓: Float32 ∘ var"#2#3"() (function of type ComposedFunction{Type{Float32}, var"#2#3"})
    outer: primitive type Float32 <: AbstractFloat
    inner: #2 (function of type var"#2#3")
  𝑀: Float32 44.01f0
  Tmin: Float32 273.0f0
  Tref: Float32 298.0f0
  Tmax: Float32 1800.0f0
  uref: Float32 6885.0f0
  sref: Float32 213.685f0
  𝑅: Float32 8.31447f0

julia> typeof(Float32(C).𝑓)
ComposedFunction{Type{Float32}, var"#2#3"}

julia> Float32(C).𝑓(300)
36.983765f0
```

Julia function composition is used, not only to perform the intended conversions, but also, to
render multiple conversions lossless, i.e., if a `SpecificHeat{Float64}` is converted to a
`Float32` precision, and then back to `Float64`, it preserves the intrinsic precision of the
original model:

```julia
julia> a = [ C.𝑓, Float32(C).𝑓, Float64(Float32(C)).𝑓 ]
3-element Vector{Function}:
 #2 (generic function with 1 method)
 Float32 ∘ var"#2#3"()
 #2 (generic function with 1 method)

julia> a[1] === a[3]
true
```



