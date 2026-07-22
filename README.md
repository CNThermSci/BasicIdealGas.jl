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

julia> C = SpecificHeat(
    :cubic,             # model ID
    # molar cp(T) model
    T -> 22.26 +5.891e-2*T -3.501e-5*T^2 +7.469e-9*T^3,
    44.01,              # Molecular weight in kg/kmol
    273,                # Minimum T in K
    298,                # Reference T in K
    1800,               # Maximum T in K
    6885,               # Ref internal energy in kJ/kmol
    213.685             # Ref entropy in kJ/kmol/K
    # Omitted molar gas constant (defaults to universal one)
    )
cubic cp₆₄(T)

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

julia> C.𝑓(300)
36.983763

julia> typeof(ans)
Float64
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
cubic cp₃₂(T)

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

*Usage:*

```julia
julia> C.view
SpecificHeat{Float64}(:cubic, var"#2#3"(), 44.01, 273.0, 298.0, 1800.0, 6885.0, 213.685, 8.31447)
          +--------------------------------+               
      1.4 |⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣀⣀⣀⣠| ⠤⠤⠤⠤ [kJ/kg·K]
          |⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⠤⠤⠔⠒⠒⠒⠉⠉⠉⠉⠁⠀⠀⠀⠀⠀|               
   cp (T) |⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⠤⠒⠊⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀|               
          |⠀⠀⠀⠀⠀⣀⠤⠊⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀|               
          |⠀⠀⢀⠔⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀|               
      0.8 |⡠⠊⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀|               
          +--------------------------------+               
          ⠀273⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀T [K]⠀⠀⠀⠀⠀⠀⠀⠀⠀1 800⠀               
julia> cp(C, 1800) # Defaults to mass base
1.327534832992502

julia> C.cp(1800, :MO) # Molar base
58.424808000000006

julia> C.cv(1800, :MO)
50.110338000000006

julia> C.ga(1800) # γ = cp/cv
1.1659232472149759

julia> C.u(1800) # Specific internal energy, mass base
1647.0341409742273

julia> C.h(1800) # Specific enthalpy, mass base
1987.0942636736138

julia> C.s0(1800) # Ideal gas partial entropy, mass base
6.850566852042051

julia> C.Pr(1800) # Relative pressure, Pr = 1 at reference temperature
38596.5956535214

julia> C.vr(1800) # Relative volume
0.04663623745882819
```

### Example 2 – `IdealGas`

`IdealGas` objects adds formula, name, and reference pressure data beyond the specific heat model, thus allowing for ideal gas $P$-$T$-$v$ and entropy calculations. Since these calculations require multiple input parameters, keyword argument versions are provided:

```julia
julia> CO2 = IdealGas("CO2", "Carbon Dioxide", C)
CO2₆₄ gas, cubic cp₆₄(T)

julia> CO2.s(P=100, T=300)
3.9909694845958117

julia> CO2.Pref
1.0

julia> CO2.P(T=300, v=1.2, B=:MA) # v taken in mass base
47.23057259713702

julia> CO2.P(T=300, v=1.2) # If omitted, base defaults to mass
47.23057259713702

julia> CO2.P(T=300, v=1.2, B=:MO) # v taken in molar base
2078.6175

julia> CO2.v(P=47, T=300)
1.2058869599269026
```

### Example 3 – `IdealState`

`IdealState` objects adds state information to `IdealGas`.

Currently only $(P, T)$, positional constructors are implemented:

```julia
julia> st1 = IdealState(CO2, 100, 300)
CO2₆₄ gas, cubic cp₆₄(T) @(100 kPa, 300 K)
```

Since the state is already known, user-facing convenience accessors are implemented for all
the usual thermodynamic state function (thermodynamic properties) through julia properties
"syntactic sugar", such as `st1.v` (mass-based specific volume) and `st1.vMO` (molar-based
specific volume):

```julia
julia> st1.<tab>
ID    M     P     Pr    Pref  R     RMA   T     Tmax
Tmin  Tref  cp    cpMO  cv    cvMO  f     fMA   form
ga    gas   h     hMO   hmod  name  s     s0    s0MO
sMO   sref  u     uMO   uref  v     vMO   vr    γ
ρ     ρMO   𝐺     𝑀     𝑃     𝑅     𝑇     𝑓
```

The user-facing convenience accessors through julia properties return amounts with units, while
"raw" object fields are returned as stored:

```julia
julia> sample_properties = [ st1.v, st1.vMO, st1.u, st1.s ]
4-element Vector{Quantity{Float64}}:
   0.5667668711656442 m^3 kg^-1
  24.94341 m^3 kmol^-1
 157.74275549365427 kJ kg^-1
   3.9909694845958117 kJ kg^-1 K^-1

julia> sample_fields = [ st1.𝑀, st1.𝑅, st1.uref, st1.sref ]
4-element Vector{Float64}:
   44.01
    8.31447
 6885.0
  213.685
```

## Author

Prof. C. Naaktgeboren, PhD. [Lattes](http://lattes.cnpq.br/8621139258082919).

Hermann von Helmholtz Energy Research Group
[DGP](http://dgp.cnpq.br/dgp/espelhogrupo/8462486184187645).

Federal University of Technology, Paraná
[(site)](https://www.utfpr.edu.br/english), Guarapuava Campus.

`NaaktgeborenC <dot!> PhD {at!} gmail [dot!] com`


## License

This project is [licensed](https://github.com/CNThermSci/BasicIdealGas.jl/blob/main/LICENSE)
under the MIT license.


## Citations

Please, refer to the `CITATION.bib` file on how to cite this project.
