# BasicIdealGas

Basic ideal gas models in engineering thermodynamics.

## Description

`BasicIdealGas.jl` is a package developed in the context of mechanical engineering education. It
provides types for basic ideal gas functionality from a hierarchy of `Type`s:

- `SpecificHeat{ℙ <: Base.IEEEFloat}`: A generic ideal gas specific heat model parameterized by
  the precision `ℙ <: Base.IEEEFloat`. `SpecificHeat` methods are the ultimate fallback for
  ideal gas calculations that are solely dependent on temperature, as well as for storing and
  retrieving gas constants, such as the molecular weight and gas constant;

- `IdealGas{ℙ <: Base.IEEEFloat}`: A precision-parametric type for basic ideal gas EoS and
  entropy calculations. `IdealGas{ℙ}` objects include a `SpecificHeat{ℙ}` member data.
  `IdealGas` introduces Equation of State calculations atop of the ones covered by the
  `SpecificHeat` data member, including the ideal gas $P-T-v$ behavior, as well as entropy,
  $s:s(P, T)$, ones.

- `IdealState{ℙ <: Base.IEEEFloat}`: A precision-parametric type for an ideal gas at a
  determined `(P, T)` state. Since here the state is known, `IdealState` object instances are
  able to return ideal gas properties (in the Thermodynamic sense) through properties (in the
  julia language sense).

## Common Design Choices

- All dimensional data fields are stored with units as `Quantity{ℙ} where ℙ <: Base.IEEEFloat`
  types;

- Specific heat functions are dimensionless and molar-base normalized by the universal gas
  constant;

- Stored value units are in the kSI system, and specific quantities in the molar base, `:MO`,
  rather than in the mass base, `:MA`, i.e., energy in $kJ$, temperatures in $K$, pressure in
  $kPa$, specific internal energies in $kJ/kmol$, and specific entropies in $kJ/kmol/K$;

- User-facing outputs are accessed through fields and properties (in the julia langauge sense).
  When the amount is based, a `Symbol`ic base argument—whether `:MO`, or `:MA`, respectively for
  molar or mass base—can be optionally specified, with the mass base being the default one.

- Constructors accept any unambiguous combination of `Real` and `Quantity{<:Real}` arguments.
  Specifying unitless values of reference specific internal energy and entropy is ambiguous.
  Unitless molecular mass is assumed to be in $kg/kmol$, and temperatures in $K$.

## Examples

### Example 1 – `SpecificHeat`

*Instantiation:*

```julia
julia> using BasicIdealGas

julia> C = SpecificHeat(
    :cubic,             # model ID
    # molar cp(T)/Ru model
    T -> (22.26 +5.891e-2*T -3.501e-5*T^2 +7.469e-9*T^3) / Ru.val,
    44.01,              # Molecular weight, kg/kmol
    273,                # Minimum T, K
    298,                # Reference T, K
    1800,               # Maximum T, K
    6885u"kJ/kmol",     # Ref internal energy
    213.685u"kJ/kmol/K" # Ref entropy
    # Omitted molar gas constant (defaults to universal one)
    )
cubic cp₆₄(T)

julia> typeof(C)
SpecificHeat{Float64}
```

The constructor wraps the passed function into one that only accepts temperature arguments, and
returns values as the precision parameter of `SpecificHeat{ℙ} where ℙ <: Base.IEEEFloat`:

```julia
julia> C.f┆R(300u"K")
4.448124274352035

julia> typeof(ans)
Float64

julia> C.f┆R
#2 (generic function with 1 method)

julia> C.R
8.31446261815324 kJ K^-1 kmol^-1
```

It is worth noting that (i) each specific heat model may have it's own gas constant—this is so
due to legacy databases such as NASA Glenn coefficients employing the universal gas constant of
CODATA 1986; (ii) the model function, i.e., the wrapped `f┆R` field is composed with the
`SpecificHeat` object's precision parameter if the return type of the function passed upon
construction is different, which isn't the case for the `SpecificHeat{Float64}` object—since the
provided function already returns a `Float64` value—but is the case for the converted
`SpecificHeat{Float32}` object below:

*Precision conversion:*

```julia
julia> Float32(C)
cubic cp₃₂(T)

julia> typeof(Float32(C).f┆R)
ComposedFunction{Type{Float32}, BasicIdealGas.var"#2#3"{Float32, var"#11#12"}}

julia> Float32(C).f┆R(300u"K")
4.4481244f0
```

Julia function composition is used, not only to perform the intended conversions, but also, to
render multiple _function_ conversions lossless, i.e., if a `SpecificHeat{Float64}` is converted to a
`Float32` precision, and then back to `Float64`, it preserves the intrinsic precision of the
original model function (but not the other parameters!):

```julia
julia> a = [ C.f┆R, Float32(C).f┆R, Float64(Float32(C)).f┆R ]
3-element Vector{Function}:
 #2 (generic function with 1 method)
 Float32 ∘ BasicIdealGas.var"#2#3"{Float32, var"#11#12"}(var"#11#12"())
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
CO2 gas, cubic cp₆₄(T)

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

### Example 3 - `PropPair`

`PropPair` designates a thermodynamic $(P, T)$ property pair, and serves to determine the state
of an ideal gas inside a simple compressible system.

```julia
julia> p = PropPair(100, 300)
@₆₄(100 kPa, 300 K)

julia> dump(p)
PropPair{Float64}
  𝑃: Float64 100.0
  𝑇: Float64 300.0

julia> Float32(p)
@₃₂(100 kPa, 300 K)
```

`PropPair` has some rough edges to be trimmed on upcomming releases.

### Example 4 - `Interact`

`Interact` represents simple compressible system interactions of heat and work.

```julia
julia> i = Interact(-10.0, -20.0)
Interact{Float64}(-10.0, -20.0)

julia> dump((i, Float32(i)))
Tuple{Interact{Float64}, Interact{Float32}}
  1: Interact{Float64}
    𝑞: Float64 -10.0
    𝑤: Float64 -20.0
  2: Interact{Float32}
    𝑞: Float32 -10.0f0
    𝑤: Float32 -20.0f0

julia> i.q
-10.0 kJ kg^-1

julia> i.w
-20.0 kJ kg^-1
```

`Interact` is very incipient, and has significant rough edges to be trimmed on upcomming releases.

### Example 5 – `IdealState`

`IdealState` objects adds state (through a property pair, `PropPair`) information to `IdealGas`.

```julia
julia> st1 = IdealState(CO2, PropPair(100, 300))
CO2 gas, cubic cp₆₄(T) @₆₄(100 kPa, 300 K)
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
ρ     ρMO   𝐺     𝑀     𝑅     𝑓     𝑝
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

### Example 6 - Processes

Some incipient ideal gas processes _functions_ are available:

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

julia> CO2 = IdealGas("CO2", "Carbon Dioxide", C)
CO2 gas, cubic cp₆₄(T)

julia> st1 = IdealState(CO2, PropPair(100, 300))
CO2 gas, cubic cp₆₄(T) @₆₄(100 kPa, 300 K)

julia> st2 = isoP(st1, T = 400)             # Isobaric process up to T = 400 (K)
CO2 gas, cubic cp₆₄(T) @₆₄(100 kPa, 400 K)

julia> st3 = isoT(st2, P = 150)             # Isothermal process up to P = 150 (kPa)
CO2 gas, cubic cp₆₄(T) @₆₄(150 kPa, 400 K)

julia> st3.u
227.21082809079132 kJ kg^-1

julia> st4 = isov(st3, u = (250, :MA))      # Isochoric process up to u = 250 (kJ/kg)
CO2 gas, cubic cp₆₄(T) @₆₄(161.43 kPa, 430.48 K)

julia> st4.s
4.222878505947657 kJ kg^-1 K^-1

julia> st5 = isos(st4, P = st1.𝑃)           # Isentropic process up to st1 pressure
CO2 gas, cubic cp₆₄(T) @₆₄(100 kPa, 390.68 K)

julia> st5.s
4.222878505947658 kJ kg^-1 K^-1
```

Process interactions are not yet being calculated and returned from process functions.

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
