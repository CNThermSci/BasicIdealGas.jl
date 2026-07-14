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


