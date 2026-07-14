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
  determined `(P, T)` state.

## Common Design Choices

- All data fields are stored as plain `ℙ <: Base.IEEEFloat` types;

- Instantiated `SpecificHeat` objects store caloric data in the molar base, `:MO`;

- `SpecificHeat` functions that return based amounts default to the mass base, `:MA`, due to
  the envisioned engineering applications; however, based amount calculating functions accept a
  `Symbol`ic base argument---either `:MA`, or `:MO`.

- 
