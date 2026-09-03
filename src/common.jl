# Type aliasing
# -------------

# IEEE-754 normalized floating point types of half, single, and double precision
FLOAT = Base.IEEEFloat

# Thermodynamic state function Quantity type alias - dimension set (for arguments)
const PRES = Quantity{ℙ, dimension(u"kPa")} where {ℙ <: Real}
const TEMP = Quantity{ℙ, dimension(u"K")} where {ℙ <: Real}
const MOLW = Quantity{ℙ, dimension(u"kg/kmol")} where {ℙ <: Real}
const VOLU = Union{
    Quantity{ℙ, dimension(u"m^3/kg")},
    Quantity{ℙ, dimension(u"m^3/kmol")},
} where {ℙ <: Real}
const ENER = Union{
    Quantity{ℙ, dimension(u"kJ/kg")},
    Quantity{ℙ, dimension(u"kJ/kmol")},
} where {ℙ <: Real}
const ENTR = Union{
    Quantity{ℙ, dimension(u"kJ/kg/K")},
    Quantity{ℙ, dimension(u"kJ/kmol/K")},
} where {ℙ <: Real}
const DENS = Union{
    Quantity{ℙ, dimension(u"kg/m^3")},
    Quantity{ℙ, dimension(u"kmol/m^3")},
} where {ℙ <: Real}

# Termodynamic base Unions
const MASS = Union{
    Quantity{ℙ, dimension(u"m^3/kg")},
    Quantity{ℙ, dimension(u"kJ/kg")},
    Quantity{ℙ, dimension(u"kJ/kg/K")},
    Quantity{ℙ, dimension(u"kg/m^3")},
} where {ℙ <: Real}
const MOLR = Union{
    Quantity{ℙ, dimension(u"m^3/kmol")},
    Quantity{ℙ, dimension(u"kJ/kmol")},
    Quantity{ℙ, dimension(u"kJ/kmol/K")},
    Quantity{ℙ, dimension(u"kmol/m^3")},
} where {ℙ <: Real}

# Thermodynamic state function Quantity type alias - unit set (for type defs)
const PRESty{ℙ <: FLOAT} = Quantity{ℙ, dimension(u"kPa"), typeof(u"kPa")}
const TEMPty{ℙ <: FLOAT} = Quantity{ℙ, dimension(u"K"), typeof(u"K")}
const VOLUty{ℙ <: FLOAT} = Quantity{ℙ, dimension(u"m^3/kmol"), typeof(u"m^3/kmol")}
const MOLWty{ℙ <: FLOAT} = Quantity{ℙ, dimension(u"kg/kmol"), typeof(u"kg/kmol")}
const ENERty{ℙ <: FLOAT} = Quantity{ℙ, dimension(u"kJ/kmol"), typeof(u"kJ/kmol")}
const ENTRty{ℙ <: FLOAT} = Quantity{ℙ, dimension(u"kJ/kmol/K"), typeof(u"kJ/kmol/K")}
const DENSty{ℙ <: FLOAT} = Quantity{ℙ, dimension(u"kmol/m^3"), typeof(u"kmol/m^3")}

# Thermodynamic unit conversion/stripping
kSI(x::Real) = x
kSI(x::PRES) = uconvert(u"kPa", x).val
kSI(x::TEMP) = uconvert(u"K", x).val
kSI(x::MOLW) = uconvert(u"kg/kmol", x).val

function kSI(x::MASS)
    return if x isa VOLU
        uconvert(u"m^3/kg", x).val
    elseif x isa ENER
        uconvert(u"kJ/kg", x).val
    elseif x isa ENTR
        uconvert(u"kJ/kg/K", x).val
    elseif x isa DENS
        uconvert(u"kg/m^3", x).val
    end
end

function kSI(x::MOLR)
    return if x isa VOLU
        uconvert(u"m^3/kmol", x).val
    elseif x isa ENER
        uconvert(u"kJ/kmol", x).val
    elseif x isa ENTR
        uconvert(u"kJ/kmol/K", x).val
    elseif x isa DENS
        uconvert(u"kmol/m^3", x).val
    end
end

# Constants
# ---------

universal_R = 8.31447

# Function in/out type utilities
# ------------------------------

function type_table(𝑓::Function; hint = 300)
    ret = Tuple{DataType, DataType}[]
    for 𝕌 in (1, u"K")
        for 𝕀 in (Float64, Float32, Float16, Int64)
            𝕚 = 𝕀(hint) * 𝕌
            try
                𝕠 = 𝑓(𝕚)
                push!(ret, typeof.((𝕚, 𝕠)))
            catch
                push!(ret, (typeof(𝕚), Exception))
            end
        end
    end
    return ret
end

is_qty_dim(T, dim) = T <: Quantity && dimension(T) == dim

function profile_user_func(𝑓::Function; hint = 300)
    # Initializations
    iwrap = Dict{Symbol, Union{Function, Nothing}}(:uless => nothing, :ufull => nothing)
    owrap = Dict{Symbol, Union{Function, Nothing}}(:uless => nothing, :ufull => nothing)
    # Helper functions
    Qprec(Q::Type{Quantity{ℙ, 𝔻, 𝕌}}) where {ℙ, 𝔻, 𝕌} = ℙ
    Qprec(Q::Type{Quantity{ℙ, 𝔻}}) where {ℙ, 𝔻} = ℙ
    Qprec(Q::Type{Quantity{ℙ}}) where {ℙ} = ℙ
    # Full i/o type table for function
    table = type_table(𝑓, hint = hint)
    # Exception-filtered
    valid = [ ti for ti in table if ti[2] != Exception ]
    # Separate by Input units
    uless_input = [ ti for ti in valid if ti[1] <: Union{FLOAT, Integer} ]
    ufull_input = [ ti for ti in valid if ti[1] <: Quantity ]
    if length(uless_input) == 0
        # No valid unitless input
        if length(ufull_input) > 0
            # Input must be unitful
            iwrap[:uless] = x -> x * u"K"
            table = type_table(𝑓 ∘ iwrap[:uless], hint = hint)
            valid = [ ti for ti in table if ti[2] != Exception ]
            uless_input = [ ti for ti in valid if ti[1] <: Union{FLOAT, Integer} ]
        else
            # No valid wrappers
            return (iwrap = iwrap, owrap = owrap)
        end
    elseif length(ufull_input) == 0
        iwrap[:ufull] = x -> uconvert(u"K", x).val
    end
    # Homogeneous I/O precision type sets
    same_prec_ul = Set(p[1] for p in [ ti for ti in uless if ti[1] == ti[2] ])
    same_prec_uf = Set(Qprec(p[1]) for p in [ ti for ti in ufull if Qprec(ti[1]) == Qprec(ti[2]) ])

    # Same unitless I/O types for the following:
    same_ul_ty = Set(p[1] for p in [ ti for ti in uless if ti[1] == ti[2] ])
    # Function adds units for the following unitless types:
    addu4 = [ ti for ti in uless if ti[2] <: Quantity ]
    f_adds_units = any(t -> t[2] <: Quantity, uless)
end

# Utilities
# ---------

# Precision Composition Simplification
⊚(p::Type{ℙ}, f::Function) where {ℙ <: FLOAT} = f(300) isa ℙ ? f : p ∘ f

# Chained Precision Composition Simplification
⊚(
    p::Type{ℙ},
    c::ComposedFunction{<:Union{Type{ℚ}, typeof(float)}}
) where {ℙ <: FLOAT, ℚ <: FLOAT} = ⊚(p, c.inner)

# Auxiliary methods
function subscript(x::Int)
    asSub(c::Char) = Char(Int(c) - Int('0') + Int('₀'))
    return map(asSub, "$(x)")
end

pDeco(::Type{Float16}) = subscript(16)
pDeco(::Type{Float32}) = subscript(32)
pDeco(::Type{Float64}) = subscript(64)

# Numerical integrator
# --------------------

function ∫(
        𝑔::Function,
        a::Union{Float32, Float64, Integer, Rational},
        b::Union{Float32, Float64, Integer, Rational},
    )
    ℙ = typeof(promote(a, b)[1])
    return quadgk(𝑔, a, b, rtol = eps(ℙ) * 2 << 6)[1]
end

function ∫(
        𝑔::Function,
        a::Union{Float16, Integer, Rational},
        b::Union{Float16, Integer, Rational},
    )
    a32, b32 = Float32.((a, b))
    n = max(Int(ceil((b32 - a32) / 0.25f0)), 32)
    x32 = range(a32, step = (b32 - a32) / n, length = n + 1) |> collect
    y32 = map(𝑔, x32)
    return Float16(integrate(x32, y32, Trapezoidal()))
end
