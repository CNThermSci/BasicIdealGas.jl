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

function user_func_trywrap(𝑓::Function; hint = 300)
    # Initializations
    wrap = Dict{Symbol, Union{Function, Nothing}}(:i => nothing, :o => nothing)
    # Helper functions
    Qprec(Q::Type{Quantity{ℙ, 𝔻, 𝕌}}) where {ℙ, 𝔻, 𝕌} = ℙ
    Qprec(Q::Type{Quantity{ℙ, 𝔻}}) where {ℙ, 𝔻} = ℙ
    Qprec(Q::Type{Quantity{ℙ}}) where {ℙ} = ℙ
    # Full i/o type table for function
    table = type_table(𝑓, hint = hint)
    # Exception-filtered
    valid = [ ti for ti in table if ti[2] != Exception ]
    # INPUT WRAPPING
    valid_uless_input = [ ti for ti in valid if ti[1] <: Union{FLOAT, Integer} ]
    valid_ufull_input = [ ti for ti in valid if ti[1] <: Quantity ]
    if length(valid_uless_input) == 0
        # No valid unitless input
        if length(valid_ufull_input) > 0
            # 𝑓 only accepts Quantity / wrapper's input must be Real (no conversion)
            tmpi(x::Real) = x * u"K"
            tmpi(x::Quantity) = x
            wrap[:i] = tmpi
        else
            # No valid wrappers: return Dict of nothings
            return wrap
        end
    elseif length(valid_ufull_input) == 0
        # No valid unitfull input
        tmpj(x::Real) = x
        tmpj(x::Quantity{𝕋, dimension(u"K")}) where 𝕋 = uconvert(u"K", x).val
        wrap[:i] = tmpj
    end
    if !isnothing(wrap[:i])
        table = type_table(𝑓 ∘ wrap[:i], hint = hint)
        valid = [ ti for ti in table if ti[2] != Exception ]
        valid_uless_input = [ ti for ti in valid if ti[1] <: Union{FLOAT, Integer} ]
        valid_ufull_input = [ ti for ti in valid if ti[1] <: Quantity ]
    end
    # OUTPUT WRAPPING
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
