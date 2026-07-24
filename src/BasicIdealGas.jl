module BasicIdealGas

# Imports
using Reexport
@reexport using Unitful

using UnicodePlots
using Printf
using Roots
using QuadGK
using NumericalIntegration

# Includes
include("common.jl")
include("cpModel.jl")
include("idealEoS.jl")
include("propPair.jl")
include("idealState.jl")
include("idealProcs.jl")

end
