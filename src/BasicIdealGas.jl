module BasicIdealGas

# Imports
using Reexport
@reexport using Unitful
using QuadGK

# Includes
include("common.jl")
include("cpModel.jl")
include("idealEoS.jl")

end
