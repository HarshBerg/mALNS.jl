module mALNS

# TODO: Import packages here
using CSV
using Plots
using Revise
using DataFrames

# TODO: Include files here
include("datastructure.jl")
include("functions.jl")
include("operations.jl")
include("initialize.jl")
include("visualize.jl")

# TODO: Export function here
export build, initialize, visualize

end
