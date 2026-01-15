module mALNS

# TODO: Import packages here
using CSV
using Plots
using Random
using Revise
using StatsBase
using DataFrames
using Distributions

# TODO: Include files here
include("datastructure.jl")
include("functions.jl")
include("operations.jl")
include("initialize.jl")
include("remove.jl")
include("visualize.jl")
include("insert.jl")
include("localsearch.jl")

# TODO: Export function here
export build, initialize, f, visualize

end
