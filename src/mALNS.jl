module mALNS

using CSV
using Plots
using Random
using Revise
using StatsBase
using DataFrames
using ProgressMeter
using Distributions

include("datastructure.jl")
include("functions.jl")
include("initialize.jl")
include("operations.jl")
include("remove.jl")
include("insert.jl")
include("localsearch.jl")
include("parameters.jl")
include("ALNS.jl")
include("visualize.jl")

export  Node, Arc, Vehicle, Graph, Solution,
        build, initialize,
        randomnode!, randomarc!, randomsegment!, randomvehicle!,
        relatednode!, relatedarc!, relatedsegment!, relatedvehicle!,
        worstnode!, worstarc!, worstsegment!, worstvehicle!,
        bestprecise!, bestperturb!,
        greedyprecise!, greedyperturb!,
        regret2precise!, regret2perturb!,
        regret3precise!, regret3perturb!,
        intramove!, intermove!,
        intraswap!, interswap!,
        intraopt!, interopt!,
        vectorize, sol, isfeasible, f, h,
        ALNSparameters, ALNS, 
        visualize, pltcnv

end
