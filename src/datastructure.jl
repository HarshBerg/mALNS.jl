
# TODO: Add definition
mutable struct Node
    i::Int                   # Customer index
    x::Float64                # x-Coordinate of the customer
    y::Float64                # y-Coordinate of the customer
    q::Float64               # Customer demand
    t::Int                   # Index of the tail node (depot)
    h::Int                   # Index of the head node (depot)
    v::Int                   # Vehicle index serving this customer
    # TODO: Define the new customer
    Node(# some inputs) = new(# all input values)
end

""" 
    Arc(iᵗ::Int, iʰ::Int, l::Float64)
    An 'Arc' is a directed edge from node iᵗ (tail) to node iʰ (head) with length l.
"""
# TODO: Update the definition
struct Arc
    i::Int                   # Index of the tail node
    j::Int                   # Index of the head node
    l::Float64                 # Length of the arc                
end

"""
    Vehicle(iᵛ::Int, qᵛ::Float64, π::Float64, R::Vector{Route})
    A 'Vehicle' has index iᵛ, capacity qᵛ, operating cost π per unit distance, and a vector of routes R.
"""
# TODO: Update the definition
mutable struct Vehicle
    i::Int                   # Vehicle index
    # TODO: Add start and end node Index
    q::Float64               # Vehicle capacity
    n::Int                    # Number of customers served
    l::Float64                # Total demand served by the route
    x::Vector{Float64}        # x-Coordinates of the route
    y::Vector{Float64}        # y-Coordinates of the route
    # TODO: Define the new vehicle
end

"""
    Solution(C::Vector{Customer}, D::Vector{Depot}, A::Vector{Arc}, c::Float64)
    A 'Solution' consists of vectors of customers C, depots D, and arcs A with cost c.
"""
# TODO: Update Definition
mutable struct Solution
    N::Vector{Node}           # Vector of Nodes
    A::Matrix{Arc}            # Matrix of arcs 
    c::Float64                # Cost
end

  



    
