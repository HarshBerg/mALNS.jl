
"""
    Node(i::Int, x::Float64, y::Float64, q::Float64, t::Int, h::Int, v::Int)
    A 'Node' has index i, coordinates (x, y), demand q, tail node index t, head node index h,
    and vehicle index v serving this customer.
"""
mutable struct Node
    i::Int                   # node index
    x::Float64               # x-Coordinate of the node
    y::Float64               # y-Coordinate of the node
    q::Float64               # node demand
    t::Int                   # Index of the tail node (depot)
    h::Int                   # Index of the head node (depot)
    v::Int                   # Vehicle index serving this customer
    # TODO: Define the new customer
    #Node (some inputs) = new(# all input values)
end

""" 
    Arc(i::Int, j::Int, l::Float64)
    An 'Arc' is a directed edge from node i (tail) to node j (head) with length l.
"""
struct Arc
    i::Int                   # Index of the tail node
    j::Int                   # Index of the head node
    l::Float64                 # Length of the arc                
end

"""
    Vehicle(i::Int, s::Int, e::Int, q::Float64, n::Int, l::Float64, x::Vector{Float64}, y::Vector{Float64})
    A 'Vehicle' has index i, starts at node s, ends at node e, has capacity q, serves n customers,
    total demand l, and follows a route defined by coordinates x and y.
"""
mutable struct Vehicle
    i::Int                    # Vehicle index
    s::Int                    # Start node index
    e::Int                    # End node index
    q::Float64                # Vehicle capacity
    n::Int                    # Number of customers served
    l::Float64                # Total demand served by the route
    
end

"""
    Solution(N::Vector{Node}, A::Matrix{Arc}, c::Float64)
    A 'Solution' consists of a vector of nodes N, a matrix of arcs A, and a total cost c.
"""
mutable struct Solution
    N::Vector{Node}           # Vector of Nodes
    A::Matrix{Arc}            # Matrix of arcs 
    c::Float64                # Cost
end

  



    
