
"""
    Node(i::Int, x::Float64, y::Float64, q::Float64, t::Int, h::Int, v::Int)
    A 'Node' with index `i`, coordinates `(x, y)`, demand `q`, tail node index `t`, head node index `h`,
    and vehicle index `v` serving this node.
"""
mutable struct Node
    i::Int                   # Index
    x::Int                   # Abcissa (x-coordinate) 
    y::Int                   # Ordinate (y-coordinate)
    q::Int                   # Demand
    t::Int                   # Index of the tail node
    h::Int                   # Index of the head node
    v::Int                   # Vehicle index serving this node
    Node(i, x, y, q) = new(i, x, y, q, 0, 0, 0)
end

"""
    Arc(t::Int, h::Int, c::Float64)
    An 'Arc' is a directed edge from node `t` (tail) to node `h` (head) with cost `c`.
"""
struct Arc
    t::Int                    # Index of the tail node
    h::Int                    # Index of the head node
    c::Float64                # Cost of the arc    
end

"""
    Vehicle(i::Int, q::Int, s::Int, e::Int, n::Int, l::Int, x::Float64, y::Float64, c::Float64)
    A 'Vehicle' with index `i`, capacity `q`, start node index `s`, end node index `e`, number of customers served `n`,
    total vehicle load `l`, centroid coordinates `(x, y)`, and cost `c
"""
mutable struct Vehicle
    i::Int                    # Vehicle index
    q::Int                    # Vehicle capacity
    s::Int                    # Start node index
    e::Int                    # End node index
    n::Int                    # Number of customers served
    l::Int                    # Total vehicle load
    x::Float64                # Centroid abscissa (x-coordinate) 
    y::Float64                # Centroid ordinate (y-coordinate)
    Vehicle(i, q) = new(i, q, 1, 1, 0, 0, 0., 0.)
end

"""
    Solution(N::Vector{Node}, A::Matrix{Arc}, V::Vector{Vehicle})
A 'Graph' consists of a vector of nodes `N`, a matrix of arcs `A`, and a vector of vehicles `V`.
"""
struct Graph
    N::Vector{Node}           # nodes
    A::Matrix{Arc}            # arcs
    V::Vector{Vehicle}        # vehicles
end

"""
    Solution(G::Graph, c::Float64, p::Float64)
    A 'Solution' pertains to a graph `G` with total cost `c` and penalty `p`.
"""
mutable struct Solution
    G::Graph                  # graph
    c::Float64                # Cost
    p::Float64                # Penalty
    Solution(G) = new(G, 0., 0.)
end
