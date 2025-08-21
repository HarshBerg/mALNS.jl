# Arc is line path joining two points
# Tailnode and headnode are the two points

struct Arc
    iᵗ::Int                   # Index of the tail node
    iʰ::Int                   # Index of the head node
    l::Float64                 # Length of the arc                
end

# Route is one loop from depot to depot; includes vehicle index, coordinates of the route, total customers served, length total
mutable struct Route
    iᵛ::Int                   # Vehicle index
    x::Vector{Float64}        # x-Coordinates of the route
    y::Vector{Float64}        # y-Coordinates of the route
    n::Int                    # Number of customers served
    qᶜ::Float64               # Total demand served by the route
    l::Float64                # Length of the route
end

# Vehicle has index, capacity, operating cost, and a vector of routes
mutable struct Vehicle
    iᵛ::Int                   # Vehicle index
    qᵛ::Float64               # Vehicle capacity
    π::Float64                # Vehicle operating cost per unit distance
    R::Vector{Route}          # Vector of routes
end

# node is a point in the graph
abstract type Node end

# Depot is a node with coordinates
struct Depot <: Node
    x::Float64                # x-Coordinate of the depot
    y::Float64                # y-Coordinate of the depot
    V::Vector{Vehicle}        # Vector of vehicles at the depot
end

# Customer is a node with coordinates, index and demand
mutable struct Customer <: Node
    iⁿ::Int                   # Customer index
    x::Float64                # x-Coordinate of the customer
    y::Float64                # y-Coordinate of the customer
    qᶜ::Float64               # Customer demand
    iᵗ::Int                   # Index of the tail node (depot)
    iʰ::Int                   # Index of the head node (depot)
    iʳ::Int                   # Route index serving this customer
    iᵛ::Int                   # Vehicle index serving this customer
end

# solution is a collection of routes and vehicles
mutable struct Solution
    C::Vector{Customer}       # Vector of customers

  



    
