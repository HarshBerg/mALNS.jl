""" 
    Arc(iᵗ::Int, iʰ::Int, l::Float64)
    An 'Arc' is a directed edge from node iᵗ (tail) to node iʰ (head) with length l.
"""
struct Arc
    iᵗ::Int                   # Index of the tail node
    iʰ::Int                   # Index of the head node
    l::Float64                 # Length of the arc                
end


"""
    Route(iᵛ::Int, x::Vector{Float64}, y::Vector{Float64}, n::Int, qᶜ::Float64, l::Float64)
    A 'Route' represents a loop from depot to depot served by vehicle iᵛ, with coordinates (x,y), number of customers n, total demand qᶜ, and length l.
"""
mutable struct Route
    iᵛ::Int                   # Vehicle index
    x::Vector{Float64}        # x-Coordinates of the route
    y::Vector{Float64}        # y-Coordinates of the route
    n::Int                    # Number of customers served
    qᶜ::Float64               # Total demand served by the route
    l::Float64                # Length of the route
end


"""
    Vehicle(iᵛ::Int, qᵛ::Float64, π::Float64, R::Vector{Route})
    A 'Vehicle' has index iᵛ, capacity qᵛ, operating cost π per unit distance, and a vector of routes R.
"""
# Vehicle has index, capacity, operating cost, and a vector of routes
mutable struct Vehicle
    iᵛ::Int                   # Vehicle index
    qᵛ::Float64               # Vehicle capacity
    π::Float64                # Vehicle operating cost per unit distance
    R::Vector{Route}          # Vector of routes
end


"""
    Node
    An abstract type 'Node' representing a node in the network.
"""
abstract type Node end


"""
    Depot(x::Float64, y::Float64, V::Vector{Vehicle})
    A 'Depot' is a node with coordinates (x,y) and a vector of vehicles V.
"""
struct Depot <: Node
    x::Float64                # x-Coordinate of the depot
    y::Float64                # y-Coordinate of the depot
    V::Vector{Vehicle}        # Vector of vehicles at the depot
end


"""
    Customer(iⁿ::Int, x::Float64, y::Float64, qᶜ::Float64, iᵗ::Int, iʰ::Int, iʳ::Int, iᵛ::Int)
    A 'Customer' has index iⁿ, coordinates (x,y), demand qᶜ, indices of tail and head depots (iᵗ, iʰ), route index iʳ, and vehicle index iᵛ.
"""
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


"""
    Solution(C::Vector{Customer}, D::Vector{Depot}, A::Vector{Arc})
    A 'Solution' consists of vectors of customers C, depots D, and arcs A.
"""
mutable struct Solution
    C::Vector{Customer}       # Vector of customers
    D::Vector{Depot}          # Vector of depots
    A::Vector{Arc}            # Vector of arcs 
end

  



    
