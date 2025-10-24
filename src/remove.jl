"""
   randomnode!(rng::AbstractRNG, k::Int, s::Solution)

   returns solution 's' after removing k random nodes from solution s.
"""
function randomnode!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    # Weighting: depot nodes get weight 0, customer nodes get weight 1
    # uniform weights for all customer nodes
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    # loop: remove exactly k nodes
    for _ in 1:k
        i = sample(rng, 1:length(N), Weights(W))
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        W[i] = 0
    end
    # return solution
    return s
end
"""
    randomarc!(rng::AbstractRNG, k::Int, s::Solution)

    returns solution 's' after removing k random arcs from solution s.
"""
function randomarc!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    A = s.G.A
    V = s.G.V
    # arc indices 
    C = CartesianIndices(G.A)
    # uniform weights for all arcs that connect customer nodes
    W = [isequal(N[a.t].h, N[a.h]) ? 1 : 0 for a ∈ A]
    for _ in 1:k
        # sample an arc based on uniform weights
        i = sample(rng, 1:length(A), Weights(W))
        a = A[C[i]]
        # remove tail node 
        n = N[a.t]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
        end
        # remove head node
        n = N[a.h]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
        end
        # update arc weight to 0 to avoid reselection
        W[i] = 0
    end
    # return solution
    return s
end

function randomsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end
"""
    randomvehicle!(rng::AbstractRNG, k::Int, s::Solution)

    returns solution 's' after removing atleast 'k' nodes from random vehicles.
"""
function randomvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    # uniform weights for all vehicles
    W = ones(Int, length(V))
    # loop: remove atleast k nodes
    c = 0
    while c < k
        # sample a vehicle based on uniform weights
        i = sample(rng, 1:length(V), Weights(W))
        v = V[i]
        # remove all nodes from the selected vehicle
        for _ ∈ 1:v.n
            n = N[v.s]
            t = N[n.t]
            h = N[n.h]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # update vehicle weight to 0 to avoid reselection
        W[i] = 0
    end
    # return solution
    return s
end

"""
    relatednode!(rng::AbstractRNG, k::Int, s::Solution)
    
    returns solution 's' after removing k nodes based on relatednessto a pivot node.
"""
function relatednode!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    W = zeros(Float64, length(N))
    # choose a pivot node at random
    p = N[rand(rng, eachindex(N))]
    # compute relatedness weights based on Manhattan distance to pivot node
    for i ∈ eachindex(N)
        n = N[i]
        if isdepot(n) continue end
        d = abs(n.x - p.x) + abs(n.y - p.y)
        W[i] = 1 / (d + 1e-3)
    end
    # remove k nodes based on relatedness weights
    for _ in 1:k
        i = sample(rng, 1:length(N), Weights(W))
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        W[i] = 0.
    end
    # return solution
    return s
end
"""
    relatedarc!(rng::AbstractRNG, k::Int, s::Solution)

    returns solution 's' after removing k arcs based on relatedness to a pivot arc.
"""
function relatedarc!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    A = s.G.A
    V = s.G.V
    # arc indices
    C = CartesianIndices(G.A)
    W = zeros(Float64, length(A))
    # choose a pivot arc from solution at random
    p = sample(rng, A, Weights([isequal(N[a.t].h, N[a.h]) ? 1 : 0 for a ∈ A]))
    pₜ = N[p.t]
    pₕ = N[p.h]
    # compute relatedness weights based on Manhattan distance to pivot arc
    for i ∈ eachindex(A)
        a = A[C[i]]
        aₜ = N[a.t]
        aₕ = N[a.h]
        d = abs((aₜ.x + aₕ.x) - (pₜ.x + pₕ.x)) + abs((aₜ.y + aₕ.y) - (pₜ.y + pₕ.y))
        W[i] = isequal(aₜ.h, aₕ) ? (1 / (d + 1e-3)) : 0.
    end
    for _ in 1:k
        # sample an arc based on relatedness weights
        i = sample(rng, 1:length(A), Weights(W))
        a = A[C[i]]
        # remove tail node
        n = N[a.t]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
        end
        # remove head node
        n = N[a.h]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
        end
        # update arc weight to 0 to avoid reselection
        W[i] = 0.
    end
    # return solution
    return s
end

function relatedsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end
"""
    relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)

    returns solution 's' after removing atleast 'k' nodes from vehicles related to a pivot vehicle.
"""
function relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    # choose a pivot vehicle at random
    p = V[rand(rng, eachindex(V))]
    # compute relatedness weights based on Manhattan distance to pivot vehicle
    W = zeros(Float64, I)
    for i ∈ eachindex(V)
        v = G.V[i]
        d = abs(v.x - p.x) + abs(v.y - p.y)
        W[i] = 1 / (d + 1e-3)
    end
    # remove atleast k nodes based on relatedness weights
    c = 0
    while c < k
        i = sample(rng, 1:length(V), Weights(W))
        v = V[i]
        for _ ∈ 1:v.n
            n = N[v.s]
            t = N[n.t]
            h = N[n.h]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # update vehicle weight to 0 to avoid reselection
        W[i] = 0
    end
    # return solution
    return s
end

"""
    worstnode!(rng::AbstractRNG, k::Int, s::Solution)

    returns solution 's' after removing k nodes that contribute the most to the cost.
"""
function worstnode!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    A = s.G.A
    W = zeros(Float64, length(N))
    # compute cost contribution weights for each node
    for i ∈ eachindex(N)
        if isone(i) continue end
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        W[i] = (A[t.i, n.i].c + A[n.i, h.i].c) - A[t.i, h.i].c
    end
    # remove k nodes based on weights
    for _ in 1:k
        i = sample(rng, 1:length(N), Weights(W))
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        # set weight to 0 to avoid reselection
        W[i] = 0.
    end  
    # return solution
    return s
end
"""
    worstarc!(rng::AbstractRNG, k::Int, s::Solution)

    returns solution 's' after removing k arcs that contribute the most to the cost.
"""
function worstarc!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    A = s.G.A
    V = s.G.V
    # arc indices
    C = CartesianIndices(G.A)
    # cost contribution weights for each arc
    W = [isequal(N[a.t].h, N[a.h]) ? a.c : 0. for a ∈ A]
    for _ in 1:k
        # sample an arc based on weights
        i = sample(rng, 1:length(A), Weights(W))
        a = A[C[i]]
        # remove tail node
        n = N[a.t]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
        end
        # remove head node
        n = N[a.h]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
        end
        # update arc weight to 0 to avoid reselection
        W[i] = 0
    end
    # return solution
    return s
end

function worstsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end
"""
    worstvehicle!(rng::AbstractRNG, k::Int, s::Solution)

    returns solution 's' after removing atleast 'k' nodes from the less utilized vehicles.
"""
function worstvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    # utilization based weights for each vehicle
    W = [1 - v.l / v.q for v ∈ V]
    # remove atleast k nodes from less utilized vehicles
    c = 0
    while c < k
        i = sample(rng, 1:length(V), Weights(W))
        v = V[i]
        for _ ∈ 1:v.n
            n = N[v.s]
            t = N[n.t]
            h = N[n.h]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # update vehicle weight to 0 to avoid reselection
        W[i] = 0
    end
    # return solution
    return s
end
