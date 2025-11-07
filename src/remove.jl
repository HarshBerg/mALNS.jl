"""
   randomnode!(rng::AbstractRNG, k::Int, s::Solution)

   returns solution 's' after removing exactly k random nodes from solution s.
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

    returns solution 's' after removing exactly k random arcs from solution s.
"""
function randomarc!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    N = s.G.N
    A = s.G.A
    V = s.G.V
    # arc indices 
    I = eachindex(A)
    W = zeros(Int, I)
    # uniform weights for all arcs that connect customer nodes
    for i ∈ I
        a = A[i]
        t = N[a.t]
        h = N[a.h]
        W[i] = isequal(t.h, h.i) ? 1 : 0
    end
    # loop: remove at least k nodes
    c = 0
    while c < k
        # sample an arc based on uniform weights
        i = sample(rng, I, Weights(W))
        a = A[i]
        # remove tail node 
        n = N[a.t]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # remove head node
        n = N[a.h]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # update arc weight to 0 to avoid reselection
        W[i] = 0
    end
    # return solution
    return s
end

function randomsegment!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    N = s.G.N
    A = s.G.A
    V = s.G.V
    W = ones(Int, length(V))
    c = 0 
    while c < k
        # sample a vehicle based on uniform weights
        i = sample(rng, 1:length(V), Weights(W))
        v = V[i] 
        println(v)
        if v.n < 2
            W[i] = 0
            continue
        end
        x = v.n ÷ (rand(rng) * 1.2 + 1.3)   # number of nodes to remove in that vehicle 
        println(x)
        y = v.n - x # possible ways to choose the segment of x nodes
        println(y)
        c2 = 0
        n = N[v.s]
        println(n)
        z = rand(1:y)
        println(z)
        while c2 < x
            c3 = 0
            if c3 < z
                n = N[n.h]
                c3 += 1
            else
                t = N[n.t]
                h = N[n.h]
                removenode!(n, t, h, v, s)
                c3 += 1
            end
            c2 += 1
        end
        W[i] = 0
        c += 1
    end
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
    
    returns solution 's' after removing exactly k nodes based on relatednessto a pivot node.
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

    returns solution 's' after removing exactly k arcs based on relatedness to a pivot arc.
"""
function relatedarc!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    A = s.G.A
    V = s.G.V
    G = s.G
    # arc indices
    I = eachindex(A)
    # choose a pivot arc from solution at random
    p = A[sample(rng, I, Weights([isequal(N[A[i].t].h, N[A[i].h].i) ? 1 : 0 for i ∈ I]))]
    pₜ = N[p.t]
    pₕ = N[p.h]
    # compute relatedness weights based on Manhattan distance to pivot arc
    W = zeros(Float64, I) 
    for i ∈ I
        a = A[i]
        aₜ = N[a.t]
        aₕ = N[a.h]
        d = abs((aₜ.x + aₕ.x) - (pₜ.x + pₕ.x)) + abs((aₜ.y + aₕ.y) - (pₜ.y + pₕ.y))
        W[i] = isequal(aₜ.h, aₕ.i) ? (1 / (d + 1e-3)) : 0.
    end
    c = 0
    while c < k
        # sample an arc based on relatedness weights
        i = sample(rng, 1:length(A), Weights(W))
        a = A[i]
        # remove tail node
        n = N[a.t]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # remove head node
        n = N[a.h]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # update arc weight to 0 to avoid reselection
        W[i] = 0.
    end
    # return solution
    return s
end

function relatedsegment!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    A = s.G.A
    G = s.G
    return s
end
"""
    relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)

    returns solution 's' after removing atleast 'k' nodes from vehicles related to a pivot vehicle.
"""
function relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    G = s.G
    # choose a pivot vehicle at random
    p = V[rand(rng, eachindex(V))]
    # compute relatedness weights based on Manhattan distance to pivot vehicle
    W = zeros(Float64, length(V))
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

    returns solution 's' after removing exactly k nodes that contribute the most to the cost.
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

    returns solution 's' after removing exactly k arcs that contribute the most to the cost.
"""
function worstarc!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    A = s.G.A
    V = s.G.V
    # arc indices
    I = eachindex(A)
    # cost contribution weights for each arc
    W = [isequal(N[A[i].t].h, N[A[i].h].i) ? A[i].c : 0. for i ∈ I]
    c = 0
    while c < k
        # sample an arc based on weights
        i = sample(rng, I, Weights(W))
        a = A[i]
        # remove tail node
        n = N[a.t]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # remove head node
        n = N[a.h]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
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
