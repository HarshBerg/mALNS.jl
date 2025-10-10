# Random Removal
function randomnode!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    for _ in 1:k
        i = sample(rng, 1:length(N), Weights(W))
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        W[i] = 0
    end
    return s
end

function randomarc!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function randomsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function randomvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    W = ones(Int, length(V))
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
        W[i] = 0
    end
    return s
end

# Related Removal
function relatednode!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    W = zeros(Float64, length(N))
    p = sample(rng, 2:length(N))
    for i ∈ eachindex(N)
        if isone(i) continue end
        W[i] = relatedness(N[p], N[i])
    end
    for _ in 1:k
        i = sample(rng, 1:length(N), Weights(W))
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        W[i] = 0.
    end
    return s
end

function relatedarc!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function relatedsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    p = sample(rng, 1:length(V))
    W = [relatedness(V[p], V[i]) for i ∈ eachindex(V)]
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
        W[i] = 0
    end
    return s
end

# Worst Removal
function worstnode!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    A = s.G.A
    W = zeros(Float64, length(N))
    for i ∈ eachindex(N)
        if isone(i) continue end
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        W[i] = (A[t.i, n.i].c + A[n.i, h.i].c) - A[t.i, h.i].c
    end
    for _ in 1:k
        i = sample(rng, 1:length(N), Weights(W))
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        W[i] = 0.
    end  
    return s
end

function worstarc!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function worstsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function worstvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    W = [1 - v.l / v.q for v ∈ V]
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
        W[i] = 0
    end
    return s
end
