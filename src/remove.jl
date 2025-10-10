# Random Removal
function randomnode!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    V = s.G.V
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

    return s
end

function randomsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function randomvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

# Related Removal
function relatednode!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    X = fill(-Inf, length(N))
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    iᵖ = sample(rng, 1:length(N), Weights(W))
    for i ∈ eachindex(W)
        X[i] = isone(W[i]) ? relatedness(N[iᵖ], N[i]) : -Inf end

    for _ in 1:k
        i = argmax(X)
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        W[i] = 0
        X[i] = -Inf
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
    # TODO
    return s
end

# Worst Removal
function worstnode!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    W = zeros(Float64, length(N))
    V = s.G.V
    z = f(s)
    for i ∈ eachindex(N)
        if isone(i) continue end
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        z' = f(s)
        W[i] = z' - z
        insertnode!(n, t, h, v, s)
    end
    #OR 
    A = s.G.A
    W = zeros(Float64, length(N))
    for i ∈ eachindex(N)
        if isone(i) continue end
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        W[i] = (A[t.i, n.i].c + A[n.i, h.i].c) - A[t.i, h.i].c
    end
    for _ in 1:k
        i = argmax(W)
        removenode!(n, t, h, v, s)
        z' = f(s)
        W[i] = z - z'
        insertnode!(n, t, h, v, s)
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
    # TODO
    return s
end
