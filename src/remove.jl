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
    N = s.G.N
    A = s.G.A
    V = s.G.V
    C = CartesianIndices(G.A)
    W = [isequal(N[a.i].h, N[a.j]) ? 1 : 0 for a ∈ A]
    for _ in 1:k
        i = sample(rng, 1:length(A), Weights(W))
        a = A[C[i]]
        n = N[a.i]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
        end
        n = N[A.j]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
        end
        W[i] = 0
    end
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
    p = N[rand(rng, eachindex(N))]
    for i ∈ eachindex(N)
        n = G.N[i]
        if isdepot(n) continue end
        d = abs(n.x - p.x) + abs(n.y - p.y)
        W[i] = 1 / (d + 1e-3)
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
    N = s.G.N
    A = s.G.A
    V = s.G.V
    C = CartesianIndices(G.A)
    W = zeros(Float64, length(A))
    p = A[rand(rng, C)]
    pᵢ = N[p.i]
    pⱼ = N[p.j]
    for i ∈ eachindex(A)
        a = A[C[i]]
        aᵢ = N[a.i]
        aⱼ = N[a.j]
        d = abs((aᵢ.x + aⱼ.x) - (pᵢ.x + pⱼ.x)) + abs((aᵢ.y + aⱼ.y) - (pᵢ.y + pⱼ.y))
        W[i] = isequal(aᵢ.h, aⱼ) ? (1 / (d + 1e-3)) : 0.
    end
    for _ in 1:k
        i = sample(rng, 1:length(A), Weights(W))
        a = A[C[i]]
        n = N[a.i]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
        end
        n = N[A.j]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
        end
        W[i] = 0.
    end
    return s
end

function relatedsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    V = s.G.V
    p = V[rand(rng, eachindex(V))]
    W = zeros(Float64, I)
    for i ∈ eachindex(V)
        v = G.V[i]
        d = abs(v.x - p.x) + abs(v.y - p.y)
        W[i] = 1 / (d + 1e-3)
    end
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
    N = s.G.N
    A = s.G.A
    V = s.G.V
    C = CartesianIndices(G.A)
    W = [isequal(N[a.i].h, N[a.j]) ? a.c : 0. for a ∈ A]
    for _ in 1:k
        i = sample(rng, 1:length(A), Weights(W))
        a = A[C[i]]
        n = N[a.i]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
        end
        n = N[A.j]
        if iscustomer(n) && isclose(n) 
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
        end
        W[i] = 0
    end
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
