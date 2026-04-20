"""
insert!(rng::AbstractRNG, s::Solution; mode::Symbol)

Insert open nodes into the solution `s` using the specified insertion heuristic.

available methods include:
- `bestprecise!` : best insertion with precise evaluation
- `bestperturb!` : best insertion with perturbed evaluation
- `greedyprecise!` : greedy insertion with precise evaluation
- `greedyperturb!` : greedy insertion with perturbed evaluation
- `regret2precise!` : regret-2 insertion with precise evaluation
- `regret2perturb!` : regret-2 insertion with perturbed evaluation
- `regret3precise!` : regret-3 insertion with precise evaluation
- `regret3perturb!` : regret-3 insertion with perturbed evaluation
"""
insert!(rng::AbstractRNG, s::Solution; method::Function)::Solution = method(rng, s)

function best!(rng::AbstractRNG, s::Solution; mode::Symbol)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # intialize
    φ = isequal(mode, :perturb)
    L = [n for n ∈ N if iscustomer(n) && isopen(n)] # set of open nodes
    I = eachindex(L)
    J = eachindex(V)
    W = ones(Int, I)
    C = fill(Inf, (I,J))                            # C[i,j]: best insertion cost of node L[i] in vehicle route V[j]
    P = fill((0,0), (I,J))                          # P[i,j]: best insertion position on node L[i] in vehicle route V[j]
    for _ ∈ I
        z = f(s)
        i = sample(rng, I, Weights(W))
        n = L[i]
        for j ∈ J
            v = V[j]
            t = N[1]
            h = N[v.s]
            for _ ∈ 0:v.n
                insertnode!(n, t, h, v, s)
                z′ = f(s) * (1 + φ * rand(rng, Uniform(-0.2, 0.2)))
                Δ  = z′ - z
                if Δ < C[i,j] C[i,j], P[i,j] = Δ, (t.i, h.i) end
                removenode!(n, t, h, v, s)
                t = h
                h = N[t.h]
            end
        end
        j = argmin(C[i,:])
        p = P[i,j]
        t = N[p[1]]
        h = N[p[2]]
        v = V[j]
        insertnode!(n, t, h, v, s)
        W[i] = 0
        C[i,:] .= Inf
        P[i,:] .= ((0,0),)
    end
    return s
end

bestprecise!(rng::AbstractRNG, s::Solution) = best!(rng, s; mode=:precise)
bestperturb!(rng::AbstractRNG, s::Solution) = best!(rng, s; mode=:perturb)

function greedy!(rng::AbstractRNG, s::Solution; mode::Symbol)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    A = G.A
    # intialize
    φ = isequal(mode, :perturb)
    L = [n for n ∈ N if iscustomer(n) && isopen(n)] # set of open nodes
    I = eachindex(L)
    J = eachindex(V)
    W = ones(Int, I)
    C = fill(Inf, (I,J))                            # C[i,j]: best insertion cost of node L[i] in vehicle route V[j]
    P = fill((0,0), (I,J))                          # P[i,j]: best insertion position on node L[i] in vehicle route V[j]
    ϕ = ones(Int, J)                                #
    # loop until all nodes are inserted
    for _ ∈ I
        z = f(s)
        for i ∈ I
            n = L[i]
            if isclose(n) continue end
            for j ∈ J
                if iszero(ϕ[j]) continue end
                v = V[j]
                t = N[1]
                h = N[v.s]
                for _ ∈ 0:v.n
                    insertnode!(n, t, h, v, s)
                    z′ = f(s) * (1 + φ * rand(rng, Uniform(-0.2, 0.2)))
                    Δ  = z′ - z
                    if Δ < C[i,j] C[i,j], P[i,j] = Δ, (t.i, h.i) end
                    removenode!(n, t, h, v, s)
                    t = h
                    h = N[t.h]
                end
            end
        end
        # insert the best node found
        i, j = Tuple(argmin(C))
        p = P[i,j]
        n = L[i]
        t = N[p[1]]
        h = N[p[2]]
        v = V[j]
        insertnode!(n, t, h, v, s)
        ϕ .= 0
        ϕ[j] = 1
        C[i,:] .= Inf
        P[i,:] .= ((0,0),)
        C[:,j] .= Inf
        P[:,j] .= ((0,0),)
    end
    return s
end

greedyprecise!(rng::AbstractRNG, s::Solution) = greedy!(rng, s; mode=:precise)
greedyperturb!(rng::AbstractRNG, s::Solution) = greedy!(rng, s; mode=:perturb)

function regretk!(rng::AbstractRNG, s::Solution; K::Int, mode::Symbol)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # intialize
    φ = isequal(mode, :perturb)
    L = [n for n ∈ N if iscustomer(n) && isopen(n)] # set of open nodes
    I = eachindex(L)
    J = eachindex(V)
    C = fill(Inf, (I,J))                            # C[i,j]    : best insertion cost of node L[i] in vehicle route V[j]
    P = fill((0,0), (I,J))                          # P[i,j]    : best insertion position on node L[i] in vehicle route V[j]
    Cₖ= fill(Inf, (I,K))                            # Δₖ[i,k]   : k-th least insertion cost of node L[i]
    R = zeros(I)                                    # R[i]      : regret cost of node L[i]
    ϕ = ones(Int, J)                                
    for _ ∈ I
        z = f(s)
        for (i,n) ∈ enumerate(L)
            if isclose(n) continue end
            for (j,v) ∈ enumerate(V)
                if iszero(ϕ[j]) continue end
                t = N[1]
                h = N[v.s]
                for _ ∈ 0:v.n
                    insertnode!(n, t, h, v, s)
                    z′ = f(s) * (1 + φ * rand(rng, Uniform(-0.2, 0.2)))
                    Δ  = z′ - z
                    if Δ < C[i,j] C[i,j], P[i,j] = Δ, (t.i, h.i) end
                    # revise k least insertion cost
                    for (k,Δₖ) ∈ enumerate(Cₖ[i,:])
                        if Δ < Δₖ
                            Cₖ[i,k] = Δ
                            Δ = Δₖ
                        end
                    end
                    removenode!(n, t, h, v, s)
                    t = h
                    h = N[t.h]
                end
            end
            # compute the regret cost
            Δₒ = Cₖ[i,1]
            for (k,Δₖ) ∈ enumerate(Cₖ[i,:]) R[i] += Δₖ - Δₒ end
        end
        # insert the best node found (only consider open nodes)
        i = argmax([isopen(L[ii]) ? R[ii] : -Inf for ii ∈ eachindex(L)])
        j = argmin(C[i,:])
        p = P[i,j]
        n = L[i]
        t = N[p[1]]
        h = N[p[2]]
        v = V[j]
        insertnode!(n, t, h, v, s)
        ϕ .= 0
        ϕ[j] = 1
        C[i,:] .= Inf
        P[i,:] .= ((0,0),)
        C[:,j] .= Inf
        P[:,j] .= ((0,0),)
        Cₖ[i,:] .= Inf
        R .= 0
    end
    return s
end

regret2precise!(rng::AbstractRNG, s::Solution) = regretk!(rng, s; K=2, mode=:precise)
regret2perturb!(rng::AbstractRNG, s::Solution) = regretk!(rng, s; K=2, mode=:perturb)
regret3precise!(rng::AbstractRNG, s::Solution) = regretk!(rng, s; K=3, mode=:precise)
regret3perturb!(rng::AbstractRNG, s::Solution) = regretk!(rng, s; K=3, mode=:perturb)