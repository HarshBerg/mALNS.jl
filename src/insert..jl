function best!(rng::AbstractRNG, s::Solution; mode::Symbol)
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

# TODO: Run any removal operator and then test best insertion.
# TODO: @code_warntype best!(rng, s; mode=:default).

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
    # loop until all nodes are inserted
    for _ ∈ I
        z = f(s)
        for i ∈ I
            n = L[i]
            if isclosed(n) continue end
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
        end
        i,j = argmin(C)
        p = P[i,j]
        n = L[i]
        t = N[p[1]]
        h = N[p[2]]
        v = V[j]
        insertnode!(n, t, h, v, s)
        for i ∈ I
            pivot = L[i]
            if isopen(n) continue end
            t = N[p[1]]
            h = n
            insertnode!(pivot, t, h, v, s)
            z′′ = f(s) 
            Δ = z′′ - z′
            if Δ < C[i,j] C[i,j], P[i,j] = Δ, (t.i, h.i) end
            removenode!(pivot, t, h, v, s)
            t = n
            h = N[p[2]]
            insertnode!(pivot, t, h, v, s)
            z′′ = f(s)
            Δ = z′′ - z′
            if Δ < C[i,j] C[i,j], P[i,j] = Δ, (t.i, h.i) end
            removenode!(pivot, t, h, v, s)
        end
        C[i,:] .= Inf
        P[i,:] .= ((0,0),)    
    end
    return s
end

# TODO: Complete the greedy insertion function.
# NOTE: You don't want to compute Δ multiple times in greedy insertion for routes that have not chaged in the last insertion.

#fprecise(rng::AbstractRNG, s::Solution) = f(rng, s; mode=:precise)
#fpertrub(rng::AbstractRNG, s::Solution) = f(rng, s; mode=:perturb)  