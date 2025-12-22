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
    W = ones(Int, I)
    C = fill(Inf, (I,J))                            # C[i,j]    : best insertion cost of node L[i] in vehicle route V[j]
    P = fill((0,0), (I,J))                          # P[i,j]    : best insertion position on node L[i] in vehicle route V[j]
    Cₖ= fill(Inf, (I,K))                            # Cₖ[i,k]   : k-th least insertion cost of node L[i]
    R = zeros(I)                                    # R[i]      : regret cost of node L[i]
    ϕ = ones(Int, J)                                #
    for _ ∈ I
        z = f(s)
        i = sample(rng, I, Weights(W))
        n = L[i]
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
                # revise k least insertion cost
                for (k,Δₖ) ∈ enumerate(Cₖ[i,:])
                    if Δ < Δₖ
                        # do this
                    else
                        # do this
                    end
                end
                removenode!(n, t, h, v, s)
                t = h
                h = N[t.h]
            end
            # compute the regret cost
            Δₒ = Cₖ[i,1]
            for (k,Δₖ) ∈ enumerate(Cₖ[i,:]) R[i] += Δₖ - Δₒ end
        end
        # insert the best node found
        i = argmax(R)
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
        R .= 0
    end
    return s
end

#fprecise(rng::AbstractRNG, s::Solution) = f(rng, s; mode=:precise)
#fpertrub(rng::AbstractRNG, s::Solution) = f(rng, s; mode=:perturb)  