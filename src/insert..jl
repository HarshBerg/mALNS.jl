function best!(rng::AbstractRNG, s::Solution; mode::Symbol)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    A = G.A
    # intialize
    L = [n for n ∈ N if iscustomer(n) && isopen(n)] # set of open nodes
    I = eachindex(L)
    J = eachindex(V)
    C = fill(Inf, (I,J))                            # C[i,j]: best insertion cost of node L[i] in vehicle route V[j]
    P = fill((0,0), (I,J))                          # P[i,j]: best insertion position on node L[i] in vehicle route V[j]
    for _ ∈ I
        i = rand(rng, I)
        p = L[i]
        for j ∈ J
            v = V[j]
            n = v.s
            h = n.h 
            for _ ∈ 0:v.n
                n = h
                t = isdepot(n) ? N[v.e] : N[n.t]
                h = isdepot(n) ? N[v.s] : N[n.h]
                insertnode!(p, t, h, v, s)
                Δ = (A[t.i, p.i].c + A[p.i, h.i].c - A[t.i, h.i].c)
                if Δ < C[i, j]
                    C[i, j] = Δ
                    P[i, j] = (t.i, h.i)
                end
                removenode!(p, t, h, v, s)
            end
        end
        j = argmin(C[i, :])
        if isfinite(C[i, j])
            p = L[i]
            tᵢ, hᵢ = P[i, j]
            t = N[tᵢ]
            h = N[hᵢ]
            v = V[j]
            insertnode!(p, t, h, v, s)
        end
    end
    return s
end

#fprecise(rng::AbstractRNG, s::Solution) = f(rng, s; mode=:precise)
#fpertrub(rng::AbstractRNG, s::Solution) = f(rng, s; mode=:perturb)  