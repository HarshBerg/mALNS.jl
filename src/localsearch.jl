# move, swap, opt.
# inter & intra route.

function move!(rng::AbstractRNG, k::Int, s::Solution; scope::Symbol)
    G = s.G
    N = G.N
    V = G.V
    I = eachindex(N)
    Wₙ = [isdepot(n) ? 0 : 1 for n ∈ N]
    Wᵥ = [[isequal(n.v, v.i) ? isequal(scope, :intra) : isequal(scope, :inter) for v ∈ V] for n ∈ N]
    for _ ∈ 1:k
        i = sample(rng, I, Weights(Wₙ))
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        c = 0.
        p = (t.i, h.i, v.i)
        z = f(s)
        v = sample(rng, V, Weights(Wᵥ[n.i]))
        t = N[1]
        h = N[v.s]
        for _ ∈ 0:v.n
            insertnode!(n, t, h, v, s)
            z′ = f(s)
            Δ  = z′ - z
            if Δ < c c, p = Δ, (t.i, h.i, v.i) end
            removenode!(n, t, h, v, s)
            t = h
            h = N[t.h]
        end
        # re-insert at best position
        t = N[p[1]]
        h = N[p[2]]
        v = V[p[3]]
        insertnode!(n, t, h, v, s)
    end
    return s
end

function swap!(rng::AbstractRNG, k::Int, s::Solution; scope::Symbol)
    G = s.G
    N = G.N
    V = G.V
    I = eachindex(N)
    Wₙ = [isdepot(n) ? 0 : 1 for n ∈ N]
    Wᵥ = [[isequal(n.v, v.i) ? isequal(scope, :intra) : isequal(scope, :inter) for v ∈ V] for n ∈ N]
    for _ ∈ 1:k
        i = sample(rng, I, Weights(Wₙ))
        n = N[i]
        tₙ = N[n.t]
        hₙ = N[n.h]
        vₙ = V[n.v]
        removenode!(n, tₙ, hₙ, vₙ, s)
        c = 0.
        p = (tₙ.i, hₙ.i, vₙ.i)
        z = f(s)
        vₘ = sample(rng, V, Weights(Wᵥ[n.i]))
        m  = N[vₘ.s]
        tₘ = N[1]
        hₘ = N[m.h]
        for _ ∈ 0:vₘ.n
            # t-n-n₂-h
            if isequal(n, m.t)
                insertnode!(n, m, hₘ, vₘ, s)
                z′ = f(s)
                Δ  = z′ - z
                if Δ < c c, p = Δ, (tₘ.i, hₘ.i, vₘ.i) end
                removenode!(n, m, hₘ, vₘ, s)
            # t-n₂-n-h
            elseif isequal(m, n.t)
                insertnode!(n, tₘ, m, vₘ, s)
                z′ = f(s)
                Δ  = z′ - z
                if Δ < c c, p = Δ, (tₘ.i, hₘ.i, vₘ.i) end
                removenode!(n, tₘ, m, vₘ, s)
            else 
                removenode!(m, tₘ, hₘ, vₘ, s)
                insertnode!(n, tₘ, hₘ, vₘ, s)
                insertnode!(m, tₙ, hₙ, vₙ, s)
                z′ = f(s)
                Δ  = z′ - z
                if Δ < c c, p = Δ, (tₘ.i, hₘ.i, vₘ.i) end
                removenode!(n, tₘ, hₘ, vₘ, s)
                removenode!(m, tₙ, hₙ, vₙ, s)
                insertnode!(m, tₘ, hₘ, vₘ, s)
            end
            tₘ = m
            m = hₘ
            hₘ = N[m.h]
        end
        # re-insert at best position
        tₘ = N[p[1]]
        hₘ = N[p[2]]
        vₘ = V[p[3]]
        m  = N[hₘ.t]
        if isequal(n, m.t)
            insertnode!(n, m, hₘ, vₘ, s)
        elseif isequal(m, n.t)
            insertnode!(n, tₘ, m, vₘ, s)
        else 
            removenode!(m, tₘ, hₘ, vₘ, s)
            insertnode!(n, tₘ, hₘ, vₘ, s)
            insertnode!(m, tₙ, hₙ, vₙ, s)
        end
    end
    return s
end

function opt!(rng::AbstractRNG, k::Int, s::Solution; scope::Symbol)
    G = s.G
    N = G.N
    V = G.V
    A = G.A
    I = eachindex(N)
    Wₙ = [isdepot(n) ? 0 : 1 for n ∈ N]
    Wᵥ = [[isequal(n.v, v.i) ? isequal(scope, :intra) : isequal(scope, :inter) for v ∈ V] for n ∈ N]
    for _ ∈ 1:k
        i = sample(rng, I, Weights(Wₙ))
        t = N[i]
        n = N[t.h]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s) 
        c = 0.
        p = (t.i, h.i, v.i)
        z = f(s)
        v = sample(rng, V, Weights(Wᵥ[n.i]))
        m = N[v.s]
        tₘ = N[1]
        hₘ = N[n.h]
        for _ ∈ 0:v.n
            if isequal(n, m.h)
                insertnode!(n, tₘ, m, v, s)
            elseif isequal(n, m.t)
                insertnode!(n, m, hₘ, v, s)
            else
                removenode!(n, tₘ, hₘ, v, s)
                insertnode!(m, t, tₘ, v, s)
                insertnode!(n, t, hₘ, v, s)
                z′ = f(s)
                Δ  = z′ - z
                if Δ < c c, p = Δ, (tₘ.i, hₘ.i, v.i) end
                removenode!(m, t, tₘ, v, s)
                removenode!(n, t, hₘ, v, s)
                insertnode!(n, tₘ, hₘ, v, s)
                tₘ = m
                m = hₘ
                hₘ = N[m.h]
            end
        end
        # re-insert at best position