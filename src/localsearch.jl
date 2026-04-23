# move, swap, opt.
# inter & intra route.
"""
    localsearch!(rng::AbstractRNG, s::Solution; method::Symbol)

Apply a local search method to the solution `s` using random number generator `rng`. The local search method is specified by `method`.
"""
localsearch!(rng::AbstractRNG, k::Int, s::Solution; method::Function)::Solution = method(rng, k, s)

function move!(rng::AbstractRNG, k::Int, s::Solution; scope::Symbol)
    G = s.G
    N = G.N
    V = G.V
    I = eachindex(N)
    Wₙ = [isdepot(n) || isopen(n) ? 0 : 1 for n ∈ N]
    Wᵥ = [[isequal(n.v, v.i) ? isequal(scope, :intra) : (isequal(scope, :inter) && isopt(v)) for v ∈ V] for n ∈ N]
    if iszero(sum(Wₙ)) return s end
    for _ ∈ 1:k
        i = sample(rng, I, Weights(Wₙ))
        n = N[i]
        if iszero(sum(Wᵥ[n.i])) continue end
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

intramove!(rng::AbstractRNG, k::Int, s::Solution) = move!(rng, k, s; scope=:intra)
intermove!(rng::AbstractRNG, k::Int, s::Solution) = move!(rng, k, s; scope=:inter)

function swap!(rng::AbstractRNG, k::Int, s::Solution; scope::Symbol)
    G = s.G
    N = G.N
    V = G.V
    I = eachindex(N)
    Wₙ = [isdepot(n) || isopen(n) ? 0 : 1 for n ∈ N]
    Wᵥ = [[isequal(n.v, v.i) ? (isequal(scope, :intra) && v.n ≥ 2) : (isequal(scope, :inter) && isopt(v)) for v ∈ V] for n ∈ N]
    Wₘ = 
    if iszero(sum(Wₙ)) return s end
    for _ ∈ 1:k
        i = sample(rng, I, Weights(Wₙ))
        n = N[i]
        if iszero(sum(Wᵥ[n.i])) continue end
        tₙ = N[n.t]
        hₙ = N[n.h]
        vₙ = V[n.v]
        c = 0.
        p = (tₙ.i, hₙ.i, vₙ.i)
        z = f(s)
        vₘ = sample(rng, V, Weights(Wᵥ[n.i]))
        m  = N[vₘ.s]
        tₘ = N[1]
        hₘ = N[m.h]
        for _ ∈ 0:vₘ.n
            if isequal(n, m) continue 
            # t-n-m-h
            elseif isequal(n, m.t)
                removenode!(n, tₙ, hₙ, vₙ, s)
                insertnode!(n, m, hₘ, vₘ, s)
                z′ = f(s)
                Δ  = z′ - z
                if Δ < c c, p = Δ, (tₘ.i, hₘ.i, vₘ.i) end
                removenode!(n, m, hₘ, vₘ, s)
                insertnode!(n, tₙ, hₙ, vₙ, s)
            # t-n₂-n-h
            elseif isequal(m, n.t)
                removenode!(n, tₙ, hₙ, vₙ, s)
                insertnode!(n, tₘ, m, vₘ, s)
                z′ = f(s)
                Δ  = z′ - z
                if Δ < c c, p = Δ, (tₘ.i, hₘ.i, vₘ.i) end
                removenode!(n, tₘ, m, vₘ, s)
                insertnode!(n, tₙ, hₙ, vₙ, s)
            else 
                removenode!(n, tₙ, hₙ, vₙ, s)
                removenode!(m, tₘ, hₘ, vₘ, s)
                insertnode!(n, tₘ, hₘ, vₘ, s)
                insertnode!(m, tₙ, hₙ, vₙ, s)
                z′ = f(s)
                Δ  = z′ - z
                if Δ < c c, p = Δ, (tₘ.i, hₘ.i, vₘ.i) end
                removenode!(n, tₘ, hₘ, vₘ, s)
                removenode!(m, tₙ, hₙ, vₙ, s)
                insertnode!(n, tₙ, hₙ, vₙ, s)
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

interswap!(rng::AbstractRNG, k::Int, s::Solution) = swap!(rng, k, s; scope=:inter)
intraswap!(rng::AbstractRNG, k::Int, s::Solution) = swap!(rng, k, s; scope=:intra)

"""
function intraopt!(rng::AbstractRNG, k::Int, s::Solution)
    G = s.G
    N = G.N
    V = G.V
    # initialize
    Wₙ = [isdepot(n) ? 0 : 1 for n ∈ N]
    Wₘ = [[isequal(n.v, m.v) ? (isequal(n, m) ? 0 : 1) : 0 for m ∈ N] for n ∈ N] 
    for _ ∈ 1:k
        n = sample(rng, N, Weights(Wₙ))
        m = sample(rng, N, Weights(Wₘ[n.i]))
        p = N[m.h]
        while !isdepot(p)
            if isequal(n, p)
                removenode!(n, N[n.t], N[n.h], V[n.v], s)
                insertnode!(n, m.t, m, V[m.v], s)
                break
            end
            p = N[p.h]
        end
        if isequal(n, m) continue end
        v = V[n.v]
        removenode!(n, N[n.t], N[n.h], v, s)
"""