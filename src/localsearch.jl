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
# TODO: @code_warntype move!(rng, k, s; scope=:inter & :intra).

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
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        c = 0.
        p = (t.i, h.i, v.i)
        z = f(s)
        v₂ = sample(rng, V, Weights(Wᵥ[n.i]))
        n₂ = N[v₂.s]
        t₂ = N[1]
        h₂ = N[n₂.h]
        for _ ∈ 0:v₂.n
            removenode!(n₂, t₂, h₂, v₂, s)
            insertnode!(n, t₂, h₂, v₂, s)
            insertnode!(n₂, t, h, v, s)
            z′ = f(s)
            Δ  = z′ - z
            if Δ < c c, p = Δ, (t₂.i, h₂.i, v₂.i) end
            removenode!(n, t₂, h₂, v₂, s)
            removenode!(n₂, t, h, v, s)
            insertnode!(n₂, t₂, h₂, v₂, s)
            t₂ = n₂
            n₂ = h₂
            h₂ = N[n₂.h]
        end
        # re-insert at best position
        t₂ = N[p[1]]
        h₂ = N[p[2]]
        v₂ = V[p[3]]
        n₂ = N[h₂.t]
        removenode!(n₂, t₂, h₂, v₂, s)
        insertnode!(n, t₂, h₂, v₂, s)
        insertnode!(n₂, t, h, v, s)
    end
    return s
end

function opt!(rng::AbstractRNG, k::Int, s::Solution; scope::Symbol)
end