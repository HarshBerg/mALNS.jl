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
            # TODO: update t and h
        end
        # TODO: re-insert at best position
    end
    return s
end
# TODO: @code_warntype move!(rng, k, s; scope=:inter & :intra).

function swap!(rng::AbstractRNG, k::Int, s::Solution; scope::Symbol)
end

function opt!(rng::AbstractRNG, k::Int, s::Solution; scope::Symbol)
end