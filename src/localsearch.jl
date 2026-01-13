# move, swap, opt.
# inter & intra route.

function intramove!(rng::AbstractRNG, s::Solution, k::Int)
    G = s.G
    N = G.N
    V = G.V
    I = eachindex(N)
    W = zeros(Int, I)
    C = Inf
    for i ∈ I
        n = N[i]
        if isdepot(n) continue end
        W[i] = 1
    end
    c = 0
    while c < k
        i = sample(rng, I, Weights(W))
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        W[i] = 0
        c += 1
        # find best insertion
        z = f(s)
        v = V[n.v]
        t = N[1]
        h = N[v.s]
        for _ ∈ 0:v.n
            insertnode!(n, t, h, v, s)
            z′ = f(s)
            Δ  = z′ - z
            if Δ < 0 
end
