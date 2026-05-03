"""
    remove!([rng::AbstractRNG], k::Int, s::Solution, method::Function)

Returns solution removing `k` nodes from solution s using the given `method`.

Available methods include,
- Random Node Removal       : `randomnode!`
- Random Arc Removal        : `randomarc!`
- Random Segment Removal    : `randomsegment!`
- Random Vehicle Removal    : `randomvehicle!`
- Related Node Removal      : `relatednode!`
- Related Arc Removal       : `relatedarc!`
- Related Segment Removal   : `relatedsegment!`
- Related Vehicle Removal   : `relatedvehicle!`
- Worst Node Removal        : `worstnode!`
- Worst Arc Removal         : `worstarc!`
- Worst Segment Removal     : `worstsegment!`
- Worst Vehicle Removal     : `worstvehicle!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
remove!(rng::AbstractRNG, k::Int, s::Solution, method::Function)::Solution = method(rng, k, s)
remove!(k::Int, s::Solution, method::Function) = remove!(Random.GLOBAL_RNG, k, s, method)

"""
    randomnode!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing exactly `k` nodes 
selected randomly.
"""
function randomnode!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # node indices
    I = eachindex(N)
    # set node weights: uniform
    W = zeros(Int, I)
    for i ∈ I
        n = N[i]
        if isdepot(n) continue end
        W[i] = 1
    end
    # loop: remove exactly k sampled nodes
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
    end
    # return solution
    return s
end

"""
    relatednode!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing exactly `k` nodes 
selected based on relatedness to a pivot node.
"""
function relatednode!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    A = G.A
    V = G.V
    # node indices
    I = eachindex(N)
    # randomize a pivot node
    p = N[rand(rng, I)]
    xₚ, yₚ = p.x, p.y
    rₚ, θₚ = p.r, p.θ
    # set node weights: relatedness
    W = zeros(Float64, I)
    φ = rand(rng, Uniform(0., 1.))
    for i ∈ I
        n = N[i]
        if isdepot(n) continue end
        xₙ, yₙ = n.x, n.y
        rₙ, θₙ = n.r, n.θ
        d = φ * hypot(xₙ - xₚ, yₙ - yₚ) + (1 - φ) * abs(rₙ * θₙ - rₚ * θₚ) + 1e-3
        r = 1 / d
        W[i] = r
    end
    # loop: remove exactly k sampled nodes
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
    end
    # return solution
    return s
end

"""
    worstnode!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing exactly `k` nodes 
selected based on removal cost.
"""
function worstnode!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    A = G.A
    V = G.V
    # node indices
    I = eachindex(N)
    # set node weights: cost
    W = zeros(Int, I)
    for i ∈ I
        n = N[i]
        if isdepot(n) continue end
        t = N[n.t]
        h = N[n.h]
        W[i] = (A[t.i, n.i].c + A[n.i, h.i].c) - A[t.i, h.i].c
    end
    # loop: remove exactly k sampled nodes
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
    end
    # return solution
    return s
end

"""
    randomarc!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing at least `k` nodes 
from arcs selected randomly.
"""
function randomarc!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    A = G.A
    V = G.V
    # arc indices
    I = eachindex(A)
    # set arc weights: uniform
    W = zeros(Int, I)
    for i ∈ I
        a = A[i]
        t = N[a.t]
        h = N[a.h] 
        W[i] = isequal(t.h, h.i) ? 1 : 0
    end
    # loop: until at least k nodes are removed
    c = 0
    while c < k
        # sample an arc
        i = sample(rng, I, Weights(W))
        a = A[i]
        # remove tail node
        n = N[a.t]
        if iscustomer(n) && isclose(n)
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # remove head node
        n = N[a.h]
        if iscustomer(n) && isclose(n)
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # update arc weight
        W[i] = 0
    end
    # return solution
    return s
end

"""
    relatedarc!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing at least `k` nodes f
rom arcs selected based on relatedness to s pivot arc.
"""
function relatedarc!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    A = G.A
    V = G.V
    # arc indices
    I = eachindex(A)
    # randomize a pivot arc
    p = A[sample(rng, I, Weights([isequal(N[A[i].t].h, N[A[i].h].i) ? 1 : 0 for i ∈ I]))]
    xₚ, yₚ = (N[p.t].x + N[p.h].x) / 2, (N[p.t].y + N[p.h].y) / 2
    rₚ, θₚ = (N[p.t].r + N[p.h].r) / 2, (N[p.t].θ + N[p.h].θ) / 2
    # set arc weights: relatedness
    W = zeros(Float64, I)
    φ = rand(rng, Uniform(0., 1.))
    for i ∈ I
        a = A[i]
        t = N[a.t]
        h = N[a.h]
        xₐ, yₐ = (N[a.t].x + N[a.h].x) / 2, (N[a.t].y + N[a.h].y) / 2
        rₐ, θₐ = (N[a.t].r + N[a.h].r) / 2, (N[a.t].θ + N[a.h].θ) / 2
        d = φ * hypot(xₐ - xₚ, yₐ - yₚ) + (1 - φ) * abs(rₐ * θₐ - rₚ * θₚ) + 1e-3
        r = 1 / d
        W[i] = isequal(t.h, h.i) ? r : 0.
    end
    # loop: until at least k nodes are removed
    c = 0
    while c < k
        # sample an arc
        i = sample(rng, I, Weights(W))
        a = A[i]
        # remove tail node
        n = N[a.t]
        if iscustomer(n) && isclose(n)
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # remove head node
        n = N[a.h]
        if iscustomer(n) && isclose(n)
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # update arc weight
        W[i] = 0
    end
    # return solution
    return s
end

"""
    worstarc!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing at least `k` nodes 
from arcs selected based on removal cost.
"""
function worstarc!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    A = G.A
    V = G.V
    # arc indices
    C = CartesianIndices(A)
    I = eachindex(A)
    # set arc weights: cost
    W = zeros(Int, I)
    for i ∈ I
        a = A[i]
        t = N[a.t]
        h = N[a.h] 
        W[i] = isequal(t.h, h.i) ? a.c : 0
    end
    # loop: until at least k nodes are removed
    c = 0
    while c < k
        # sample an arc
        i = sample(rng, I, Weights(W))
        a = A[i]
        # remove tail node
        n = N[a.t]
        if iscustomer(n) && isclose(n)
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # remove head node
        n = N[a.h]
        if iscustomer(n) && isclose(n)
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # update arc weight
        W[i] = 0
    end
    # return solution
    return s
end

"""
    randomsegment!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution 's' after removing at least 'k' nodes 
from segments selected randomly.
"""
function randomsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    A = G.A
    V = G.V
    # node indices
    I = eachindex(N)
    # set node weights: uniform
    W = zeros(Int, I)
    for i ∈ I
        n = N[i]
        if isdepot(n) continue end
        W[i] = 1
    end
    # loop: until at least k nodes are removed
    c = 0 
    while c < k
        i = sample(rng, I, Weights(W))
        n = N[i]
        v = V[n.v]
        # determine segment size and divide it into two halves around node n
        cₒ = sample(rng, 1:v.n)
        cₕ = cₒ ÷ 2
        cₜ = cₒ - cₕ - 1
        # remove cₕ nodes after n
        n = N[i]
        h = N[n.h]
        for _ in 1:cₕ
            n = h
            t = isdepot(n) ? N[v.e] : N[n.t]
            h = isdepot(n) ? N[v.s] : N[n.h]
            if isdepot(n) continue end
            removenode!(n, t, h, v, s)
            W[n.i] = 0
            c += 1
        end
        # remove cₜ nodes before n
        n = N[i]
        t = N[n.t]
        for _ in 1:cₜ
            n = t
            t = isdepot(n) ? N[v.e] : N[n.t]
            h = isdepot(n) ? N[v.s] : N[n.h]
            if isdepot(n) continue end
            removenode!(n, t, h, v, s)
            W[n.i] = 0
            c += 1
        end
        # remove node n
        n = N[i]
        h = N[n.h]
        t = N[n.t]
        removenode!(n, t, h, v, s)
        W[n.i] = 0
        c += 1
    end
    # return solution
    return s
end

"""
    relatedsegment!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing exactly `k` nodes 
from segments selected based on relatedness to a pivot segment.
"""
function relatedsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    A = G.A
    V = G.V
    # node indices
    I = eachindex(N)
    # randomize a pivot node
    p = N[rand(rng, I)]
    xₚ, yₚ = p.x, p.y
    rₚ, θₚ = p.r, p.θ
    # set node weights: relatedness
    W = zeros(Float64, I)
    φ = rand(rng, Uniform(0., 1.))
    for i ∈ I
        n = N[i]
        if isdepot(n) continue end
        xₙ, yₙ = n.x, n.y
        rₙ, θₙ = n.r, n.θ
        d = φ * hypot(xₙ - xₚ, yₙ - yₚ) + (1 - φ) * abs(rₙ * θₙ - rₚ * θₚ) + 1e-3
        r = 1 / d
        W[i] = r
    end
    # loop: until at least k nodes are removed
    c = 0 
    while c < k
        i = sample(rng, I, Weights(W))
        n = N[i]
        v = V[n.v]
        # determine segment size and divide it into two halves around node n
        cₒ = sample(rng, 1:v.n)
        cₕ = cₒ ÷ 2
        cₜ = cₒ - cₕ - 1
        # remove cₕ nodes after n
        n = N[i]
        h = N[n.h]
        for _ in 1:cₕ
            n = h
            t = isdepot(n) ? N[v.e] : N[n.t]
            h = isdepot(n) ? N[v.s] : N[n.h]
            if isdepot(n) continue end
            removenode!(n, t, h, v, s)
            W[n.i] = 0
            c += 1
        end
        # remove cₜ nodes before n
        n = N[i]
        t = N[n.t]
        for _ in 1:cₜ
            n = t
            t = isdepot(n) ? N[v.e] : N[n.t]
            h = isdepot(n) ? N[v.s] : N[n.h]
            if isdepot(n) continue end
            removenode!(n, t, h, v, s)
            W[n.i] = 0
            c += 1
        end
        # remove node n
        n = N[i]
        h = N[n.h]
        t = N[n.t]
        removenode!(n, t, h, v, s)
        W[n.i] = 0
        c += 1
    end
    # return solution
    return s
end

"""
    worstsegment!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution 's' after removing atleast 'k' nodes 
from worst segments.
"""
function worstsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    A = G.A
    V = G.V
    # node indices
    I = eachindex(N)
    # set node weights: cost
    W = zeros(Int, I)
    for i ∈ I
        n = N[i]
        if isdepot(n) continue end
        t = N[n.t]
        h = N[n.h]
        W[i] = (A[t.i, n.i].c + A[n.i, h.i].c) - A[t.i, h.i].c
    end
    # loop: until at least k nodes are removed
    c = 0 
    while c < k
        i = sample(rng, I, Weights(W))
        n = N[i]
        v = V[n.v]
        # determine segment size and divide it into two halves around node n
        cₒ = sample(rng, 1:v.n)
        cₕ = cₒ ÷ 2
        cₜ = cₒ - cₕ - 1
        # remove cₕ nodes after n
        n = N[i]
        h = N[n.h]
        for _ in 1:cₕ
            n = h
            t = isdepot(n) ? N[v.e] : N[n.t]
            h = isdepot(n) ? N[v.s] : N[n.h]
            if isdepot(n) continue end
            removenode!(n, t, h, v, s)
            W[n.i] = 0
            c += 1
        end
        # remove cₜ nodes before n
        n = N[i]
        t = N[n.t]
        for _ in 1:cₜ
            n = t
            t = isdepot(n) ? N[v.e] : N[n.t]
            h = isdepot(n) ? N[v.s] : N[n.h]
            if isdepot(n) continue end
            removenode!(n, t, h, v, s)
            W[n.i] = 0
            c += 1
        end
        # remove node n
        n = N[i]
        h = N[n.h]
        t = N[n.t]
        removenode!(n, t, h, v, s)
        W[n.i] = 0
        c += 1
    end
    # return solution
    return s
end

"""
    randomvehicle!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing at least `k` nodes 
from vehicles selected randomly.
"""
function randomvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # vehicle indices
    I = eachindex(V)
    # set vehicle weights: uniform
    W = ones(Int, I)
    # loop: until at least k nodes are removed
    c = 0
    while c < k
        # sample a vehicle
        i = sample(rng, I, Weights(W))
        v = V[i]
        # remove all associated nodes
        for _ ∈ 1:v.n
            n = N[v.s]
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # update vehicle weight
        W[i] = 0
    end
    # return solution
    return s
end

"""
    relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing at least `k` nodes 
from vehicles selected based on relatedness to a pivot vehicle.
"""
function relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # vehicle indices
    I = eachindex(V)
    # randomize a pivot vehicle
    p = V[rand(rng, I)]
    xₚ, yₚ = p.x, p.y
    rₚ, θₚ = hypot(p.x, p.y), atan(p.y / p.x)
    # set vehicle weights: relatedness
    W = zeros(Float64, I)
    φ = rand(rng, Uniform(0., 1.))
    for i ∈ I
        v = V[i]
        xₙ, yₙ = v.x, v.y
        rₙ, θₙ = hypot(v.x, v.y), atan(v.y / v.x)
        d = φ * hypot(xₙ - xₚ, yₙ - yₚ) + (1 - φ) * abs(rₙ * θₙ - rₚ * θₚ) + 1e-3
        r = 1 / d
        W[i] = r
    end
    # loop: until at least k nodes are removed
    c = 0
    while c < k
        # sample a vehicle
        i = sample(rng, I, Weights(W))
        v = V[i]
        # remove all associated nodes
        for _ ∈ 1:v.n
            n = N[v.s]
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # update vehicle weight
        W[i] = 0
    end
    # return solution
    return s
end

"""
    worstvehicle!(rng::AbstractRNG, k::Int, s::Solution)

Returns solution `s` after removing at least `k` nodes 
from vehicles selected based on utilization.
"""
function worstvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    # initialiaze
    G = s.G
    N = G.N
    V = G.V
    # vehicle indices
    I = eachindex(V)
    # set vehicle weights: utilization
    W = zeros(Float64, I)
    for i ∈ I
        v = V[i]
        W[i] = 1 - v.l / v.q + 1e-3
    end
    # loop: until at least k nodes are removed
    c = 0
    while c < k
        # sample a vehicle
        i = sample(rng, I, Weights(W))
        v = V[i]
        # remove all associated nodes
        for _ ∈ 1:v.n
            n = N[v.s]
            t = N[n.t]
            h = N[n.h]
            v = V[n.v]
            removenode!(n, t, h, v, s)
            c += 1
        end
        # update vehicle weight
        W[i] = 0
    end
    # return solution
    return s
end