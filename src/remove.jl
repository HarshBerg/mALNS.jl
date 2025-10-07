# Random Removal
function randomnode!(rng::AbstractRNG, k::Int, s::Solution)
    N = s.G.N
    W = [isdepot(n) ? 0 : 1 for n ∈ N]
    V = s.G.V
    for _ in 1:k
        i = sample(rng, 1:length(N), Weights(W))
        n = N[i]
        t = N[n.t]
        h = N[n.h]
        v = V[n.v]
        removenode!(n, t, h, v, s)
        W[i] = 0
    end
    return s
end

function randomarc!(rng::AbstractRNG, k::Int, s::Solution)

    return s
end

function randomsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function randomvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

# Related Removal
function relatednode!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function relatedarc!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function relatedsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function relatedvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

# Worst Removal
function worstnode!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function worstarc!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function worstsegment!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end

function worstvehicle!(rng::AbstractRNG, k::Int, s::Solution)
    # TODO
    return s
end