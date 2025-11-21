# best
# greedy
# regret

function f(rng::AbstractRNG, s::Solution; mode::Symbol)
    # pre-initialize
    G = s.G
    N = G.N
    V = G.V
    # intialize
    L = [n for n ∈ N if iscustomer(n) && isopen(n)] # set of open nodes
    I = eachindex(L)
    J = eachindex(V)
    C = fill(Inf, (I,J))                            # C[i,j]: best insertion cost of node L[i] in vehicle route V[j]
    P = fill((0,0), (I,J))                          # P[i,j]: best insertion position on node L[i] in vehicle route V[j]
    for _ ∈ 1:I 
    
    end
    return s
end
fprecise(rng::AbstractRNG, s::Solution) = f(rng, s; mode=:precise)
fpertrub(rng::AbstractRNG, s::Solution) = f(rng, s; mode=:perturb)