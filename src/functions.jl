
@inline isdepot(n::Node) = isone(n.i)

@inline iscustomer(n::Node) = !isone(n.i) 

@inline isopen(n::Node) = iszero(n.v)

@inline isclose(n::Node) = !iszero(n.v)

@inline Base.isequal(n::Node, m::Node) = isequal(n.i, m.i)

@inline isopt(v::Vehicle) = !iszero(v.n)

@inline Base.isequal(u::Vehicle, v::Vehicle) = isequal(u.i, v.i)

function vectorize(s::Solution)
    Z = Int[]
    G = s.G
    N = G.N
    V = G.V
    for v in V
        n = N[v.s]
        push!(Z, 1)  # depot index
        for _ in 1:v.n 
            push!(Z, n.i)  # customer node index
            n = N[n.h]
        end
        push!(Z, 1)  # depot index
    end
    return Z
end

@inline f(s::Solution) = s.c + s.p * 10 ^ ceil(log10(s.c))

@inline h(s::Solution) = hash(s) # TODO: Alternatively, check hashing with Node Vector instead

@inline Base.deepcopy_internal(G::Graph, dict::IdDict) = Graph(Base.deepcopy_internal(G.N, dict), G.A, Base.deepcopy_internal(G.V, dict))

@inline Base.deepcopy_internal(s::Solution, dict::IdDict) = Solution(Base.deepcopy_internal(s.G, dict), s.c, s.p)

"""
    relatedness(n₁::Node, n₂::Node)

Returns the relatedness between nodes `n₁` and `n₂`, calculated based on 
their spatial proximity.
"""
@inline function relatedness(n₁::Node, n₂::Node)
    ϵ = 1e-3
    d = abs(n₁.x - n₂.x) + abs(n₁.y - n₂.y)
    r = 1 / (d + ϵ)
    return r
end

"""
    relatedness(a₁::Arc, a₂::Arc)

Returns the relatedness between arcs `a₁` and `a₂`, calculated based on 
their spatial proximity. 
"""
@inline function relatedness(a₁::Arc, a₂::Arc)
    ϵ = 1e-3
    d = abs(a₁.x - a₂.x) + abs(a₁.y - a₂.y)
    r = 1 / (d + ϵ)
    return r
end

"""
    relatedness(v₁::Vehicle, v₂::Vehicle)

Returns the relatedness between vehicles `v₁` and `v₂`, calculated based on 
their spatial proximity. 
"""
@inline function relatedness(v₁::Vehicle, v₂::Vehicle)
    ϵ = 1e-3
    d = abs(v₁.x - v₂.x) + abs(v₁.y - v₂.y)
    r = 1 / (d + ϵ)
    return r
end