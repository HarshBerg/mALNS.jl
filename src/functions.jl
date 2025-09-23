
@inline isdepot(n::Node) = isone(n.i)

@inline iscustomer(n::Node) = !isone(n.i) 

@inline isopen(n::Node) = iszero(n.v)

@inline isclose(n::Node) = !iszero(n.v)

@inline isopt(v::Vehicle) = !iszero(v.n)

function vectorize(s::Solution)
    Z = Int[]
    # TODO: add depot node index before start node index
    # TODO: loop over vehicles adding customer node indicies
    # TODO: add depot node index after end node index
    return Z
end

@inline f(s::Solution) = s.c + s.p * 10 ^ ceil(log10(s.c))

@inline h(s::Solution) = hash(s) # TODO: Alternatively, check hashing with Node Vector instead

@inline Base.deepcopy_internal(G::Graph, dict::IdDict) = Graph(Base.deepcopy_internal(G.N, dict), G.A, Base.deepcopy_internal(G.V, dict))

@inline Base.deepcopy_internal(s::Solution, dict::IdDict) = Solution(Base.deepcopy_internal(s.G, dict), s.c, s.p)