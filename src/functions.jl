
isdepot(n::Node) = isone(n.i)

iscustomer(n::Node) = !isone(n.i) 

isopen(n::Node) = iszero(n.v)

isclose(n::Node) = !iszero(n.v)

isopt(v::Vehicle) = !iszero(v.n)

function vectorize(s::Solution)
    Z = Int[]
    # TODO: add depot node index before start node index
    # TODO: loop over vehicles adding customer node indicies
    # TODO: add depot node index after end node index
    return Z
end

f(s::Solution) = s.c + s.p * 10 ^ ceil(log10(s.c))
h(s::Solution) = hash(s) # TODO: Alternatively, check hashing with Node Vector instead