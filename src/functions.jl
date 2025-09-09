# TODO: Complete the following functions
isdepot(n::Node) = 
iscustomer(n::Node) = 

isopen(n::Node) = 
isclose(n::Node) = 

isopt(v::Vehicle) = 

function vectorize(s::Solution)
    Z = Int[]
    # TODO: add depot node index before start node index
    # TODO: loop over vehicles adding customer node indicies
    # TODO: add depot node index after end node index
    return Z
end

f(s::Solution) = s.c
h(s::Solution) = hash(s) # TODO: Alternatively, check hashing with Node Vector instead