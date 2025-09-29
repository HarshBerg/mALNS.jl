function build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    # read instance file
    df = CSV.read("$dir/$instance.vrp", DataFrame, silencewarnings=true)
    # fetch key indices
    k₁ = findfirst(contains("DIMENSION"), df[:,1])
    k₂ = findfirst(contains("CAPACITY"), df[:,1])
    k₃ = findfirst(contains("NODE_COORD_SECTION"), df[:,1])
    k₄ = findfirst(contains("DEMAND_SECTION"), df[:,1])
    # nodes
    n = parse(Int, df[k₁,2])
    N = Vector{Node}(undef, n)
    for i ∈ 1:n
        x = parse(Int, split(df[k₃+i,1])[2])
        y = parse(Int, split(df[k₃+i,1])[3])
        q = parse(Int, split(df[k₄+i,1])[2])
        N[i] = Node(i, x, y, q)
    end
    # arcs
    A = Matrix{Arc}(undef, n, n)
    for i ∈ 1:n
        xᵢ = N[i].x
        yᵢ = N[i].y 
        for j in 1:n
            xⱼ = N[j].x
            yⱼ = N[j].y 
            c = ((xⱼ-xᵢ)^2 + (yⱼ-yᵢ)^2)^0.5
            A[i,j] = Arc(i, j, c)
        end
    end
    # vehicles
    q = parse(Int, df[k₂,2])
    V = Vector{Vehicle}(undef, n-1)
    for i ∈ 1:(n-1) V[i] = Vehicle(i, q) end
    # create graph
    G = Graph(N, A, V)
    return G
end

"""

"""
function initialize(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    G = build(instance; dir=dir)
    s = Solution(G)
    G = s.G
    N = G.N
    A = G.A
    V = G.V
    K = eachindex(N)
    # initialize one-to-one routes
    d = N[1]
    for k ∈ K
        n = N[k]
        if isdepot(n) continue end
        v = V[k-1]
        insertnode!(n, d, d, v, s)
    end
    # calculate savings for each pair of customers
    Δ = zeros(Float64, (K,K))
    for i ∈ K
        if isone(i) continue end
        for j ∈ K
            if isone(j) continue end
            if isequal(i,j) continue end
            if j > i continue end
            δ = (A[i,1].c + A[1,j].c) - A[i,j].c
            Δ[i,j] = δ
        end
    end
    # perform feasible greedy merge
    P = sortperm(vec(Δ), rev=true)
    T = Tuple.(CartesianIndices(Δ)[P])
    for (i,j) ∈ T
        if iszero(Δ[i,j]) break end
        nᵢ = N[i]   
        nⱼ = N[j]
        if nᵢ.h ≠ 1 || nⱼ.t ≠ 1 continue end # Node i must be the last in its route and Node j must be the first in its route
        vᵢ = V[nᵢ.v]
        vⱼ = V[nⱼ.v]

        if vᵢ.i == vⱼ.i continue end # Vehicles must be different

        if vᵢ.l + vⱼ.l > vᵢ.q continue end # Capacity constraint

        removenode!(nᵢ, N[nᵢ.t], d, vᵢ, s)
        insertnode!(nᵢ, N[nᵢ.t], nᵢ, vᵢ, s)  #error here
        removenode!(nⱼ, d, N[nⱼ.h], vⱼ, s)
        insertnode!(nⱼ, nᵢ, N[nᵢ.h], vᵢ, s)
    end
    return s

end